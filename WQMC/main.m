% Water Quality Simulation
%%
clear; close('all'); clc;
%% Load EPANET MATLAB TOOLKIT
start_toolkit;
%% run EPANET MATLAB TOOLKIT to obtain data
Network = 3; % Don't use case 2
switch Network
    case 1
        % Quality Timestep = 1 min, and  Global Bulk = -0.3, Global Wall= -0.0
        NetworkName = 'Threenode-cl-2.inp';
    case 2
        % Don't not use one: Quality Timestep = 5 min, and  Global Bulk = -0.3, Global Wall=
        % -1.0
        NetworkName = 'tutorial8node.inp';
    case 3
        % Quality Timestep = 1 min, and  Global Bulk = -0.3, Global Wall= -0.0
        NetworkName = 'tutorial8node1.inp';
    case 4
        % Quality Timestep = 1 min, and  Global Bulk = -0.3, Global Wall=
        % -0.0; initial value: J2 = 0.5 mg/L, J6 = 1.2 mg/L, R1 = 0.8 mg/L;
        % segment = 1000;
        NetworkName = 'tutorial8node1inital.inp';
    case 5
        % Quality Timestep = 1 min, and  Global Bulk = -0.5, Global Wall=
        % -0.0; 
        NetworkName = 'Net1-1min.inp';
    case 6
        % The initial value is slightly different
        NetworkName = 'Net1-1mininitial.inp';
    otherwise
        disp('other value')
end

switch Network
    case 1
        SavedData = 'PreparedData3node.mat';
    case {2,3}
        SavedData = 'PreparedData8node.mat';
    case 4
        SavedData = 'PreparedData8nodeInitial.mat';
    case 5
        SavedData = 'PreparedDataNet1.mat';
    case 6
        SavedData = 'PreparedDataNet1Initial.mat';
    otherwise
        disp('other value')
end

if exist(SavedData, 'file') == 2
    disp('file exists, loading...')
    load('PreparedDataNet1.mat')
else
    if ispc %Code to run on Windows platform
        disp('file does not exist, simulating...')
        PrepareData
    end
end

% if ismac
%     % Code to run on Mac platform
% elseif isunix
%     % Code to run on Linux platform
% elseif ispc
%     % Code to run on Windows platform
% else
%     disp('Platform not supported')
% end
%% Initialization

% initialize concentration at nodes
x0 = zeros(NumberofX,1);
C0 = [QsN(1,:) QsL(1,:)];

Head0 = Head(1,:);
x0 = InitialConcentration(x0,C0,MassEnergyMatrix,Head0,IndexInVar,ElementCount);

%% Start iteration
CurrentTime = 0; % change this later;
kC = 0;
CurrentNodeTankVolume = NodeTankVolume(1,:);

SimutionTimeInMinute = Constants4Concentration.SimutionTimeInMinute;
SetTimeParameter = SimutionTimeInMinute*Constants4Concentration.MinInSecond/double(TimeQualityStep);
maxiumStorage = SimutionTimeInMinute*Constants4Concentration.MinInSecond/0.1; % assume delta_t is 0.1, in fact, much less than it
%X = zeros(NumberofX,maxiumStorage);
X = sparse(NumberofX,maxiumStorage);
TimeRecord = [];
X(:,1) = x0;
iteration = 1;

%% INITIALIZE BOOSTER
% flow of Booster, assume we put booster at each nodes, so the size of it
% should be the number of nodes.

nodeCount = double(JunctionCount+ReservoirCount+TankCount);
Location_B = {'22'}; % NodeID here;
flowRate_B = [200]; % unit: GPM
C_B = [1]; % unit: mg/L % Concentration of booster
NodeID = Variable_Symbol_Table(1:nodeCount,1);
[q_B,C_B] = InitialBoosterFlow(nodeCount,Location_B,flowRate_B,NodeID,C_B);

tic
while(kC <= SetTimeParameter)
%profile on
TimeRecord = [TimeRecord CurrentTime];
CurrentTime
kC = floor(CurrentTime/double(TimeQualityStep))+1;
CurrentVelocityPipe = VelocityPipe(kC,:);
CurrentNodeNetFlowTank = NodeNetFlowTank(kC,:);
CurrentFlow = Flow(kC,:);
CurrentFlow = CurrentFlow';

delta_t = LinkLengthPipe./NumberofSegment./CurrentVelocityPipe;
delta_t = min(delta_t); % the minium step length for all pipes

% for junctions
[A_J, B_J] = ConstructMatrixForJunction(CurrentFlow,JunctionMassMatrix,ElementCount,IndexInVar,q_B);
% for reservoirs
A_R = ConstructMatrixForReservoir(ElementCount,IndexInVar);
% for tanks                  
[A_TK,B_TK,CurrentNodeTankVolume] = ConstructMatrixForTank(delta_t,CurrentFlow,CurrentNodeTankVolume,TankMassMatrix,ElementCount,IndexInVar,aux,q_B);
% for Pipes
EnergyMatrixPipe= MassEnergyMatrix(PipeIndex,:);
A_P = ConstructMatrixForPipe(delta_t,CurrentFlow,EnergyMatrixPipe,ElementCount,IndexInVar,aux);
% for Pumps
EnergyMatrixPump= MassEnergyMatrix(PumpIndex,:);
A_M = ConstructMatrixForPump(EnergyMatrixPump,ElementCount,IndexInVar);
% for Valves
EnergyMatrixValve= MassEnergyMatrix(ValveIndex,:);
A_W = ConstructMatrixForValve(EnergyMatrixValve,ElementCount,IndexInVar);
% construct A;
A = [A_J;A_R;A_TK;A_P;A_M;A_W];
% construct B;
B_R = sparse(ReservoirCount,nodeCount);
B_P = sparse(NumberofSegment*PipeCount,nodeCount);
B_M = sparse(PumpCount,nodeCount);
B_W = sparse(ValveCount,nodeCount);
B = [B_J;B_R;B_TK;B_P;B_M;B_W];

%[A,NextNodeTankVolume] = ConstructMatrix(delta_t,CurrentVelocityPipe,CurrentNodeTankVolume,CurrentNodeNetFlowTank,aux,IndexInVar,ElementCount);
x0 = A*x0 + B*C_B;
iteration = iteration + 1;
X(:,iteration) = x0;
CurrentTime = CurrentTime + delta_t;
%profile viewer
end

toc

disp('Simulation Done!! Start to organize data')
TimeRecord_Min = TimeRecord/double(TimeQualityStep);
TimeRecord_Min = floor(TimeRecord_Min);
% find min level data;
maxSimulationTime = max(TimeRecord_Min);
X_Min = [];
X_Min = [X_Min X(:,1)];
for i = 1:maxSimulationTime
    IndexofEachMin = min(find(TimeRecord_Min==i));
    X_Min = [X_Min X(:,IndexofEachMin)];
end

% find average data;
[m,n] = size(X_Min);
X_Min_Average = zeros(NumberofElement,n);
basePipeCIndex = min(Pipe_CIndex);
First = basePipeCIndex:basePipeCIndex+NumberofSegment-1;
for i = 1:n
    X_Min_Average(JunctionIndexInOrder,i) = X_Min(Junction_CIndex,i);
    X_Min_Average(ReservoirIndexInOrder,i) = X_Min(Reservoir_CIndex,i);
    X_Min_Average(TankIndexInOrder,i) = X_Min(Tank_CIndex,i);
    for j = 1:PipeCount
        Indexrange = (j-1)*NumberofSegment + First;
        X_Min_Average(PipeIndexInOrder(j),i) = mean(X_Min(Indexrange,i));
    end
    X_Min_Average(PumpIndexInOrder,i) = X_Min(Pump_CIndex-PipeCount*NumberofSegment,i);
    X_Min_Average(ValveIndexInOrder,i) = X_Min(Valve_CIndex-PipeCount*NumberofSegment,i);
end

% EPANET's result;
Qs = [QsN QsL]';

disp('Data organized')

