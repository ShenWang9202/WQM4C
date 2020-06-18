% Main Program to do MPC control considering Unknown, Demand, Parameter
% Uncertainties

% System: only on Windows
% Author: Shen Wang
% Date: 3/7/2020

% In order to do LDE without input u, just (1) uncomment x_estimated = A *
% x_estimated + B * U;  and make it as x_estimated = A * x_estimated; in
% EstimateState_XX.m source file.
% (2) don't apply control action to epanet (3) remove all uncertainty


clear all
clc
close all
%% Load EPANET MATLAB TOOLKIT
start_toolkit;

% check this example Toolkit_EX3_Minimum_chlorine_residual.m
%% run EPANET MATLAB TOOLKIT to obtain data
Kb_uncertainty = 0.1;
Kw_uncertainty = 0.1;

COMPARE = 1;

Network = 1; % Don't use case 2
Network = 4; % Don't use case 2
switch Network
    case 1
        % Quality Timestep = 1 min, and  Global Bulk = -0.3, Global Wall= -0.0
        % NetworkName = 'Threenode-cl-2-paper.inp'; pipe flow direction
        % never changed, and the result is perfectly matched with EPANET
        %NetworkName = 'Threenode-cl-3-paper.inp'; % Pipe flow direction changes
        NetworkName = 'Threenode-cl-2-paper.inp'; % topogy changes
        Unknown_Happen_Time = 200;
        PipeID_Cell = {'P1'};
        JunctionID_Cell = {'J2'};
        Sudden_Concertration = 1.0; % Suddenly the concentration jumps to this value for no reason
        filename = 'Three-node_1day.mat';
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
        %         NetworkName = 'tutorial8node1inital2.inp';
        filename = '8node_1day.mat';
    case 5
        % Quality Timestep = 1 min, and  Global Bulk = -0.5, Global Wall=
        % -0.0;
        NetworkName = 'Net1-1min.inp';
    case 6
        % The initial value is slightly different
        NetworkName = 'Net1-1mininitial.inp';
    case 7
        % Quality Timestep = 1 min, and  Global Bulk = -0.3, Global Wall= -0.0
        NetworkName = 'Net1-1min-new-demand-pattern.inp';
        Unknown_Happen_Time = 3000; % Unknow Disturbance happened at the 3000-th minutes
        PipeID_Cell = {'P11','P21','P31'};
        JunctionID_Cell = {'J11','J21','J31'};
        Sudden_Concertration = 0.5;
        filename = 'Net1_4days.mat';
    case 8
        % Quality Timestep = 1 min, and  Global Bulk = -0.3, Global Wall= -0.0
        NetworkName = 'Fournode-Cl-As-1.inp';
    case 9
        NetworkName = 'Net3-NH2CL-24hour.inp';
    otherwise
        disp('other value')
end

%% Prepare constants data for MPC
PrepareData4Control

%% initialize concentration at nodes

nx = NumberofX; % Number of states

% initialize BOOSTER
% flow of Booster, assume we put booster at each nodes, so the size of it
% should be the number of nodes.
JunctionCount = double(JunctionCount);
ReservoirCount = double(ReservoirCount);
TankCount = double(TankCount);
nodeCount = JunctionCount + ReservoirCount + TankCount;

switch Network
    case 1
        Location_B = {'J2'}; % NodeID here;
        flowRate_B = [10]; % unit: GPM
        Price_B = [1];
    case {2,3,4}
        Location_B = {'J3','J7'}; % NodeID here;
        flowRate_B = [10,10]; % unit: GPM
        Price_B = [1,1];
    case {5,6,7}
        Location_B = {'J11','J22','J31'}; % NodeID here;
        flowRate_B = [10,10,10]; % unit: GPM
        Price_B = [1,1,1];
    case 8
        Location_B = {'J2'}; % NodeID here;
        flowRate_B = [100]; % unit: GPM
        Price_B = [1];
    case 9
        Location_B = {'J10'}; % NodeID here;
        flowRate_B = [100]; % unit: GPM
        Price_B = [1];
    otherwise
        disp('other value')
end
NodeID = Variable_Symbol_Table(1:nodeCount,1);
[q_B,Price_B,BoosterLocationIndex,BoosterCount] = InitialBooster(nodeCount,Location_B,flowRate_B,NodeID,Price_B);

% Compute Quality without MSX
% (This function contains events) Don't uncomment this commands!!! Crash
% easily
qual_res = d.getComputedQualityTimeSeries; %Value x Node, Value x Link
LinkQuality = qual_res.LinkQuality;
NodeQuality = qual_res.NodeQuality;

% Initial Concentration
C0 = [NodeQuality(1,:) LinkQuality(1,:)];

%% Construct aux struct

aux = struct('NumberofSegment',NumberofSegment,...
    'LinkLengthPipe',LinkLengthPipe,...
    'LinkDiameterPipe',LinkDiameterPipe,...
    'TankBulkReactionCoeff',TankBulkReactionCoeff,...
    'TankMassMatrix',TankMassMatrix,...
    'JunctionMassMatrix',JunctionMassMatrix,...
    'MassEnergyMatrix',MassEnergyMatrix,...
    'flowRate_B',flowRate_B,...
    'q_B',q_B,...
    'Price_B',Price_B,...
    'NodeNameID',{NodeNameID},...
    'LinkNameID',{LinkNameID},...
    'NodesConnectingLinksID',{NodesConnectingLinksID},...
    'COMPARE',COMPARE);


%% Start MPC control

QsN_Control = [];
QsL_Control = [];
NodeSourceQuality = [];
T = [];
PreviousSystemDynamicMatrix = [];
UeachMin = [];
X_estimated = [];
PreviousDelta_t = [];
ControlActionU = [];
JunctionActualDemand = [];
Head = [];
Flow = [];
XX_estimated = [];
ControlActionU_LDE = [];
Magnitude = [];

Hq_min = Constants4Concentration.Hq_min;% I need that all concention 5 minutes later are  in 0.2 mg 4 mg
SimutionTimeInMinute = Constants4Concentration.SimutionTimeInMinute;

PreviousValue = struct('PreviousDelta_t',PreviousDelta_t,...
    'PreviousSystemDynamicMatrix',PreviousSystemDynamicMatrix,...
    'X_estimated',X_estimated,...
    'U_C_B_eachStep',0,...
    'UeachMinforEPANET',0);

d.openHydraulicAnalysis;
d.openQualityAnalysis;
d.initializeHydraulicAnalysis;
d.initializeQualityAnalysis;

tleft=1;
tInMin = 0;
delta_t = 0;


% profile on
tic
while (tleft>0 && tInMin < SimutionTimeInMinute && delta_t <= 60)
    t1 = d.runHydraulicAnalysis;
    t=d.runQualityAnalysis;
    
    % Obtain the actual Concentration
    QsN_Control=[QsN_Control; d.getNodeActualQuality];
    QsL_Control=[QsL_Control; d.getLinkQuality];
    Head=[Head; d.getNodeHydaulicHead];
    Flow=[Flow; d.getLinkFlows];
    TempDemand = d.getNodeActualDemand;
    JunctionActualDemand = [JunctionActualDemand; TempDemand(NodeJunctionIndex)];
    
    % Calculate Control Action
    tInMin = t/60;
    if(mod(tInMin,Hq_min)==0)
        % 5 miniute is up, Calculate the New Control Action
        disp('Current time')
        tInMin
        tInHour = tInMin/60
        CurrentVelocity = d.getLinkVelocity;
        
        CurrentVelocityPipe = CurrentVelocity(:,PipeIndex);
        PipeReactionCoeff = CalculatePipeReactionCoeff(CurrentVelocityPipe,LinkDiameterPipe,Kb_all,Kw_all,PipeIndex);
        % the minium step length for all pipes
        delta_t = LinkLengthPipe./NumberofSegment./CurrentVelocityPipe;
        delta_t = min(delta_t);
        delta_t = MakeDelta_tAsInteger(delta_t)
        
        CurrentFlow = d.getLinkFlows;
        CurrentHead = d.getNodeHydaulicHead;
        
        Volume = d.getNodeTankVolume;
        CurrentNodeTankVolume = Volume(NodeTankIndex);
        CurrentHead = d.getNodeHydaulicHead;
        
        % Estimate Hp of concentration; basciall 5 mins = how many steps
        SetTimeParameter = Hq_min*Constants4Concentration.MinInSecond/delta_t;
        Np = round(SetTimeParameter)
        
        % Collect the current value
        CurrentValue = struct('CurrentVelocityPipe',CurrentVelocityPipe,...
            'CurrentNodeTankVolume',CurrentNodeTankVolume,...
            'CurrentFlow',CurrentFlow,...
            'CurrentHead',CurrentHead,...
            'delta_t',delta_t,...
            'PipeReactionCoeff',PipeReactionCoeff,...
            'Np',Np);
        
        % Esitmate the concentration in all elements according to the
        % system dynamics each 5 mins
        %[x_estimated,xx_estimated] = EstimateState_XX(CurrentValue,IndexInVar,aux,ElementCount,q_B,tInMin,C0,PreviousValue);
        xx_estimated = EstimateState_XX_SaveMem(CurrentValue,IndexInVar,aux,ElementCount,q_B,tInMin,C0,PreviousValue);
        x_estimated = xx_estimated(:,end);
        
        XX_estimated = [XX_estimated xx_estimated];
        
        [A,B,C] = ObtainDynamicNew(CurrentValue,IndexInVar,aux,ElementCount,q_B);
        
        PreviousSystemDynamicMatrix = struct('A',A,...
            'B',B,...
            'C',C);
        
        % Calculate all of the control actions at each min
        
        PreviousDelta_t = [PreviousDelta_t delta_t];
        
        %PreviousDelta_t = [PreviousDelta_t delta_t];
        PreviousValue.PreviousDelta_t = delta_t;
        PreviousValue.PreviousSystemDynamicMatrix = PreviousSystemDynamicMatrix;
        PreviousValue.X_estimated = xx_estimated(:,end);
    end
    
    T=[T; t];
    tstep1 = d.nextHydraulicAnalysisStep;
    tstep = d.nextQualityAnalysisStep;
end
runningtime = toc
d.closeQualityAnalysis;
d.closeHydraulicAnalysis;

%% Start to plot
disp('Done!! Start to organize data')

% find average data;
X_Min = XX_estimated;
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
    
    % Epanet's way to define the pump and valve concnetration, wthich is
    % the average value of the upstream and downstream nodes
    
    NodeIndex4EachLink = findIndexofNode_Link(MassEnergyMatrix);
    NodeIndex4EachPump = NodeIndex4EachLink(PumpIndex,:);
    NodeIndex4EachValve = NodeIndex4EachLink(ValveIndex,:);
    for ithPump = 1:PumpCount
        X_Min_Average(PumpIndexInOrder,i) = (X_Min(NodeIndex4EachPump(ithPump,1),i) +  X_Min(NodeIndex4EachPump(ithPump,2),i))*0.5;
    end
    for ithValve = 1:ValveCount
        X_Min_Average(ValveIndexInOrder,i) = (X_Min(NodeIndex4EachValve(ithValve,1),i) +  X_Min(NodeIndex4EachValve(ithValve,2),i))*0.5;
    end
end
close all

NodeIndex = d.getNodeIndex;
LinkIndex = nodeCount+d.getLinkIndex;

X_Min_Average = X_Min_Average';
X_node_control_result =  X_Min_Average(:,NodeIndex);
X_link_control_result =  X_Min_Average(:,LinkIndex);
X_Junction_control_result =  X_Min_Average(:,NodeJunctionIndex);
NodeID4Legend = Variable_Symbol_Table2(NodeIndex,1);
LinkID4Legend = Variable_Symbol_Table2(LinkIndex,1);

figure
plot(QsN_Control);
legend(NodeID4Legend)
xlabel('Time (minute)')
ylabel('Concentrations at junctions (mg/L)')

figure
plot(QsL_Control);
legend(LinkID4Legend)
xlabel('Time (minute)')
ylabel('Concentrations in links (mg/L)')

figure
plot(X_node_control_result);
legend(NodeID4Legend)
xlabel('Time (minute)')
ylabel('Concentrations at junctions (mg/L)')

figure
plot(X_link_control_result);
legend(LinkID4Legend)
xlabel('Time (minute)')
ylabel('Concentrations in links (mg/L)')

figure
plot(JunctionActualDemand)
xlabel('Time (minute)')
ylabel('Demand at junctions (GPM)')
legend(NodeID4Legend)

figure
plot(Flow)
legend(LinkID4Legend)
xlabel('Time (minute)')
ylabel('Flow rates in links (GPM)')

disp('Summary:')
disp(['Compare is: ',num2str(COMPARE)]);

save(filename)

% remove all uncertainty and then resimulate to get the error between LDE
% and EPANET
epanetResult1 = [NodeQuality LinkQuality];
epanetResult = [QsN_Control QsL_Control];
% note that the above are both from EPAENT, their result are a little bit
% different. I believe this because they use different method to implement
% this.

LDEResult = [X_node_control_result X_link_control_result];
LDEResult1 = updateLDEResult(LDEResult,IndexInVar,MassEnergyMatrix);

Calculate_Error_EPANET_LDE(epanetResult1',LDEResult1');

%plotSensorSelectionResult(sensorSelectionResult,NodeID4Legend);
[m,n] = size(X_Min);
basePipeCIndex = min(Pipe_CIndex);
First = basePipeCIndex:basePipeCIndex+NumberofSegment-1;
for j = 1:PipeCount
    Indexrange = (j-1)*NumberofSegment + First;
    figure
    imagesc(X_Min(Indexrange,:));
    pipeIDTtile = LinkID4Legend{j}
    title(pipeIDTtile)
    colorbar;
end



