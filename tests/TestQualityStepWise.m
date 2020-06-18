clear all
clear; clc;  close all
%% Load EPANET MATLAB TOOLKIT
start_toolkit;
%% run EPANET MATLAB TOOLKIT to obtain data

NetworkName = 'Fournode-Cl-As-2.inp';
%MSXName = 'Arsenite.msx';
MSXName = 'Threenode-cl-3.msx';% Load MSX file, only Chlorine

%% Prepare constants data for MPC
PrepareData4Control
d.loadMSXFile(MSXName); % Load MSX file,mulit-species
%% initialize concentration at nodes

nx = NumberofX; % Number of states

% initialize BOOSTER
% flow of Booster, assume we put booster at each nodes, so the size of it
% should be the number of nodes.
JunctionCount = double(JunctionCount);
ReservoirCount = double(ReservoirCount);
TankCount = double(TankCount);
nodeCount = JunctionCount + ReservoirCount + TankCount;

NodeID = Variable_Symbol_Table(1:nodeCount,1);

% Compute Quality without MSX
% (This function contains events) Don't uncomment this commands!!! Crash
% easily
qual_res = d.getComputedQualityTimeSeries; %Value x Node, Value x Link
LinkQuality = qual_res.LinkQuality;
NodeQuality = qual_res.NodeQuality;

C0 = [NodeQuality(1,:) LinkQuality(1,:)];



%% Start MPC control
T = [];  

JunctionActualDemand = []; Head = []; Flow = []; 

Hq_min = Constants4Concentration.Hq_min;% I need that all concention 5 minutes later are  in 0.2 mg 4 mg
SimutionTimeInMinute = Constants4Concentration.SimutionTimeInMinute;

NodeCount=1:d.getNodeCount;%index node
LinkCount=1:d.getLinkCount;%index link
SpeciesCount=1:d.getMSXSpeciesCount;

uu = SpeciesCount
ss = NodeCount
value.NodeQuality = cell(1, length(NodeCount));
value.LinkQuality = cell(1, length(LinkCount));
% Obtain a hydraulic solution   Don't uncomment this commands!!! Crash
d.solveMSXCompleteHydraulics

k=1; tleft=1;t=0;
value.Time(k, :)=0;
time_step = d.getMSXTimeStep;
timeSmle=d.getTimeSimulationDuration;%bug at time

tleft=1;
tInMin = 0;
delta_t = 0;

d.openHydraulicAnalysis;
d.openQualityAnalysis;
d.initializeHydraulicAnalysis(0);
d.initializeQualityAnalysis(0);
tic
while (tleft>0 && d.Errcode==0 && timeSmle~=t && delta_t <= 60)
    d.runHydraulicAnalysis
    [t, tleft]=d.stepMSXQualityAnalysisTimeLeft;
    
    Head=[Head; d.getNodeHydaulicHead];
    Flow=[Flow; d.getLinkFlows];
    TempDemand = d.getNodeActualDemand;
    JunctionActualDemand = [JunctionActualDemand; TempDemand(NodeJunctionIndex)];
    
    [NodeQuality,LinkQuality] = ObtainSpeciesConcentration(d,t,time_step,NodeCount,LinkCount,SpeciesCount);
    
    % Obtain the actual Concentration
    for i = NodeCount
        value.NodeQuality{i}(k,:) = NodeQuality{i};
    end
    for i = LinkCount
        value.LinkQuality{i}(k,:) = LinkQuality{i};
    end
    
    if k>1
        value.Time(k, :)=t;
    end
    k = k + 1;
    T=[T; t];
    
    
    tsstep = d.nextHydraulicAnalysisStep
    %    qststep = d.nextQualityAnalysisStep
end
d.closeQualityAnalysis;
d.closeHydraulicAnalysis;
toc

QNode = d.getMSXComputedQualityNode
QLink = d.getMSXComputedQualityLink

close all
NodeID4Legend = Variable_Symbol_Table2(d.getNodeIndex,1);
LinkID4Legend = Variable_Symbol_Table2(nodeCount+d.getLinkIndex,1);


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