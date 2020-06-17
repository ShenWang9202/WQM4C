% Create scenarios based on node count of inp file EPANET
% For more check out the S-PLACE Toolkit (Sensor placement Toolkit)
% https://github.com/KIOS-Research/splace-toolkit


% clear; close('all'); clc;
% start_toolkit;
% 
% % Load a network
% d = epanet('Threenode-cl-booster-mass.inp');
% % Compute Quality without MSX
% % (This function contains events) Don't uncomment this commands!!! Crash
% % easily
% qual_res = d.getComputedQualityTimeSeries; %Value x Node, Value x Link
% LinkQuality = qual_res.LinkQuality;
% TimeQualityStep = d.getTimeQualityStep;
% %% get info from EPANET
% % get index info from EPANET
% %PipeIndex = 1:d.getLinkPipeCount;
% PipeIndex = d.getLinkPipeIndex;
% PumpIndex = d.getLinkPumpIndex;
% ValveIndex = d.getLinkValveIndex;
% NodeJunctionIndex = d.getNodeJunctionIndex;
% ReservoirIndex = d.getNodeReservoirIndex;
% NodeTankIndex = d.getNodeTankIndex;
% 
% % get LinkDiameter from EPANET
% LinkDiameter = d.getLinkDiameter;
% LinkDiameterPipe = LinkDiameter(PipeIndex);
% LinkLength = d.getLinkLength;
% LinkLengthPipe = LinkLength(PipeIndex);
% 
% LinkQualityPipe = LinkQuality(:,PipeIndex);
% %%
% 
% d.openQualityAnalysis
% d.initializeQualityAnalysis
% tleft=1; P=[];T=[];QsN=[]; QsL=[]; Velocity=[]; Head=[];
% Flow=[]; JunctionActualDemand=[]; NodeTankVolume = [];NodeNetFlowTank = [];
% while (tleft>0)
%     %Add code which changes something related to quality
%     t=d.runQualityAnalysis;
%     P=[P; d.getNodePressure];
%     Head=[Head; d.getNodeHydaulicHead];
%     Flow=[Flow; d.getLinkFlows];
%     Velocity = [Velocity; d.getLinkVelocity];
%     TempDemand = d.getNodeActualDemand;
%     JunctionActualDemand = [JunctionActualDemand; TempDemand(NodeJunctionIndex)];
%     NodeNetFlowTank = [NodeNetFlowTank; TempDemand(NodeTankIndex)];
%     Volume = d.getNodeTankVolume;
%     NodeTankVolume = [NodeTankVolume; Volume(NodeTankIndex)];
%     QsN=[QsN; d.getNodeActualQuality];
%     QsL=[QsL; d.getLinkQuality];
%     T=[T; t];
%     tleft = d.stepQualityAnalysisTimeLeft;
% end
% 
% % verify C_junction = [sum(m_pipes_in) + m_booster]/[sum(flow_in)]
% % unit C_junction: mg/L 
% %      m_booster: mg/min
% %      flow_in: LperMin  (1GPM = Constants4Concentration.Gallon2Liter LPM)
% 
% 
% PatternIndex = d.getNodeSourcePatternIndex;
% PatternIndexofJunction2 = PatternIndex(1);
% Patterns = d.getPattern;
% PatternofJunction2 = Patterns(PatternIndexofJunction2,:);
% NodeSourceQuality = d.getNodeSourceQuality;
% SourceQualityofJunction2 = NodeSourceQuality(1);
% 
% MassAtJunction2 = SourceQualityofJunction2*PatternofJunction2;
% ConcentrationAtJunction2 = [];
% for i = 0:1440-1
%     concentrationAtJunction2 = MassAtJunction2(floor(i/60)+1)/(Flow(i+1,2)* Constants4Concentration.Gallon2Liter);
%     ConcentrationAtJunction2 = [ConcentrationAtJunction2 concentrationAtJunction2];
% end
% 


%% 

clear; close('all'); clc;

% Load a network
d = epanet('Threenode-cl-booster-mass0.inp');
% Compute Quality without MSX
% (This function contains events) Don't uncomment this commands!!! Crash
% easily
qual_res = d.getComputedQualityTimeSeries; %Value x Node, Value x Link
LinkQuality = qual_res.LinkQuality;
TimeQualityStep = d.getTimeQualityStep;
%% get info from EPANET
% get index info from EPANET
%PipeIndex = 1:d.getLinkPipeCount;
PipeIndex = d.getLinkPipeIndex;
PumpIndex = d.getLinkPumpIndex;
ValveIndex = d.getLinkValveIndex;
NodeJunctionIndex = d.getNodeJunctionIndex;
ReservoirIndex = d.getNodeReservoirIndex;
NodeTankIndex = d.getNodeTankIndex;

% get LinkDiameter from EPANET
LinkDiameter = d.getLinkDiameter;
LinkDiameterPipe = LinkDiameter(PipeIndex);
LinkLength = d.getLinkLength;
LinkLengthPipe = LinkLength(PipeIndex);

LinkQualityPipe = LinkQuality(:,PipeIndex);
%%


d.setNodeSourceType(1, 'MASS'); %Junction2's index is 1; we set it as mass booster
SourcePattern = d.getNodeSourcePatternIndex;
SourcePattern = [3 0 0]; % set the third pattern
d.setNodeSourcePatternIndex(SourcePattern); 
d.getNodeSourcePatternIndex
TmpNodeSourceQuality = d.getNodeSourceQuality;
TmpNodeSourceQuality(1) = 5000;
d.setNodeSourceQuality(TmpNodeSourceQuality)


d.openQualityAnalysis
d.initializeQualityAnalysis
tleft=1; P=[];T=[];QsN=[]; QsL=[]; Velocity=[]; Head=[];
Flow=[]; JunctionActualDemand=[]; NodeTankVolume = [];NodeNetFlowTank = [];
while (tleft>0)
    %Add code which changes something related to quality
    t=d.runQualityAnalysis;
    P=[P; d.getNodePressure];
    Head=[Head; d.getNodeHydaulicHead];
    Flow=[Flow; d.getLinkFlows];
    Velocity = [Velocity; d.getLinkVelocity];
    TempDemand = d.getNodeActualDemand;
    JunctionActualDemand = [JunctionActualDemand; TempDemand(NodeJunctionIndex)];
    NodeNetFlowTank = [NodeNetFlowTank; TempDemand(NodeTankIndex)];
    Volume = d.getNodeTankVolume;
    NodeTankVolume = [NodeTankVolume; Volume(NodeTankIndex)];
    QsN=[QsN; d.getNodeActualQuality];
    QsL=[QsL; d.getLinkQuality];
    T=[T; t];
    tleft = d.stepQualityAnalysisTimeLeft;
end

% verify C_junction = [sum(m_pipes_in) + m_booster]/[sum(flow_in)]
% unit C_junction: mg/L 
%      m_booster: mg/min
%      flow_in: LperMin  (1GPM = Constants4Concentration.Gallon2Liter LPM)









