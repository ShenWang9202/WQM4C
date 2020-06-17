%% Hydraulic and Quality analysis
% This function contains:
% Load a network
% Hydraulic and Quality analysis STEP-BY-STEP


% Pattern Uncertainty test

% System: only in Windows
% Author: Shen Wang
% Date: 3/4/2020
%%
clear; close('all'); clc;
start_toolkit;

% Load a network
d = epanet('Threenode-cl-2-paper.inp');
    %Uncertainty
    qunc=0.1;
% Set time hydraulic and quality steps
% etstep = 300;
% d.setTimeReportingStep(etstep);
% d.setTimeHydraulicStep(etstep);
% d.setTimeQualityStep(etstep);
% Hstep = min(Pstep,Hstep)
% Hstep = min(Rstep,Hstep)
% Hstep = min(Qstep,Hstep)
%% Prepare for the demand uncertainty


% Hydraulic analysis using the functions ENopenH, ENinit, ENrunH, ENgetnodevalue/&ENgetlinkvalue, ENnextH, ENcloseH
% (This function contains events)
%hyd_res_3 = d.getComputedHydraulicTimeSeries;

JunctionCount = d.getNodeJunctionCount;

NodeJunctionIndex = d.getNodeJunctionIndex;
NodePatternIndex = d.getNodePatternIndex;
JunctionPatternIndex = NodePatternIndex(NodeJunctionIndex);

Patterns = d.getPattern;
JunctionPattern = Patterns(JunctionPatternIndex,:);


JunctionPattern_Uncertainty = add_uncertainty(JunctionPattern,qunc);

Demand_known = d.getNodeActualDemand;
Demand_known = (Demand_known(NodeJunctionIndex));
d.setPattern(JunctionPatternIndex,JunctionPattern_Uncertainty);

% Hydraulic and Quality analysis STEP-BY-STEP
d.openHydraulicAnalysis;
d.openQualityAnalysis;
d.initializeHydraulicAnalysis(0);
d.initializeQualityAnalysis(d.ToolkitConstants.EN_NOSAVE);





tstep = 1;
T = []; P = []; F = []; QN = []; multiplier = []; Demand_all =[];
while (tstep>0)
    t = d.runHydraulicAnalysis
    d.getNodeActualDemand
    d.getNodeBaseDemands
    
    qt = d.runQualityAnalysis
    
    P = [P; d.getNodePressure];
    F = [F; d.getLinkFlows];
    
    QN = [QN; d.getNodeActualQuality];
    T = [T; t];
    
    tstep = d.nextHydraulicAnalysisStep;
    qtstep = d.nextQualityAnalysisStep;
    
    Demand_known = d.getNodeActualDemand
    Demand_all = [Demand_all; Demand_known];
%     baseDemand_cell = d.getNodeBaseDemands;                                                                                                                
%     baseDemand = baseDemand_cell{1};
%     
%     rand_a = -5;
%     rand_b = 5;
%     rand_noise_percent = (rand_b-rand_a).*rand(1,3) + rand_a;
%     rand_noise_percent = rand_noise_percent/100
%     baseDemand_noise = baseDemand.*(1+rand_noise_percent);
%     baseDemand_cell{1} = baseDemand_noise;
    %d.setNodeBaseDemands(baseDemand_cell);
    Demand_noise = d.getNodeActualDemand
    
    
%     
%     rand_a = -20;
%     rand_b = 20;
%     for i = 1:JunctionCount
%         ind = NodeJunctionIndex(i);
%         if (abs(junctionDemand(ind)) >= 1e-4)
%             mult = Demand_known(ind) / junctionDemand(ind);
%             
%             multiplier = [multiplier mult];
%         else
%             multiplier = [multiplier 0];
%         end
%     end
    
end
d.closeQualityAnalysis;
d.closeHydraulicAnalysis;
%
% P
% F
% QN

% % Unload library
% d.unload;
function value_unc = add_uncertainty(qext, qunc)
    ql=qext-qunc*qext;
    qu=qext+qunc*qext;
    value_unc=ql+rand(1,length(qext)).*(qu-ql);
end

