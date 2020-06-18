clear; clc;  close all
start_toolkit;

% Load EPANET Network and MSX
d = epanet('Threenode-cl-2.inp'); % Load EPANET Input file

sensor_id = {'2'};
sensor_index = d.getNodeIndex(sensor_id);
node_id = d.getNodeNameID;

%% Run
H = d.getComputedHydraulicTimeSeries;
Flow = H.Flow;
Head = H.Head;

Q = d.getComputedQualityTimeSeries; % Solve hydraulics and MSX quality dynamics
Q.NodeQuality
Q.LinkQuality