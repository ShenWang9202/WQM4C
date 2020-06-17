%%  Case Study for CCWI2016

close all; clear; 

%% Load EPANET Network and MSX
G = epanet('Network_1.inp'); % Load EPANET Input file
G.loadMSXFile('Arsenite.msx'); % Load MSX file


close all; clear; 
G = epanet('Network_2.inp'); % Load EPANET Input file
G.loadMSXFile('Arsenite.msx'); % Load MSX file

