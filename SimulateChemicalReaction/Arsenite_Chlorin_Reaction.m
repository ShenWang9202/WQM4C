function dmf=Arsenite_Chlorin_Reaction(t,X)
% Constant
K2 = 69.4;        %100000/24/60= Arsenite oxidation velocity;
MWChlorine = 70.9;        %mg/mmol free chlorine as Cl2
MWAsIII = 74.92;       %mg/mmol arsenite as As
MWAsV = 74.92 ;      %mg/mmol arsenate as As
% Reaction Rate
K_Cl_decay = 6.94444E-4;  %1/24/60=Free chlorine bulk decay rate == 1/day
K_AsIII_Cl = K2/MWAsIII;
K_Cl_AsIII = K2/MWChlorine;
K_AsV = K2*MWAsV/MWChlorine/MWAsIII;



% Initialization
Chlorine = X(1);
AsIII = X(2);
AsV = X(3);
% Reaction
dChlorine = -K_Cl_decay*Chlorine - K_AsIII_Cl*Chlorine*AsIII;
dAsIII = -K_Cl_AsIII*Chlorine*AsIII;
dAsV = K_AsV*Chlorine*AsIII;
% Results
dmf=[dChlorine;dAsIII;dAsV];