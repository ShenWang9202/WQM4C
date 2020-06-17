function dmf=Arsenite_Chlorin_Reaction(t,X)
KB = 6.94444E-4;  %1/24/60=Free chlorine bulk decay rate == 1/day
K2 = 69.4;        %100000/24/60= Arsenite oxidation velocity;
MWChlorine = 70.9;        %mg/mmol free chlorine as Cl2
MWAsIII = 74.92;       %mg/mmol arsenite as As
MWAsV = 74.92 ;      %mg/mmol arsenate as As
Chlorine = X(1);
AsIII = X(2);
AsV = X(3);
%dChlorine=-6.9444444E-4*Chlorine;
dChlorine = -KB*Chlorine - K2*Chlorine*AsIII/MWAsIII;
dAsIII = -K2*Chlorine*AsIII/MWChlorine;
dAsV = K2*MWAsV*Chlorine*AsIII/MWChlorine/MWAsIII;
dmf=[dChlorine;dAsIII;dAsV];