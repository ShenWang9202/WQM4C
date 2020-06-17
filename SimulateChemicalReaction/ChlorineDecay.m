function dmf=ChlorineDecay(t,X)
Chlorine=X(1);
Kb = 1/24/60;
dChlorine=-Kb*Chlorine;
dmf=[dChlorine;];