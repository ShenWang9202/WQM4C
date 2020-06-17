function dmf=enzyme(t,X)
S=X(1);
E=X(2);
ES=X(3);
P=X(4);
dS=-10*S+2*ES*S+ES;
dE=2.5*ES-10*S+2*ES*S;
dES=10*S-2*ES*S-2.5*ES;
dP=1.5*ES;
dmf=[dS;dE;dES;dP];