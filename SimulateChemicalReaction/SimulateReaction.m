[ts,data] = ode45(@enzyme,0:0.01:10,[30 5 0 0]);
plot(ts,data(:,1),'r',ts,data(:,2),'b',ts,data(:,3),'g',ts,data(:,4),'m')
legend('[S]','[E]','[ES]','[P]')
f=[ts,data];


[ts,data] = ode45(@ChlorineDecay,0:1:60*24,[1]); % 10 minus initial value 0.8
plot(ts,data)
legend('[Cl]')
f=[ts,data];



[ts,data] = ode45(@Arsenite_Chlorin_Reaction,0:1:60,[0.8,10,0]); % 10 minus initial value 0.8
plot(ts,data(:,1),'r',ts,data(:,2),'b',ts,data(:,3),'g')
legend('[Cl]','[AsIII]','[AsV]')
f=[ts,data];