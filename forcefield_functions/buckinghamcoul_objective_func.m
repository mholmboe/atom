function res = buckinghamcoul_objective_func(param,data,r,scalefactors,varargin)

param=param./scalefactors;

q1=param(1);
q2=param(2);

A1=param(3);
A2=param(4);

B1=param(5);
B2=param(6);

C1=param(7);
C2=param(8);

%% Combination rules from Gromacs documentation
A_mix=(A1*A2)^.5;
B_mix=2/(1/B1+1/B2);
C_mix=(C1*C2)^.5;

buck=A_mix * exp(-B_mix * r) - C_mix./r.^6;
coul=(1.60217646E-19)^2*6.022E+23*q1*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;

Utot=buck+coul;

res=data-Utot;
