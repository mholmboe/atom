function [r,buck,coul,Utot] = buckinghamcoul(param,varargin)

q1=param(1);
q2=param(2);

A1=param(3);
A2=param(4);

B1=param(5);
B2=param(6);

C1=param(7);
C2=param(8);

if nargin > 1
    r=varargin{1};
else
    %     r=1:length(data);
    r=.01:.0005:1.2; % nm
end


%% Combination rules from Gromacs documentation
A_mix=(A1*A2)^.5;
B_mix=2/(1/B1+1/B2);
C_mix=(C1*C2)^.5;

buck=A_mix * exp(-B_mix * r) - C_mix./r.^6;
coul=(1.60217646E-19)^2*6.022E+23*q1*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;

Utot=buck+coul;

hold on

plotmin=1000*(round(min(Utot*1.5)/1000));
if plotmin>=0
    plotmin=10000;
end
plot(r,buck,'b.');
plot(r,coul,'r.');
plot(r,Utot,'k.');
% xlabel('r [nm]');
% ylabel('U [kJ/mol]');
% xlim([0,1.2]);
% ylim(sort([-plotmin plotmin]))
 ylim(sort([-5000 5000]))