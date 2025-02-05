function [r,lj,coul,Utot] = ljcoul_12_6(param,varargin)

q1=param(1);
q2=param(2);

sig1=param(3);
sig2=param(4);

eps1=param(5);
eps2=param(6);

if numel(param)>6
    disp('Did you mean the 12-6-4 potential??')
    pause
end

if nargin > 1
    r=varargin{1};
else
    %     r=1:length(data);
    r=.01:.0005:1.2; % nm
end


e_mix=(eps1*eps2)^.5;
sig_mix=(sig1+sig2)/2;

lj=4*e_mix.*((sig_mix./r).^12-(sig_mix./r).^6);
coul=(1.60217646E-19)^2*6.022E+23*q1*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;


Utot=lj+coul;

hold on

plotmin=1000*(round2dec(min(Utot*1.5)/1000));
if plotmin>=0
    plotmin=10000;
end
plot(r,lj,'b');
plot(r,coul,'r');
plot(r,Utot,'k');
xlabel('r [nm]');
ylabel('U [kJ/mol]');
xlim([0,1.2]);
ylim(sort([-plotmin plotmin]))