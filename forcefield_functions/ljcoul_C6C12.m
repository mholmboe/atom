function [r,lj,coul,Utot] = ljcoul_C6C12(param,varargin)

q1=param(1);
q2=param(2);

c61=param(3);
c62=param(4);

c121=param(5);
c122=param(6);

if numel(param)>6
    C41=param(7);
    C42=param(7);
end

if numel(param)>7
    C41=param(7);
    C42=param(8);
end

if nargin > 1
    r=varargin{1};
else
    %     r=1:length(data);
    r=.01:.0005:1.2; % nm
end

C6_mix=(c61*c62)^.5;
C12_mix=(c121*c122)^.5;

lj=C12_mix./r.^12-C6_mix./r.^6;
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