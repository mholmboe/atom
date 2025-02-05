function [r,lj,coul,Utot] = ljcoul_2x(param,varargin)

q11=param(1);
q12=param(2);
q2=param(3);

sig11=param(4);
sig12=param(5);
sig2=param(6);

eps11=param(7);
eps12=param(8);
eps2=param(9);

if nargin > 1
    r=varargin{1};
else
    %     r=1:length(data);
    r=.01:.0005:1.2; % nm
end

e_mix1=(eps11*eps2)^.5;
e_mix2=(eps12*eps2)^.5;
sig_mix1=(sig11+sig2)/2;
sig_mix2=(sig12+sig2)/2;

lj1=4*e_mix1.*((sig_mix1./r).^12-(sig_mix1./r).^6);
lj2=4*e_mix2.*((sig_mix2./r).^12-(sig_mix2./r).^6);
coul1=(1.60217646E-19)^2*6.022E+23*q11*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;
coul2=(1.60217646E-19)^2*6.022E+23*q12*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;

lj=lj1+lj2;
coul=coul1+coul2;

Utot1=lj1+coul1;
Utot2=lj2+coul2;
Utot=Utot1+Utot2;

hold on

plotmin=1000*(round2dec(min(Utot*1.5)/1000));
if plotmin>=0
    plotmin=10000;
end
plot(r,lj1,'b--');
plot(r,lj2,'k--');
plot(r,lj,'b');
plot(r,coul,'r');
% plot(r,Utot1,'k');
% plot(r,Utot2,'k');
plot(r,Utot,'k');
xlabel('r [nm]');
ylabel('U [kJ/mol]');
xlim([0,1.2]);
ylim(sort([-plotmin plotmin]))
