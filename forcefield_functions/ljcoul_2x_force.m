function [r,lj,coul,Ftot] = ljcoul_2x_force(param,varargin)

if nargin > 1
    r=varargin{1};
else
    %     r=1:length(data);
    r=.01:.0005:1.2; % nm
end

q11=param(1);
q12=param(2);
q2=param(3);

sig11=param(4);
sig12=param(5);
sig2=param(6);

eps11=param(7);
eps12=param(8);
eps2=param(9);

e_mix1=(eps11*eps2)^.5;
e_mix2=(eps12*eps2)^.5;
sig_mix1=(sig11+sig2)/2;
sig_mix2=(sig12+sig2)/2;

lj1=4*e_mix1.*((sig_mix1./r).^12-(sig_mix1./r).^6);
lj2=4*e_mix2.*((sig_mix2./r).^12-(sig_mix2./r).^6);
coul1=(1.60217646E-19)^2*6.022E+23*q11*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;
coul2=(1.60217646E-19)^2*6.022E+23*q12*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;

lj1=-diff(lj1)/(r(2)-r(1));
lj2=-diff(lj2)/(r(2)-r(1));
coul1=-diff(coul1)/(r(2)-r(1));
coul2=-diff(coul2)/(r(2)-r(1));

Ftot=lj1+lj2+coul1+coul2;

lj=lj1+lj2;
coul=coul1+coul2;


hold on

plotmin=1000*(round2dec(min(Ftot*1.5)/1000));
if plotmin>=0
    plotmin=10000;
end

plot(r(2:end),lj,'b','LineWidth',1);
plot(r(2:end),coul,'r','LineWidth',1);
plot(r(2:end),Ftot,'k','LineWidth',1);
xlabel('r [nm]');
ylabel('F [kJ mol^-^1 nm^-^1]');
xlim([0,1.2]);
ylim(sort([-plotmin plotmin]))
