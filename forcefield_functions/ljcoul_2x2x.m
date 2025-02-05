function [r,lj,coul,Utot] = ljcoul_2x2x(param,varargin)

if nargin > 1
    r=varargin{1};
else
    %     r=1:length(data);
    r=.01:.0005:1.2; % nm
end

q11=param(1);
q12=param(2);
q21=param(3);
q22=param(4);

sig11=param(5);
sig12=param(6);
sig21=param(7);
sig22=param(8);

eps11=param(9);
eps12=param(10);
eps21=param(11);
eps22=param(12);

e_mix11=(eps11*eps21)^.5;
% e_mix21=(eps12*eps21)^.5;
% e_mix12=(eps11*eps22)^.5;
e_mix22=(eps12*eps22)^.5;
sig_mix11=(sig11+sig21)/2;
% sig_mix21=(sig12+sig21)/2;
% sig_mix12=(sig11+sig22)/2;
sig_mix22=(sig12+sig22)/2;

lj11=4*e_mix11.*((sig_mix11./r).^12-(sig_mix11./r).^6);
% lj21=4*e_mix21.*((sig_mix21./r).^12-(sig_mix21./r).^6);
% lj12=4*e_mix12.*((sig_mix12./r).^12-(sig_mix12./r).^6);
lj22=4*e_mix22.*((sig_mix22./r).^12-(sig_mix22./r).^6);
coul11=(1.60217646E-19)^2*6.022E+23*q11*q21./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;
% coul21=(1.60217646E-19)^2*6.022E+23*q12*q21./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;
% coul12=(1.60217646E-19)^2*6.022E+23*q11*q22./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;
coul22=(1.60217646E-19)^2*6.022E+23*q12*q22./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;

lj=lj11+lj22;
coul=coul11+coul22;
% Utot=lj11+lj21+lj12+lj22+coul11+coul21+coul12+coul22;
Utot=lj+coul;

hold on

plotmin=1000*(round2dec(min(Utot*1.5)/1000));
if plotmin>=0
    plotmin=1000;
end
plot(r,lj,'b');
plot(r,coul,'r');
% plot(r,Utot1,'k');
% plot(r,Utot2,'k');
plot(r,Utot,'k');
xlabel('r [nm]');
ylabel('U [kJ/mol]');
xlim([0,1.2]);
ylim(sort([-plotmin plotmin]))
