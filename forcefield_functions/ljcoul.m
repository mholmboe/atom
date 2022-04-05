function [r,lj,coul,Utot] = ljcoul(param,varargin)

% clear all;
% format compact;

% % 	Sigma nm	e kJ/mol
% % O	0.31656	    0.65017
% % M	0.154749107	0.3071056

% 0.407222	3.77671E-05

% q1=-1.05;
% q2=2.1;
% % q2=2.0;
% %
% sig1=0.31656; %nm
% %sig2=0.407222; %nm
% sig2=0.3302; %nm
% %
% eps1=0.65017;% kJ/mol
% %eps2=3.77671E-05;% kJ/mol
% eps2=7.70065E-06;% kJ/mol

q1=param(1);
q2=param(2);

sig1=param(3);
sig2=param(4);

eps1=param(5);
eps2=param(6);

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

if numel(param)<7
    
    e_mix=(eps1*eps2)^.5;
    sig_mix=(sig1+sig2)/2;
    
    lj=4*e_mix.*((sig_mix./r).^12-(sig_mix./r).^6);
    
else
    
    C61=4*eps1*sig1^6;
    C121=4*eps1*sig1^12;
    C62=4*eps2*sig2^6;
    C122=4*eps2*sig2^12;
    
    C12_mix=(C121*C122)^.5;
    C6_mix=(C61*C62)^.5;
    C4_mix=(C41*C42)^.5;
    
    lj=C12_mix./r.^12-C6_mix./r.^6-C4_mix./r.^4;
end

coul=(1.60217646E-19)^2*6.022E+23*q1*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;    

Utot=lj+coul;

hold on

plotmin=1000*(round(min(Utot*1.5)/1000));
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