function [r,lj,coul,Ftot] = ljcoul_force(param,varargin)

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


if nargin > 7
    LineWidth=varargin{2};
else
    LineWidth=0.5;
end

if numel(param)<7
    
    e_mix=(eps1*eps2)^.5;
    sig_mix=(sig1+sig2)/2;
    
    lj=4*e_mix.*((sig_mix./r).^12-(sig_mix./r).^6);
    coul=(1.60217646E-19)^2*6.022E+23*q1*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;
    
else
    
    C61=4*eps1*sig1^6;
    C121=4*eps1*sig1^12;
    C62=4*eps2*sig2^6;
    C122=4*eps2*sig2^12;
    
    C12_mix=(C121*C122)^.5;
    C6_mix=(C61*C62)^.5;
    C4_mix=(C41*C42)^.5;
    
    lj=C12_mix./r.^12-C6_mix./r.^6-C4_mix./r.^4;
    coul=(1.60217646E-19)^2*6.022E+23*q1*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;
    
end

lj=-diff(lj)/(r(2)-r(1));
coul=-diff(coul)/(r(2)-r(1));
Ftot=lj+coul;

hold on

plotmin=1000*(round2dec(min(Ftot*1.5)/1000));
if plotmin>=0
    plotmin=10000;
end

if nargin>2
    %% Plot the Total energy, electrostatic contribution, and the LJ
    
    if nargin>3
        color1=varargin{3};
        color2=varargin{3};
        color3=varargin{3};
    else
        color1='b--';
        color2='r--';
        color3='k--';
    end
    
    hold on
    plotmin=1000*(round2dec(min(Ftot*1.5)/1000));
    if plotmin>=0
        plotmin=10000;
    end
    
    
    plot(r(2:end),lj,color1,'LineWidth',1);
    plot(r(2:end),coul,color2,'LineWidth',1);
    plot(r(2:end),Ftot,color3,'LineWidth',1);
    
    xlabel('r [nm]');
    ylabel('F [kJ mol^-^1 nm^-^1]');
    xlim([0,1.2]);
    ylim(sort([-plotmin plotmin]))
    
end
