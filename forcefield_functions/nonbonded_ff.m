function [r,lj,coul,Utot,q1,q2,sig1,sig2,eps1,eps2] = nonbonded_ff(ff,atomtype,varargin) % ff and atomtype can be single variables or cell 'tuplets'

if ischar(atomtype)
    disp('atomtype must be a {1x1} or {1xN} cell')
end

if size(atomtype,2)==1
    atomtype1=atomtype;
    atomtype2=atomtype;
elseif size(atomtype,2)==2
    atomtype1=atomtype(1);
    atomtype2=atomtype(2);
end

if size(ff,2)>2
    ff1=ff;
    ff2=ff;
    atomtype1=ff(strcmpi([ff1.type],atomtype1));
    atomtype2=ff(strcmpi([ff2.type],atomtype2));
elseif size(ff,2)==2
    ff1=ff{1};
    ff2=ff{2};
    atomtype1=ff1(strcmpi([ff1.type],atomtype1));
    atomtype2=ff2(strcmpi([ff2.type],atomtype2));
end

q1=atomtype1.charge;
q2=atomtype2.charge;
sig1=atomtype1.sigma_nm;
sig2=atomtype2.sigma_nm;
eps1=atomtype1.e_kJmol;
eps2=atomtype2.e_kJmol;

[q1 q2 sig1 sig2 eps1 eps2]

r=.01:.0005:1.2; % nm

e_mix=(eps1*eps2)^.5
sig_mix=(sig1+sig2)/2

lj=4*e_mix.*((sig_mix./r).^12-(sig_mix./r).^6);
coul=(1.60217646E-19)^2*6.022E+23*q1*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;
Utot=lj+coul;

if nargin>2
    %% Plot the Total energy, electrostatic contribution, and the LJ
    
    
    hold on
    plotmin=1000*(round(min(Utot*1.5)/1000));
    if plotmin>=0
        plotmin=10000;
    end
    
    if nargin>2
        plot(r,lj,'b','LineWidth',1);
        plot(r,coul,'r','LineWidth',1);
        plot(r,Utot,'k','LineWidth',1);
    else
        plot(r,lj,'LineWidth',1);
        plot(r,coul,'LineWidth',1);
        plot(r,Utot,'LineWidth',1);
    end
    xlabel('r [nm]');
    ylabel('U [kJ/mol]');
    xlim([0,1.2]);
    ylim(sort([-plotmin plotmin]))
    
    [atomtype1.type]
    
    titlestring=strcat([atomtype1.type],'-',[atomtype2.type]);
    titlestring=strrep(titlestring,'_','-');
    title(titlestring)
    
    % figure
    % hold on
    % P=exp(-Utot./2.4789);
    % P=P./sum(P);
    % plot(r,P)
    % xlabel('r [nm]');
    % ylabel('P');
    %
    % assignin('caller','P',P);
    %
    % gaussEqn='a*exp(-((x-b)/(2^.5*c))^2)';%+d';%
    % startPoints=[0.01 0.15 .25];
    % gaussEqn='a*exp(-((x-b)/c)^2)+d'; % Orig from Matlab
    % startPoints=[0.01 0.15 .25 0];
    % fgauss=fit(r',P',gaussEqn,'Start',startPoints);
    % coeff=coeffvalues(fgauss);
    % scale=coeff(1);
    % ave_r=coeff(2);
    % fwhm=2*(2*log(2))^.5*abs(coeff(3));
    % sumP=sum(coeff(1)*exp(-((r-coeff(2))./(2^.5*coeff(3))).^2));%+coeff(4));
    % plot(fgauss,r,P)
    % legend('Location','northeast')
    % xlim([0,0.35]);
    % hold off
    
end


