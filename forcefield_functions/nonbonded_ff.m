function [r,lj,coul,Utot,q1,q2,sig1,sig2,eps1,eps2,C4] = nonbonded_ff(ff,atomtype,varargin) % ff and atomtype can be single variables or cell 'tuplets'

if ischar(atomtype)
    disp('atomtype must be a {1x1} or {1xN} cell')
end

if size(atomtype,2)==1
    atomtype1=atomtype;
    atomtype2=atomtype;
elseif size(atomtype,2)==2
    atomtype1=atomtype(1);
    atomtype2=atomtype(2);
elseif size(atomtype,2)>2
    disp('Too many atomtypes')
end

if size(ff,2)==1
    ff=ff.ff;
end

if size(ff,2)>2
    ff1=ff;
    ff2=ff;
    type1=ff1(strncmpi([ff1.type],atomtype1,2))
    if numel(type1)==0
        type1=ff(strncmpi([ff1.type],atomtype1,1))
    end
    type2=ff2(strcmpi([ff2.type],atomtype2))
    if numel(type2)==0
        type2=ff2(strncmpi([ff2.type],atomtype2,2))
    end
    if numel(type2)==0
        type2=ff2(strncmpi([ff2.type],atomtype2,1))
    end
elseif size(ff,2)==2
    ff1=ff(1);
    ff2=ff(2);
    ff1=ff1.ff;
    ff2=ff2.ff;
    assignin('base','ff1',ff1)
    assignin('base','ff2',ff2)
    assignin('base','atomtype1',atomtype1)
    assignin('base','atomtype2',atomtype2)
    type1=ff1(strncmpi([ff1.type],atomtype1,2));
    if numel(type1)==0
        type1=ff1(strncmpi([ff1.type],atomtype1,1));
    end
    type2=ff2(strcmpi([ff2.type],atomtype2));
    if numel(type2)==0
        type2=ff2(strncmpi([ff2.type],atomtype2,2));
    end
    if numel(type2)==0
        type2=ff2(strncmpi([ff2.type],atomtype2,1));
    end
end

atomtype1=type1;
atomtype2=type2;

try
    if strncmpi([atomtype2.type],'Ow',2)
        if [atomtype2.charge]==0
            [atomtype2.charge]=-0.85
        end
    end
catch
    
end

try
    if strncmpi([atomtype2.type],'Hw',2)
        if [atomtype2.charge]==0
            [atomtype2.charge]=+0.425
        end
    end
catch
    
end

q1=atomtype1.charge
q2=atomtype2.charge

sig1=atomtype1.sigma_nm;
sig2=atomtype2.sigma_nm;
eps1=atomtype1.e_kJmol;
eps2=atomtype2.e_kJmol;


if isfield(ff1,'C4_kJmolnm')
    C4=atomtype1.C4_kJmolnm
    [q1 q2 sig1 sig2 eps1 eps2 C4]
else
    C4=[];
    [q1 q2 sig1 sig2 eps1 eps2]
end

r=.01:.0005:1.2; % nm

e_mix=(eps1*eps2)^.5
sig_mix=(sig1+sig2)/2


if numel(C4)==0
    
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
    
    lj=C12_mix./r.^12-C6_mix./r.^6-C4./r.^4;
end

coul=(1.60217646E-19)^2*6.022E+23*q1*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;
Utot=lj+coul;

if nargin>2
    %% Plot the Total energy, electrostatic contribution, and the LJ
    
    if nargin>3
        color1=varargin{2};
        color2=varargin{2};
        color3=varargin{2};
    else
        color1='b';
        color2='r';
        color3='k';
    end
    
    hold on
    plotmin=1000*(round(min(Utot*1.5)/1000));
    if plotmin>=0
        plotmin=10000;
    end
    
    if numel(C4)==0
        plot(r,lj,color1,'LineWidth',1);
        plot(r,coul,color2,'LineWidth',1);
        plot(r,Utot,color3,'LineWidth',1);
    else
        plot(r,lj,color1,'LineWidth',1);
        plot(r,coul,color2,'LineWidth',1);
        plot(r,Utot,color3,'LineWidth',1);
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


