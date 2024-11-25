function [r,dispersion,coul,Utot] = dispcoul_C6C8C10(param,varargin)

C61=0;
C62=0;
C81=0;
C82=0;
C101=0;
C102=0;
Q1=0;
Q2=0;

if nargin>1
    ind=varargin{1};
    param=param(:,ind);
end

[nParams,nElements]=size(param);

if nElements>2
    disp('This function works best for only two elements');
    pause(1)
end

if nParams==4
%     C6=param(1,:);
%     C8=param(2,:);
%     C10=param(3,:);

    C61=param(1,1);
    C62=param(1,2);
    C81=param(2,1);
    C82=param(2,2);
    C101=param(3,1);
    C102=param(3,2);
end

if numel(find(param(end,:)<0))>0
    Q=param(end,:);
    if nParams>1
%         C6=param(1,:);
        C61=param(1,1);
        C62=param(1,2);
    else
        disp('Dispersion terms are missing!!!')
        pause
    end
    Q1=param(end,1);
    Q2=param(end,2);
end

if nargin > 2
    r=varargin{2};
else
    %     r=1:length(data);
    r=.01:.001:1.2; % nm
end

C6_mix=(C61*C62)^.5
C8_mix=(C81*C82)^.5;
C10_mix=(C101*C102)^.5;

dispersionC6=-C6_mix./r.^6;
dispersionC8=-C8_mix./r.^8;
dispersionC10=-C10_mix./r.^10;
dispersion=-C6_mix./r.^6-C8_mix./r.^8-C10_mix./r.^10;
coul=(1.60217646E-19)^2*6.022E+23*Q1*Q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;

% end

Utot=dispersion+coul;

hold on

plotmin=1000*(round(min(Utot*1.5)/1000));
if plotmin>=0
    plotmin=10000;
end
plot(r,dispersionC6,'r-');
plot(r,dispersionC8,'g-');
plot(r,dispersionC10,'b-');
plot(r,dispersion,'b--');
plot(r,coul,'r--');
plot(r,Utot,'k--');
xlabel('r [nm]');
ylabel('U [kJ/mol]');
xlim([0,1.2]);
ylim(sort([-plotmin plotmin]))
