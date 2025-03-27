% clear;
% close all;
figure
hold on
format short;

r=.1:.001:1.2; % nm

% q1=2.1;
% q2=-1.05;
% sig1=0.33;
% sig2=0.3165;
% eps1=7.7e-6;
% eps2=0.6504;
% C4=0.05;

hold on
ff=ff2.ff;
[rout,lj,coul,Utot,q1,q2,sig1,sig2,eps1,eps2,C4] = nonbonded_ff(ff,{'Fe2+' 'OW'});

[r,lj,coul,data] = ljcoul_C12C6C4([q1,q2,sig1,sig2,eps1,eps2,C4,C4],r);

% A1=41970.1381; % nm
% A2=41970.1381;
% B1=29.4118; % kJmol-1
% B2=29.4118;
% C1=0; % nm-1
% C2=0;
% [r,buck,coul,data] = buckinghamcoul([ q1  q2  A1 A2 B1 B2 C1 C2 ],r);

% [~, ind_start]=min(abs(Utot));
% rmin=rout(322)
% r=rmin:.001:1.2; % nm

[r,lj,coul,data] = ljcoul([q1,q2,sig1,sig2,eps1,eps2],r,1);

q1=1.20;
q2=-1.2;
sig1=.1;
sig2=.317;
eps1=.6;
eps2=.6837;
C4 = 0
%% Initial values
xinit    = [ q1  q2  sig1  sig2  eps1  eps2]
delta    = [ 1   1    .01     .01   .01    .01];

x0=xinit;
%% In order to keep parameters around 1..
scalefactors=1./x0;
scalefactors(scalefactors==0)=1;
x0=scalefactors.*x0;

%% Set lower/upper bounds
lb       = delta;
ub       = 1./delta;

ind_neg=find(lb>ub);
if numel(ind_neg)>0
    temp=lb(ind_neg);
    lb(ind_neg)=ub(ind_neg);
    ub(ind_neg)=temp;
end

%% Set Options for Optimization
objective_func_string='ljcoul_objective_func(x,data,r,scalefactors)';
eval(strcat('f=@(x)',objective_func_string,';'));
options = optimoptions(@lsqnonlin,'Diagnostics','on','Display','iter-detailed','MaxFunctionEvaluations',50000,'MaxIterations',5000,'Tolfun',1e-12);
options.Algorithm = 'levenberg-marquardt'; %'trust-region-reflective';%'levenberg-marquardt'; % or
options.StepTolerance = 1.000000000000000e-12
options.FiniteDifferenceType = 'central';
options.FiniteDifferenceStepSize = eps^(1/3);
options.CheckGradients = true;
options.OptimalityTolerance = 1E-18;

%% Run lsqnonlin
fx=lsqnonlin(f,x0,lb,ub,options);

fx=fx./scalefactors

copy(fx)

[fxr,fxlj,fxcoul,fxdata] = ljcoul(fx,r,1,'k');
% [r_plot,lj_plot,coul_plot,data_plot] = ljcoul_C12C6C4([q1,q2,sig1,sig2,eps1,eps2,C4,C4],r_plot);

% [fxr,fxbuck,fxcoul,fxdata] = buckinghamcoul([ q1  q2  A1 A2 B1 B2 C1 C2 ],r_plot);
% [r_plot,lj_plot,coul_plot,data_plot] = ljcoul(fx,r_plot);
% norm((fxdata-data_plot)/numel(data_plot))
% plot(r_plot,fxdata-data_plot)

x=fx;
eval(strcat('res=',objective_func_string,';'));
disp('Mean abs(residual):');
mean(abs(res))

% if nargin>5
% hold on
% plot(r,data,'r');
% plot(r,data+res,'g');
% plot(r,res,'b');
%     xlim([min(r) max(r)]);
%     ylim([min(data) max(data)]);
% end
