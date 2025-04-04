% clear;
close all;
hold on
format compact;

r=.12:.001:1.2; % nm

C4=0.05;

hold on
[rout,lj,coul,Utot,q1,q2,sig1,sig2,eps1,eps2] = nonbonded_ff(ff,{'Cc' 'Ob'});

[r,lj,coul,data] = ljcoul_C12C6C4([q1,q2,sig1,sig2,eps1,eps2,C4,C4],r);
%[r,lj,coul,data] = ljcoul([q1,q2,sig1,sig2,eps1,eps2],r);

% eps1=.6;
% eps2=.6;
% q1=2.1;
% q2=-1.05;
C4 = 0.05
%% Initial values
xinit    = [ q1  q2  5E8 5E8 60 60 0.007 0.007 ]
delta    = [ 1   1   .01 .01 .01 .01 .01   .01 ];

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
objective_func_string='buckinghamcoul_objective_func(x,data,r,scalefactors)';
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

r_plot=.12:.001:1.2;
[fxr,fxlj,fxcoul,fxdata] = buckinghamcoul(fx,r_plot);
[r_plot,lj_plot,coul_plot,data_plot] = ljcoul_C12C6C4([q1,q2,sig1,sig2,eps1,eps2,C4,C4],r_plot);
%[r_plot,lj_plot,coul_plot,data_plot] = buckinghamcoul([q1,q2,sig1,sig2,eps1,eps2],r);
norm((fxdata-data_plot)/numel(data_plot))
plot(r_plot,fxdata-data_plot,'g--')

% x0=fx;
% eval(strcat('res=',objective_func_string,';'));
% disp('Mean abs(residual):');
% mean(abs(res))

% if nargin>5
% hold on
% plot(r,data,'r');
% plot(r,data+res,'g');
% plot(r,res,'b');
%     xlim([min(r) max(r)]);
%     ylim([min(data) max(data)]);
% end
