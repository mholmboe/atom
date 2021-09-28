% clear;
close all;
hold on
format compact;

r=.1:.001:1.2; % nm

C4 = 0.05;

hold on

[rout,lj,coul,Utot,q1,q2,sig1,sig2,eps1,eps2] = nonbonded_ff(ff,{'Si' 'Ob'});
%[r,lj,coul,data] = ljcoul([q1,q2,sig1,sig2,eps1,eps2],r);
[r,lj,coul,data] = ljcoul_C12C6C4([q1,q2,sig1,sig2,eps1,eps2,C4,C4],r);

dvalue=1.5;
eps11=eps1*dvalue;
eps12=eps1/dvalue;
sig11=sig1*dvalue;
sig12=sig1/dvalue;
q11=2.1/2;
q12=2.1/2;

%% Initial values
xinit    = [ q11  q12  q2  sig11  sig12   sig2  eps11  eps12   eps2]
delta    = [ 1    1    1    .1     .1     1     .1       .1    1   ];

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
objective_func_string='ljcoul_2x_objective_func(x,data,r,scalefactors)';
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

r_plot=.1:.001:1.2;


copy(fx)

[fxr,fxlj,fxcoul,fxdata] = ljcoul_2x(fx,r_plot);
% [r_plot,lj_plot,coul_plot,data_plot] = ljcoul([q1,q2,sig1,sig2,eps1,eps2],r_plot);
[r_plot,lj_plot,coul_plot,data_plot] = ljcoul_C12C6C4([q1,q2,sig1,sig2,eps1,eps2,C4,C4],r_plot);
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
