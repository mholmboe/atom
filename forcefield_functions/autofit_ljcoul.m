% clear;
close all;
format compact;

r=.1:.001:1.2; % nm

% q1=2.1;
% q2=-1.05;
% sig1=0.33;
% sig2=0.3165;
% eps1=7.7e-6;
% eps2=0.6504;

hold on
[rout,lj,coul,Utot,q1,q2,sig1,sig2,eps1,eps2] = nonbonded_ff(ff,{'Lio' 'Ob'});
[r,lj,coul,data] = ljcoul(q1,q2,sig1,sig2,eps1,eps2,r);

eps1=.6;
eps2=.2;
q1=.525;
q2=-1.05;
%% Initial values
xinit    = [ q1  q2  sig1  sig2  eps1  eps2  ];
delta    = [ 1   1   .2    1     1     1];

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
options = optimset('Diagnostics','on','Display','iter','MaxFunEvals',5000,...
    'MaxIter',5000,'Tolfun',1e-14);%,'OutputFcn',custOutput);

% custOutput = @(x,optimValues,state,varargin)plot_iterations(x,optimValues,state,varargin);
% options = optimoptions(@lsqnonlin,'Diagnostics','off','Display','off','MaxFunctionEvaluations',5000,'MaxIterations',500,'Tolfun',1e-7,'PlotFcns',@plot_iterations);
% options.Algorithm = 'trust-region-reflective';%'levenberg-marquardt'; % or
% options.FiniteDifferenceStepSize = 'sqrt(eps)'; % 0.5;%
% options.DiffMinChange = 1;
% options.DiffMaxChange = 0.1;

%% Run lsqnonlin
fx=lsqnonlin(f,x0,lb,ub,options);

fx=fx./scalefactors

copy(fx)

[fxr,fxlj,fxcoul,fxdata] = ljcoul(fx(1),fx(2),fx(3),fx(4),fx(5),fx(6),r);

plot(r,data-fxdata)
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
