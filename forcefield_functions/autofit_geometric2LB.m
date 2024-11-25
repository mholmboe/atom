% clear;
close all;
hold on
format short;

r=.2:.001:1.2; % nm

%% Ref c6/c12
q1ref=0;
q2ref=0;
c61=2.26195E-03;
c62=0.0026173456; % Ow
c121=1.2769E-06;
c122=2.634129e-06; % Ow

% Ow 0.0026173456  2.634129e-06
% OP 2.26195E-03	1.77956E-06
% OO 3.59595E-04	1.11082E-07
% O  0.002261954    1.2769E-06

%% Initial guess sig/eps
q1=0;
q2=0;
sig1=0.31656;
sig2=0.31656;
eps1=0.718772678;
eps2=0.65017;

[r,lj,coul,data] = ljcoul_C6C12([q1ref,q2ref,c61,c62,c121,c122],r);

%% Initial values
xinit    = [ q1  q2  sig1  sig2  eps1  eps2 ]
delta    = [ 1   1   .1    1     1     1   ];

x0=xinit;
%% In order to keep parameters around 1..
scalefactors=1./x0;
scalefactors(scalefactors==0)=1;
scalefactors(isinf(scalefactors))=1E18;
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

% copy(fx)

r_plot=.2:.001:1.2;
[fxr,fxlj,fxcoul,fxdata] = ljcoul(fx,r_plot);

% norm((fxdata-data_plot)/numel(data_plot))
% plot(r_plot,fxdata-data_plot)

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
