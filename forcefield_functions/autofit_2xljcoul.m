clear;
close all;
hold on
format compact;

initial_ff='ions_Merz_IOD_monovalent_opc3_ff';
ref_ff1='ions_Merz_12_6_4_monovalent_opc3_ff.mat'
ref_ff2='water_models.mat'
ff=load(initial_ff);
ff1=load(ref_ff1)
ff2=load(ref_ff2)

% extraff1='ions_Merz_CM_monovalent_OPC3_ff';extraff1=load(extraff1);
% extraff2='ions_Merz_IOD_monovalent_OPC3_ff';extraff2=load(extraff2);
% extraff3='ions_Merz_HFE_monovalent_OPC3_ff';extraff3=load(extraff3);

Ion1='F-'
Ion2='OW'
Water_model='opc3'

r=.1:.001:1.2; % nm
rfull=r;

hold on

% nonbonded_ff([extraff1 ff2],{Ion1 strcat(Ion2,'_',Water_model)},1,'g');
% nonbonded_ff([extraff2 ff2],{Ion1 strcat(Ion2,'_',Water_model)},1,'m');
% nonbonded_ff([extraff3 ff2],{Ion1 strcat(Ion2,'_',Water_model)},1,'c');

[rout,lj_init,coul_init,Utot_init,q1_init,q2_init,sig1_init,sig2_init,eps1_init,eps2_init] = nonbonded_ff([ff ff2],{Ion1 strcat(Ion2,'_',Water_model)});
[rout,lj,coul,Utot,q1ref,q2ref,sig1ref,sig2ref,eps1ref,eps2ref,C4] = nonbonded_ff([ff1 ff2],{Ion1 strcat(Ion2,'_',Water_model)});
[~, ind_start]=min(abs(Utot));
rmin=r(1)%rout(ind_start)%+floor(ind_start/20))
r=rmin:.001:1.2; % nm
[r,lj,coul,data] = ljcoul_C12C6C4([q1ref,q2ref,sig1ref,sig2ref,eps1ref,eps2ref,C4,C4],r);

dvalue=1;
sig11=sig1_init/dvalue;
sig12=sig1_init*dvalue;
eps11=eps1_init/dvalue;
eps12=eps1_init*dvalue;
q11=q1_init;
q12=1E-7;

sig1=sig1ref;
sig2=sig2ref;
eps1=eps1ref;
eps2=eps2ref;
q2=q2ref;

%% Initial values
xinit    = [ q11  q12  q2   sig11  sig12   sig2  eps11  eps12  eps2]
delta    = [ 1    1    1    .5     .1      1     .5     .1     1   ];

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

% copy(fx)

[fxr,fxlj,fxcoul,fxdata] = ljcoul_2x(fx,rfull);
% [r_plot,lj_plot,coul_plot,data_plot] = ljcoul([q1,q2,sig1,sig2,eps1,eps2],r_plot);

ff1=load(ref_ff1);
ff2=load(ref_ff2);
[rout,lj,coul,Utot,q1,q2,sig1,sig2,eps1,eps2,C4] = nonbonded_ff([ff1 ff2],{Ion1 strcat(Ion2,'_',Water_model)});
[r_plot,lj_plot,coul_plot,data_plot] = ljcoul_C12C6C4([q1,q2,sig1,sig2,eps1,eps2,C4,C4],rfull);

disp('Init values')
[q1_init,q2_init,sig1_init,sig2_init,eps1_init,eps2_init]

norm((fxdata-data_plot)/numel(data_plot))
plot(r_plot,fxdata-data_plot,'g--')
disp('Final values')
x=fx
eval(strcat('res=',objective_func_string,';'));
disp('Mean abs(residual):');
mean(abs(res))

Err=100*mean(abs(fxdata(ind_start:end)-data_plot(ind_start:end)))

plotmin=1000*(round2dec(min(Utot*1.5)/1000));
if plotmin>=0
    plotmin=1500;
end
ylim(sort([-plotmin plotmin]))
