clear;
close all;
hold on
format compact;

% initial_ff='clayff.mat'
ref_ff1='ions_Merz_12_6_4_divalent_opc3_ff.mat'
ref_ff2='water_models.mat'
extraff1='ions_Merz_CM_divalent_OPC3_ff'
extraff2='ions_Merz_IOD_divalent_OPC3_ff'
extraff3='ions_Merz_HFE_divalent_OPC3_ff'

Ion1='Mg'
Ion2='OW'
Water_model='opc3'

r=.1:.001:1.2; % nm
r=r-(r(2)-r(1))/2;
rfull=r;

hold on
% load(initial_ff)
% [rout,lj,coul,Utot,q1,q2,sig1,sig2,eps1,eps2] = nonbonded_ff(ff,{Ion1 Ion2});
% [~, ind_start]=min(abs(Utot));
% rmin=rout(ind_start+floor(ind_start/5))
% r=rmin:.001:1.2; % nm
% r=r-(r(2)-r(1))/2;

%[r,lj,coul,data] = ljcoul([q1,q2,sig1,sig2,eps1,eps2],r);

ff1=load(ref_ff1)
ff2=load(ref_ff2)
[rout,lj,coul,Utot,q1ref,q2ref,sig1ref,sig2ref,eps1ref,eps2ref,C4] = nonbonded_ff([ff1 ff2],{Ion1 strcat(Ion2,'_',Water_model)});
[~, ind_start]=min(abs(Utot));
rmin=rout(ind_start+floor(ind_start/5))
r=rmin:.001:1.2; % nm
r=r-(r(2)-r(1))/2;
[r,lj,coul,data] = ljcoul_force_C12C6C4([q1ref,q2ref,sig1ref,sig2ref,eps1ref,eps2ref,C4],r);

extraff1=load(extraff1);
[rt,ljt,coult,Ut,q1t,q2t,sig1t,sig2t,eps1t,eps2t] = nonbonded_ff([extraff1 ff2],{Ion1 strcat(Ion2,'_',Water_model)});
[rt,ljt,coult,datat] = ljcoul_force([q1t,q2t,sig1t,sig2t,eps1t,eps2t],rfull,1,'g');
extraff2=load(extraff2);
[rt,ljt,coult,Ut,q1t,q2t,sig1t,sig2t,eps1t,eps2t] = nonbonded_ff([extraff2 ff2],{Ion1 strcat(Ion2,'_',Water_model)});
[rt,ljt,coult,datat] = ljcoul_force([q1t,q2t,sig1t,sig2t,eps1t,eps2t],rfull,1,'m');
extraff3=load(extraff3);
[rt,ljt,coult,Ut,q1t,q2t,sig1t,sig2t,eps1t,eps2t] = nonbonded_ff([extraff3 ff2],{Ion1 strcat(Ion2,'_',Water_model)});
[rt,ljt,coult,datat] = ljcoul_force([q1t,q2t,sig1t,sig2t,eps1t,eps2t],rfull,1,'c');

% load(initial_ff)
% [rout,lj,coul,Utot,q1,q2,sig1,sig2,eps1,eps2] = nonbonded_ff(ff,{Ion1 Ion2});
% [~, ind_start]=min(abs(Utot));
% rmin=rout(ind_start+floor(ind_start/2.5))
% r=rmin:.001:1.2; % nm
% r=r-(r(2)-r(1))/2;

dvalue=2;
sig11=sig1ref/dvalue;
sig12=sig1ref*dvalue;
eps11=eps1ref/dvalue;
eps12=eps1ref*dvalue;
q11=q1ref/2;
q12=q1ref/2;

sig1=sig1ref;
sig2=sig2ref;
eps1=eps1ref;
eps2=eps2ref;
q2=q2ref;

%% Initial values
xinit    = [ q11  q12  q2   sig11  sig12   sig2  eps11  eps12  eps2]
delta    = [ 1    1    1    1      .01     1     1      .01     1   ];

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
objective_func_string='ljcoul_2x_force_objective_func(x,data,r,scalefactors)';
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

[fxr,fxlj,fxcoul,fxdata] = ljcoul_2x_force(fx,rfull);
% [r_plot,lj_plot,coul_plot,data_plot] = ljcoul([q1,q2,sig1,sig2,eps1,eps2],r_plot);

ff1=load(ref_ff1)
ff2=load(ref_ff2)
[rout,lj,coul,Utot,q1,q2,sig1,sig2,eps1,eps2,C4] = nonbonded_ff([ff1 ff2],{Ion1 strcat(Ion2,'_',Water_model)});
[r_plot,lj_plot,coul_plot,data_plot] = ljcoul_force_C12C6C4([q1,q2,sig1,sig2,eps1,eps2,C4,C4],rfull);

norm((fxdata-data_plot)/numel(data_plot))
plot(r_plot(2:end),fxdata-data_plot,'g--')
x=fx
eval(strcat('res=',objective_func_string,';'));
disp('Mean abs(residual):');
mean(abs(res))

plotmin=1000*(round2dec(min(data_plot*1.5)/1000));
if plotmin>=0
    plotmin=1000;
end
ylim(sort([-plotmin plotmin]))
