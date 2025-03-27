%% eval_sim.m
% * This special function invokes lsqnolin, used in the forcefield optimization
% scheme called autofit_ff. This function is called by opt_ff.m
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # fx = run_opt_ff_lsqnonlin(x0,delta,dirtype,indexes); % f1,f2,f3); % etc to run with extra parameters..

function [fx,resnorm,residual,exitflag,output] = run_opt_ff_lsqnonlin(x0,delta,dirtype,varargin)

fextra=nargin-3;
if fextra>0
    for i=1:fextra
        eval(strcat('f',num2str(i),'=varargin{i};'));
    end
end

%% In order to keep parameters around 1..
scalefactors=(1./x0(1,:));
scalefactors(scalefactors==0)=1;
x0=scalefactors.*x0;

%% Initial values
disp('Initial parameters and lower/upper bounds!')
if numel(delta)>1
    x0=x0(1,:)
    lb=x0.*delta
    ub=x0./delta
else
    x0=x0(1,:)
    lb=min([(1-delta)*x0;(1+delta)*x0])
    ub=max([(1-delta)*x0;(1+delta)*x0])
end

%% Set Options for Optimization
% custOutput = @(x,optimValues,state,varargin)plot_iterations(x,optimValues,state,varargin);
options = optimoptions(@lsqnonlin,'Diagnostics','off','Display','iter','MaxFunctionEvaluations',5000,'MaxIterations',500,'Tolfun',1e-7); %,'PlotFcns',@plot_iterations);
% options = optimset('Diagnostics','on','Display','iter','MaxFunEvals',5000,'MaxIter',500,'Tolfun',1e-7);%,'OutputFcn',custOutput);
options.Algorithm = 'trust-region-reflective';%'levenberg-marquardt'; % or
options.FiniteDifferenceType = 'central'; %'forward'; %'central'; % forward
if ~mod(numel(x0),2)
    options.FiniteDifferenceStepSize = x0./50; % /10;%
    options.FiniteDifferenceStepSize(delta==1)=0.02; %0.05
else
    options.FiniteDifferenceStepSize = 0.02; % 0.05;%
end

%% Alternatively, use these settings
% options.DiffMinChange = 0.05;
% options.DiffMaxChange = Inf;


%% Run lsqnonlin
if fextra>0
    fstring=[];
    for i=1:fextra
        fstring=[fstring strcat(',f',num2str(i))];
    end
    eval(strcat('f=@(x)eval_sim(x,scalefactors,dirtype,indexes',fstring,')'));
else
    f=@(x)eval_sim(x,scalefactors,dirtype);
end
f
x0
lb
ub
options
[fx,resnorm,residual,exitflag,output]=lsqnonlin(f,x0,lb,ub,options);

%% Pass on the scale factors
assignin('caller','scalefactors',scalefactors)

end

