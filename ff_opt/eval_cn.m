%% eval_cn.m
% * This special function calculates the CoordNum from RDF and cumulative
% CN data, by finding the thicknesss of the first coordination layer
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # CoordNum = eval_cn(rdf,cn)

function CoordNum = eval_cn(rdf,cn)

rdf((rdf(:,1)<0.125),2)=0; % To remove O-H 
cn((cn(:,1)<0.125),2)=0; % To remove O-H 

x=rdf(:,1);
y=rdf(:,2);
s=x(2)-x(1);
[IOD_p,IOD_max,w,p]=findpeaks(y,x,'SortStr','descend','NPeaks',2);
max_ind=ceil(min(IOD_max)/s);
[IOD_max_value,ind]=min(IOD_max);

dcn=smooth(diff(cn(max_ind:end,2)));
[minval,min_ind]=min(dcn);
CN_ind=max_ind+min_ind;
CoordNumDist=cn(CN_ind,1);
CoordNum=cn(CN_ind,2);

end
