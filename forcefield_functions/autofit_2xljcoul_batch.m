clear all;clc;close all;

Resultfile='monovalent_pseudo_1264_IOD_SPCE_v1'
initial_ff='ions_Merz_IOD_monovalent_all_ff'; %'ions_Merz_12_6_4_monovalent_opc3_ff.mat' 
ref_ff1='ions_Merz_12_6_4_monovalent_spce_ff.mat'
ref_ff2='water_models.mat'
ff=load(initial_ff);
ff1=load(ref_ff1)
ff2=load(ref_ff2)

% extraff1='ions_Merz_CM_monovalent_OPC3_ff';
% extraff2='ions_Merz_IOD_monovalent_OPC3_ff';
% extraff3='ions_Merz_HFE_monovalent_OPC3_ff';

% Ion1='Ca'
Ion2='OW'
Water_model='spce'

% %% Initial values
%xinit    = [ q11  q12  q2   sig11  sig12   sig2  eps11  eps12  eps2]
delta    =  [ 1    1    1    1     .001      1     1     .0001     1   ];

ff=ff.ff;

Atom_labels=[ff.type];

Results=[];Error=[];
for i=1:numel(Atom_labels)
    figure
    i
    Atom_labels{i}
    [x,Err]=autofit_2xljcoul_func(Atom_labels{i},Ion2,Water_model,initial_ff,ref_ff1,ref_ff2,delta);%,varargin)
    Results(i).type=Atom_labels(i);
    [radii,type_radii]=radius_ion(Atom_labels{i});
    Results(i).ionicradii=radii;
    Results(i).type_ionicradii=type_radii;
    Results(i).sigma_nm=x(4);
    Results(i).sigma2_nm=x(5);
    Results(i).e_kJmol=x(7);
    Results(i).e2_kJmol=x(8);
    Results(i).Err=Err;
    title(Atom_labels{i})
end

figure
bar(categorical([Results.type]),[Results.Err]);
title('Error')

figure
[sig1,sig1_order]=sort([Results.sigma_nm]);
temp=Results;
for i=1:size(Results,2)
    temp(i)=Results(sig1_order(i));
end
Results_sig=temp;
Xlabels=categorical([Results_sig.type]);
Xlabels=reordercats(Xlabels,[Results_sig.type]);
bar(Xlabels,[Results_sig.Err]);
title('sigma1 vs. Error')


figure
[e1,e1_order]=sort([Results.e_kJmol]);
temp=Results;
for i=1:size(Results,2)
    temp(i)=Results(e1_order(i));
end
Results_eps=temp;
Xlabels=categorical([Results_eps.type]);
Xlabels=reordercats(Xlabels,[Results_eps.type]);
bar(Xlabels,[Results_eps.Err]);
title('epsilon1 vs. Error')

figure
[r,r_order]=sort([Results.ionicradii]);
temp=Results;
for i=1:size(Results,2)
    temp(i)=Results(r_order(i));
end
Results_radii=temp;
Xlabels=categorical([Results_radii.type]);
Xlabels=reordercats(Xlabels,[Results_radii.type]);
bar(Xlabels,[Results_radii.Err]);
title('Radii vs. Error')

save(strcat(Resultfile,'.mat'),'Results')