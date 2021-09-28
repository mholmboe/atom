format compact

%% Merz 2014
% HFE
clear
P={};
PolyvalentIons=P(1,1:2:end)
R2min_tip3p=cell2mat(P(1,2:2:end))
eps_tip3p=cell2mat(P(2,:))

R2min_spce=cell2mat(P(3,:))
eps_spce=cell2mat(P(4,:))

R2min_tip4pew=cell2mat(P(5,1:2:end))
eps_tip4pew=cell2mat(P(5,2:2:end))

save('Merz_polyvalent_HFE_2014')

%% Merz 2014
% IOD
clear
P={};
PolyvalentIons=P(1,1:2:end)
R2min_tip3p=cell2mat(P(1,2:2:end))
eps_tip3p=cell2mat(P(2,:))

R2min_spce=cell2mat(P(3,1:2:end))
eps_spce=cell2mat(P(3,2:2:end))

R2min_tip4pew=cell2mat(P(4,:))
eps_tip4pew=cell2mat(P(5,:))

save('Merz_polyvalent_IOD_2014')

%% Merz 2014
% 12-6-4
clear
P={};
PolyvalentIons=P(1,1:2:48)
R2min_tip3p=cell2mat(P(1,2:2:end))
eps_tip3p=cell2mat(P(2,1:2:48))
C4_tip3p=cell2mat(P(2,2:2:48))

R2min_spce=cell2mat(P(3,1:24))
eps_spc=cell2mat(P(4,1:2:48))
C4_spce=cell2mat(P(4,2:2:48))

R2min_tip4pew=cell2mat(P(5,1:3:end))
eps_tip4pew=cell2mat(P(5,2:3:end))
C4_tip4pew=cell2mat(P(5,3:3:end))

save('Merz_polyvalent_12-6-4_2014')

%% Merz 2020
% HFE
clear
P={};
PolyvalentIons=[P(1,1:2:22) P(2,1:2:18) P(3,1:2:8)]
R2min_OPC3=cell2mat([P(1,2:2:22) P(2,2:2:18) P(3,2:2:8)])
eps_OPC3=cell2mat(P(4,1:24))

R2min_OPC=cell2mat(P(5,1:2:48))
eps_OPC=cell2mat(P(5,2:2:48))

R2min_tip3p_fb=cell2mat(P(6,1:2:end))
eps_tip3p_fb=cell2mat(P(6,2:2:end))

R2min_tip4p_fb=cell2mat(P(7,1:2:end))
eps_tip4p_fb=cell2mat(P(7,2:2:end))

save('Merz_Divalent_HFE_2020')

%%
%% Merz 2020
% IOD
clear
P={};
PolyvalentIons=[P(1,1:2:16) P(2,1:2:12) P(3,1:2:4)]
R2min_OPC3=cell2mat([P(1,2:2:16) P(2,2:2:12) P(3,2:2:4)])
eps_OPC3=cell2mat(P(4,1:16))

R2min_OPC=cell2mat(P(5,1:2:32))
eps_OPC=cell2mat(P(5,2:2:32))

R2min_tip3p_fb=cell2mat(P(6,1:2:end))
eps_tip3p_fb=cell2mat(P(6,2:2:end))

R2min_tip4p_fb=cell2mat(P(7,1:2:end))
eps_tip4p_fb=cell2mat(P(7,2:2:end))

save('Merz_Divalent_IOD_2020')

%% Merz 2020
% CM
clear
P={};
PolyvalentIons=[P(1,1:2:22) P(2,1:2:18) P(3,1:2:8)]
R2min_OPC3=cell2mat([P(1,2:2:22) P(2,2:2:18) P(3,2:2:8)])
eps_OPC3=cell2mat(P(4,1:24))

R2min_OPC=cell2mat(P(5,1:2:48))
eps_OPC=cell2mat(P(5,2:2:48))

R2min_tip3p_fb=cell2mat(P(6,1:2:end))
eps_tip3p_fb=cell2mat(P(6,2:2:end))

R2min_tip4p_fb=cell2mat(P(7,1:2:end))
eps_tip4p_fb=cell2mat(P(7,2:2:end))

save('Merz_Divalent_CM_2020')

%% Merz 2020
% 12-6-4
clear
P={};
PolyvalentIons=[P(1,1:7:56) P(2,1:7:28) P(4,1:13:52)]
R2min_OPC3=cell2mat([P(1,2:7:56) P(2,2:7:28) P(4,2:13:52)])
eps_OPC3=cell2mat([P(1,3:7:56) P(2,3:7:28) P(4,3:13:52)])
C4_OPC3=cell2mat([P(1,4:7:56) P(2,4:7:28) P(4,4:13:52)])

R2min_OPC=cell2mat([P(1,5:7:56) P(2,5:7:28) P(4,5:13:52)])
eps_OPC=cell2mat([P(1,6:7:56) P(2,6:7:28) P(4,6:13:52)])
C4_OPC=cell2mat([P(1,7:7:56) P(2,7:7:28) P(4,7:13:52)])

R2min_tip3p_fb=cell2mat([P(3,1:6:end) P(4,8:13:52)])
eps_tip3p_fb=cell2mat([P(3,2:6:end) P(4,9:13:52)])
C4_tip3p_fb=cell2mat([P(3,3:6:end) P(4,10:13:52)])

R2min_tip4p_fb=cell2mat([P(3,4:6:end) P(4,11:13:52)])
eps_tip4p_fb=cell2mat([P(3,5:6:end) P(4,12:13:52)])
C4_tip4p_fb=cell2mat([P(3,6:6:end) P(4,13:13:52)])

save('Merz_Divalent_12-6-4_2020')

%%
%% Merz 2015
% HFE
clear
P={};
PolyvalentIons=P(1,1:2:end)
R2min_tip3p=cell2mat(P(1,2:2:end))
eps_tip3p=cell2mat(P(2,:))

R2min_spce=cell2mat(P(3,:))
eps_spce=cell2mat(P(4,:))

R2min_tip4pew=cell2mat(P(5,1:2:end))
eps_tip4pew=cell2mat(P(5,2:2:end))

save('Merz_Monovalent_HFE_2015')

%% Merz 2015
% IOD
clear
P={};
PolyvalentIons=[P(1,1:2:16) P(2,1:2:12) P(3,1:2:4)]
R2min_OPC3=cell2mat([P(1,2:2:16) P(2,2:2:12) P(3,2:2:4)])
eps_OPC3=cell2mat(P(4,1:16))

R2min_OPC=cell2mat(P(5,1:2:32))
eps_OPC=cell2mat(P(5,2:2:32))

R2min_tip3p_fb=cell2mat(P(6,1:2:end))
eps_tip3p_fb=cell2mat(P(6,2:2:end))

R2min_tip4p_fb=cell2mat(P(7,1:2:end))
eps_tip4p_fb=cell2mat(P(7,2:2:end))

save('Merz_Monovalent_IOD_2015')

%% Merz 2015
% IOD
clear
P={};
PolyvalentIons_all=[P(1:16,1)]
R2min_all=cell2mat([P(17,1:16)])
eps_all=cell2mat(P(18,1:16))

save('Merz_Monovalent_IOD_2015')

%% Merz 2015
% 12-6-4
clear
P={};
PolyvalentIons_all=[P(1:16,1)]

R2min_tip3p=cell2mat(P(17,1:3:end))
eps_tip3p=cell2mat(P(17,2:3:end))
C4_tip3p=cell2mat(P(17,3:3:end))

R2min_spce=cell2mat(P(18,1:16))
eps_spce=cell2mat(P(19,1:2:32))
C4_spce=cell2mat(P(19,2:2:32))

R2min_tip4pew=cell2mat(P(20,1:3:end))
eps_tip4pew=cell2mat(P(20,2:3:end))
C4_tip4pew=cell2mat(P(20,3:3:end))

save('Merz_Monovalent_12-6-4_2015')

