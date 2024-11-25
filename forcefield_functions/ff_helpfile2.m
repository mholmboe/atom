format long
clear all; clc; close all;
ffname='ions_Merz_IOD_monovalent_OPC_ff.mat'
load(ffname)
for n=1:size(ff,2)
    if strncmpi([ff(n).type],'F',1)
        [ff.type],'F-'
%     elseif
%    
    end     
end

% save(strcat('new_',ffname),'ff')