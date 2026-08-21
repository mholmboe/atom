%% ff_helpfile2.m
% * This helper script inspects and renames atomtypes in an already saved
% * forcefield (ff) struct.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # ff_helpfile2 % Edit ffname at the top of the script
%

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