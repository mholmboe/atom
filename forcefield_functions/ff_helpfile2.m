%% ff_helpfile2.m
% * This helper script inspects and renames atomtypes in an already saved
% * forcefield (ff) struct.
%
%% Version
% 3.0.1
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # ff_helpfile2 % Edit ffname at the top of the script
%

%% See also
% # <ff_helpfile.html ff_helpfile>
% # <ff_helpfile3.html ff_helpfile3>
% # <ff_plotfile.html ff_plotfile>
% # <run_opt_ff_lsqnonlin.html run_opt_ff_lsqnonlin>
%
format long
clear all; clc; close all;
ffname='ions_Merz_IOD_monovalent_opc_ff.mat'
load(ffname)
for n=1:size(ff,2)
    if strncmpi([ff(n).type],'F',1)
        [ff.type],'F-'
%     elseif
%    
    end     
end

% save(strcat('new_',ffname),'ff')