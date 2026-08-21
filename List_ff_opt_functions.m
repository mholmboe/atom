%% List of forcefield optimization functions
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%

%% Forcefield optimization
% # <eval_cn.html eval_cn(rdf,cn)> % Calculate the coordination number from RDF and cumulative CN data, by finding the thickness of the first shell.
% # <eval_sim.html eval_sim(param,scalefactors,dirtype,varargin)> % Objective function for lsqnonlin used in the force field optimization scheme.
% # <eval_systems.html eval_systems(dirtype,varargin)> % Run and evaluate the reference simulation systems used in the force field optimization scheme.
% # <replace_row.html replace_row(old_row,new_row,filename,outfilename)> % Replace rows in text files, used when updating force field parameter files.
% # <replace_string.html replace_string(filename_in,filename_out,old_string,new_string)> % Replace strings in files.
% # <run_opt_ff_lsqnonlin.html run_opt_ff_lsqnonlin(x0,delta,dirtype,varargin)> % Invoke lsqnonlin to optimize force field parameters against reference simulation data.

%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
