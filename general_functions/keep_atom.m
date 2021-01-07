%% keep_atom.m
% * This function removes all but resname
%
%% Version
% 2.082
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = keep_atom(atom,{'SOL'})
%
function atom = keep_atom(atom,resname)

ind=ismember([atom.resname],resname);

atom=atom(ind);

atom=update_atom(atom);

