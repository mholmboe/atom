%% keep_atom.m
% * This function removes all but resname
%
%% Version
% 2.0
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = keep_atom(atom,{'SOL'})
%
function atom = keep_atom(atom,resname)

ind=ismember([atom.resname],resname);

atom=atom(ind);

atom=update_atom(atom);

