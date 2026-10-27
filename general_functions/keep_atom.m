%% keep_atom.m
% * This function removes all but resname
%
%% Version
% 3.0.1
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = keep_atom(atom,{'SOL'})
%
%% See also
% # <update_atom.html update_atom>
% # <keep_resname.html keep_resname>
%
function atom = keep_atom(atom,resname)

ind=ismember([atom.resname],resname);

atom=atom(ind);

atom=update_atom(atom);

