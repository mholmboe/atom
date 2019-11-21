%% duplicate_atom.m
% * This function duplicates residue with molid MolID
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = duplicate_atom(atom,4)
% # atom = duplicate_atom(atom,[1 2 3])
%
function atom = duplicate_atom(atom,molID) 

new_atom=atom(ismember([atom.molid],molID));

[new_atom.molid]=deal(molID-1);

last_molID_ind=max(find(ismember([atom.molid],molID)));
atom_top=atom(1:last_molID_ind);
atom_bottom=atom(last_molID_ind+1:end);
atom=[atom_top new_atom atom_bottom];
atom=update_atom(atom);
