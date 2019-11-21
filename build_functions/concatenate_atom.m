%% concatenate_atom.m 
% * This old function concatenats atom sections. Use update_atom({atom_1 atom_2}) instead...
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = concatenate_atom(atom_1,atom_2)
%
function atom = concatenate_atom(atom_1,atom_2)

% atom=[atom_1 atom_2];

atom=update_atom({atom_1 atom_2});