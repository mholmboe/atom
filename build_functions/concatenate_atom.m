%% concatenate_atom.m 
% * This old function concatenats atom sections. Use update_atom({atom_1 atom_2}) instead...
%
%% Version
% 2.09
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = concatenate_atom(atom_1,atom_2) % Basic input arguments
%
function atom = concatenate_atom(atom_1,atom_2)
disp('This old function concatenats atom sections. Use update_atom({atom_1 atom_2}) instead...')

% atom=[atom_1 atom_2]; % This also works in principle, but does not update any indexes

atom=update_atom({atom_1 atom_2});