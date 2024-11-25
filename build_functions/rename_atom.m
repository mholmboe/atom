%% rename_atom.m
% * This function renames atoms in the atom struct
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # new_atom=rename_atom(atom,'Al','Mgo') % Basic input arguments
%
function new_atom=rename_atom(atom,old_atomtype,new_atomtype)

new_atom=atom;
% Find all original atomtypes to replace
ind_atomtype=find(strcmp([atom.type],old_atomtype));

% Rename all
[new_atom(ind_atomtype).type]=deal({new_atomtype});
[new_atom(ind_atomtype).fftype]=deal({new_atomtype});

end