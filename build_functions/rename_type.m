%% rename_type.m
% * This function renames atoms in the atom struct
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=rename_type(atom,'NA','Na')
% # atom=rename_type(atom,'NA','Na',30)

function atom=rename_type(atom,atomtype,new_atomtype,varargin)

% Find all original atomtypes to replace
ind_atomtype=find(ismember([atom.type],atomtype)); 

% Duplicate the atomtype entries into a new In_atom struct
if nargin == 4
    % Uses the num last entries...
    num=cell2mat(varargin(1));
    ind_atomtype=ind_atomtype(1:num); 
end

% Rename
[atom(ind_atomtype).type]=deal({new_atomtype});
[atom(ind_atomtype).fftype]=deal({new_atomtype});

end