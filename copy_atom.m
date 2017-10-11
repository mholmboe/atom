%% copy_atom.m
% * This function copies and translates atoms in the atom struct
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * new_atom=copy_atom(atom,'Al','Mgo','ION',[0 0 10])
% * new_atom=copy_atom(atom,'Al','Mgo','ION',[0 0 10],10)


function new_atom=copy_atom(atom,atomtype,new_atomtype,new_resname,trans_vec,varargin)
%% 

% Find all original atomtypes to replace
ind_atomtype=find(strcmp([atom.type],atomtype)); 

if length(ind_atomtype) > 2
    ind_atomtype=ind_atomtype(randperm(length(ind_atomtype)));
end

% Duplicate the atomtype entries into a new In_atom struct
if nargin == 6
    % Uses the num last entries...
    num=cell2mat(varargin(1));
    new_atom=atom(ind_atomtype(end-num+1:end)); 
else
    % Uses all entries...
    new_atom=atom(ind_atomtype); 
end

% Rename all
index=num2cell(1:size(new_atom,2));
molid=num2cell(([1:size(new_atom,2)]+[atom(end).molid]));
[new_atom.molid]=deal(molid{:});
[new_atom.index]=deal(index{:});
[new_atom.resname]=deal({new_resname});
[new_atom.type]=deal({new_atomtype});
[new_atom.fftype]=deal({new_atomtype});

new_atom = translate_atom(new_atom,trans_vec,'all')

end