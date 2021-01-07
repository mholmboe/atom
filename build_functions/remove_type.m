%% remove_type.m
% * This function removes atomtypes with types as in typescell = 
% {'OW' 'HW1' 'HW2'}
%
%% Version
% 2.09
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = remove_type(atom,{'OW' 'HW1' 'HW2'})

function atom = remove_type(atom,typescell)

rm_ind=[];
for i=1:size(typescell,2)
    rm_ind=[rm_ind find(strcmp([atom.type],typescell{i}))];
end
rm_ind = sort(rm_ind);
nAtoms=size(atom,2);
change_ind=sort(setdiff([1:nAtoms],rm_ind));
atom=atom(change_ind);
atom=atom_update(atom);

assignin('caller','rm_ind',rm_ind)
