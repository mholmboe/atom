%% overwrite_atom.m
% * This function overwrites the atom struct information with new information 
%
%% Version
% 2.082
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # Out_atom = overwrite_atom(atom,'Na','MOL')
% 
function Out_atom=overwrite_atom(In_atom,atomtype,resname)

index=num2cell(1:size(In_atom,2));
[In_atom.index]=deal(index{:});
[In_atom.resname]=deal({resname});
[In_atom.type]=deal({atomtype});
[In_atom.fftype]=deal({atomtype});

Out_atom=In_atom;
