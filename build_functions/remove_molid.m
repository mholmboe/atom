%% remove_molid.m
% * This function removes molid MolID = [1 2 3 .....], updates the MolId's,
% and returns the new atom struct.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = remove_molid(atom,[2 3 10]) % Basic input arguments

function atom = remove_molid(atom,MolID)

atom(ismember([atom.molid],MolID))=[];
atom=update_atom(atom);
composition_atom(atom);

end

