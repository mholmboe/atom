%% remove_molid.m
% * This function removes molid MolID = [1 2 3 .....], updates the MolId's,
% and returns the new atom struct.
%
%% Version
% 2.09
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = remove_molid(atom,[2 3 10])

function atom = remove_molid(atom,MolID)

atom(ismember([atom.molid],MolID))=[];
atom=update_atom(atom);
composition_atom(atom);

