%% remove_molid.m
% * This function removes molid MolID = [1 2 3 .....], updates the MolId's,
% and returns the new atom struct.
%
%% Version
% 3.0.2
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = remove_molid(atom,[2 3 10]) % Basic input arguments

%% See also
% # <composition_atom.html composition_atom>
% # <update_atom.html update_atom>
% # <molid_rotate.html molid_rotate>
% # <molid_translate.html molid_translate>
% # <position_molid.html position_molid>
%
function atom = remove_molid(atom,MolID)

atom(ismember([atom.molid],MolID))=[];
atom=update_atom(atom);
composition_atom(atom);

end

