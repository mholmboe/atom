%% remove_resname.m
% * This function removes residue with molid MolID and residue names
% resnames
%
%% Version
% 3.0.2
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = remove_resname(atom,{'SOL' 'Protein'}) % Basic input arguments
%
%% See also
% # <composition_atom.html composition_atom>
% # <remove_SOL.html remove_SOL>
% # <remove_molid.html remove_molid>
% # <remove_occypancy_atom.html remove_occypancy_atom>
% # <remove_residues.html remove_residues>
%
function atom = remove_resname(atom,resnames)

for i=1:size(resnames,2)
    atom(strcmp([atom.resname],resnames(i)))=[];
end

atom=atom_update(atom);
composition_atom(atom);

end