%% remove_occypancy_atom.m
% * This function removes all succeding particles in the atom struct that has
% * identical coordinates to a preceding particle
%
%% Similar
% fuse_atom
%
%% Version
% 3.0.2
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * atom=remove_occypancy_atom(atom) % Basic input arguments
%
%% See also
% # <fuse_atom.html fuse_atom>
% # <remove_SOL.html remove_SOL>
% # <remove_molid.html remove_molid>
% # <remove_residues.html remove_residues>
% # <remove_resname.html remove_resname>
%
function atom=remove_occypancy_atom(atom)

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];

[u,I,J] = unique(XYZ_data, 'rows', 'first');
hasDuplicates = size(u,1) < size(XYZ_data,1);
ixDupRows = setdiff(1:size(XYZ_data,1), I);
dupRowValues = XYZ_data(ixDupRows,:);
atom(ixDupRows)=[];

end