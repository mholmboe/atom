%% remove_occypancy_atom.m
% * This function removes all succeding particles in the atom struct that has
% * identical coordinates to a preceding particle
%
%% Similar
% fuse_atom
%
%% Version
% 2.09
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * atom=remove_occypancy_atom(atom)
%
function atom=remove_occypancy_atom(atom)

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];

[u,I,J] = unique(XYZ_data, 'rows', 'first');
hasDuplicates = size(u,1) < size(XYZ_data,1);
ixDupRows = setdiff(1:size(XYZ_data,1), I);
dupRowValues = XYZ_data(ixDupRows,:);
atom(ixDupRows)=[];
