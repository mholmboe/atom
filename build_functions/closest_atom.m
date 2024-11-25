%% closest_atom.m
% * This function returns the atom1 struct with the nMolId's in atom1 closest
% to the atom2 struct. nMolId can be set to be lower than the original
% number of molecules in atom1.
%
%% Function arguments
% * atom1 is the normal atom struct
% * Box_dim is the normal Box_dim variable
% * atom1 is the atom struct to be sorted with respect to distance from atom2
% * atom2 is the substrate/system atom struct
% * nMolID is the number of molecules to keep
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom1 = closest_atom(atom1,atom2,Box_dim) % Basic input arguments
% # atom1 = closest_atom(atom1,atom2,Box_dim,nMolId) % Basic input arguments
%
function atom1 = closest_atom(atom1,atom2,Box_dim,varargin)

atom1=update_atom(atom1);

if nargin==3
    nMolId=max([atom1.molid]);
else
    nMolId=varargin{1};
end

dist_matrix = dist_matrix_atom(atom1,atom2,Box_dim);

dist_tot=[];
for i=1:max([atom1.molid])
    dist_matrix_temp=dist_matrix([atom1.molid]==i,:);
    dist_ave=median(median(dist_matrix_temp));
    dist_tot=[dist_tot dist_ave];
end
[dist_tot,indMolId]=sort(dist_tot);
indMolId=indMolId(1:nMolId);

ind=[];
for i=1:numel(indMolId)
    ind=[ind find([atom1.molid]==indMolId(i))];
end

atom1=atom1(ind);

atom1=update_atom(atom1);

end

