%% neighbor_atom.m
% * This function scans the atom1 struct coordinates (containing a single site)
% * within a certain radius. It outputs neighbour index, distance and
% * coordinates of the neighbours. Triclinic support untested
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # neigh = neighbor_atom(atom1,atom2,Box_dim) % Basic input arguments
% # neigh = neighbor_atom(atom1,atom2,Box_dim,2.25) % Set the neighbour cutoff radius manually
%
function neigh = neighbor_atom(atom1,atom2,Box_dim,varargin)

if nargin < 4
    radius=2.25;
else
    radius=varargin{1};
end
if isstruct(atom1)
    XYZ_solute=single([[atom1.x]' [atom1.y]' [atom1.z]']); % use of single instead of double
    XYZ_data=single([[atom2.x]' [atom2.y]' [atom2.z]']); % use of single instead of double
else
    XYZ_solute=single(atom1); % use of single instead of double
    XYZ_data=single(atom2); % use of single instead of double
end

% Box_dim=[lx ly lz 0 0 xy 0 xz yz];
lx=Box_dim(1);ly=Box_dim(2);lz=Box_dim(3);
if numel(Box_dim)==3
    xy=0;xz=0;yz=0;
else
    xy=Box_dim(6);xz=Box_dim(8);yz=Box_dim(9);
end

% Calculate Distance Components %
rx = XYZ_data(:,1) - XYZ_solute(1);
ry = XYZ_data(:,2) - XYZ_solute(2);
rz = XYZ_data(:,3) - XYZ_solute(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_gt_ind=find(rz > lz/2);z_lt_ind=find(rz < - lz/2);
rz(z_gt_ind) = rz(z_gt_ind) - lz;
rz(z_lt_ind) = rz(z_lt_ind) + lz;
rx(z_gt_ind) = rx(z_gt_ind) - xz;
rx(z_lt_ind) = rx(z_lt_ind) + xz;
ry(z_gt_ind) = ry(z_gt_ind) - yz;
ry(z_lt_ind) = ry(z_lt_ind) + yz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_gt_ind=find(ry > ly/2);y_lt_ind=find(ry < - ly/2);
ry(y_gt_ind) = ry(y_gt_ind) - ly;
ry(y_lt_ind) = ry(y_lt_ind) + ly;
rx(y_gt_ind) = rx(y_gt_ind) - xy;
rx(y_lt_ind) = rx(y_lt_ind) + xy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_gt_ind=find(rx > lx/2);x_lt_ind=find(rx < - lx/2);
rx(x_gt_ind) = rx(x_gt_ind) - lx;
rx(x_lt_ind) = rx(x_lt_ind) + lx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist = sqrt( rx(:,1).^2 + ry(:,1).^2 + rz(:,1).^2 ); % distance calc.

% solute_index=find(dist==0);
% Find points inside circle
in=intersect(find(dist>0),find(dist<radius));
neigh.in = in;
neigh.dist = [dist(in,1)];
neigh.coords = [XYZ_data(in,1) XYZ_data(in,2) XYZ_data(in,3)];
neigh.r_vec = [rx(in) ry(in) rz(in)];


