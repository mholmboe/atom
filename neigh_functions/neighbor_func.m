%% neighbor_func.m
% * This function scans xyz data and checks who is within a certain radius 
% * It outputs neighbour index, distance and coordinates of the neighbours
%  * Triclinic support untested
%
%% Version
% 2.0
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # Neigh = neighbor_func([1 4 7 2],XYZ_solute,XYZ_data,Box_dim,2.25)
%
function Neigh = neighbor_func(solute_index,XYZ_solute,XYZ_data,Box_dim,radius)

% Box_dim=[lx ly lz 0 0 xy 0 xz yz];

rx=zeros(size(XYZ_data,1),1);ry=zeros(size(XYZ_data,1),1);rz=zeros(size(XYZ_data,1),1);

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
% Find points inside circle
in=intersect(find(dist>0),find(dist<radius));
% in=sort(unique(in));
in = in(find(in~=solute_index));
Neigh.in = [in(find(in~=solute_index))];
Neigh.dist = [dist(in,1)];
Neigh.coords = [XYZ_data(in,1) XYZ_data(in,2) XYZ_data(in,3)];
Neigh.r_vec = [rx(in) ry(in) rz(in)];


