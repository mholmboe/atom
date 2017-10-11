%% dist_matrix.m
% * This function calculates the distance matrix from the atom struct, or the distances between atom1 and atom2
% * Note that there is also and dist_matrix_atomv2 function that might be faster for large systems...
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se


%% Examples
% * dist_matrix = dist_matrix_atom(atom1,Box_dim)
% * dist_matrix = dist_matrix_atom(atom1,atom2,Box_dim)

function dist_matrix = dist_matrix_atom(atom1,varargin) % ,atom2,Box_dim); % or % ,Box_dim);

format compact;

if nargin==3
    atom2=varargin{1};
    Box_dim=varargin{2};
else
    atom2=atom1;
    Box_dim=varargin{1};
end

nAtoms1=size(atom1,2);
nAtoms2=size(atom2,2);

XYZ1=[[atom1.x]' [atom1.y]' [atom1.z]'];
XYZ2=[[atom2.x]' [atom2.y]' [atom2.z]'];

lx=Box_dim(1);ly=Box_dim(2);lz=Box_dim(3);
if length(Box_dim) > 3
    xy=Box_dim(6);xz=Box_dim(8);yz=Box_dim(9);
else
    xy=0;xz=0;yz=0;
end

dist_matrix = zeros(nAtoms2,nAtoms1);
X_dist = dist_matrix;
Y_dist = dist_matrix;
Z_dist = dist_matrix;
if size(Box_dim,2)==3
    for i = 1:size(XYZ1,1)
        %Calculate Distance Components
        rx = XYZ1(i,1) - XYZ2(:,1);
        x_gt_ind=find(rx > lx/2); x_lt_ind=find(rx < - lx/2);
        rx(x_gt_ind) = rx(x_gt_ind) - lx;
        rx(x_lt_ind) = rx(x_lt_ind) + lx;
        
        ry = XYZ1(i,2) - XYZ2(:,2);
        y_gt_ind=find(ry > ly/2); y_lt_ind=find(ry < - ly/2);
        ry(y_gt_ind) = ry(y_gt_ind) - ly;
        ry(y_lt_ind) = ry(y_lt_ind) + ly;
        
        rz = XYZ1(i,3) - XYZ2(:,3);
        z_gt_ind=find(rz > lz/2); z_lt_ind=find(rz < - lz/2);
        rz(z_gt_ind) = rz(z_gt_ind) - lz;
        rz(z_lt_ind) = rz(z_lt_ind) + lz;
        
        r = sqrt( rx(:,1).^2 + ry(:,1).^2 + rz(:,1).^2 ); % distance calc.
        dist_matrix(:,i)=r;
        X_dist(:,i)=rx;
        Y_dist(:,i)=ry;
        Z_dist(:,i)=rz;
        
        if mod(i,1000)==1
            if i > 1
                i-1
            end
        end
    end
else % and this is not actually tested yet...
    for i = 1:size(XYZ1,1)
        %Calculate Distance Components
        rx = XYZ1(i,1) - XYZ2(:,1);
        x_gt_ind=find(rx > lx/2); x_lt_ind=find(rx < - lx/2);
        rx(x_gt_ind) = rx(x_gt_ind) - lx;
        rx(x_lt_ind) = rx(x_lt_ind) + lx;
        
        ry = XYZ1(i,2) - XYZ2(:,2);
        y_gt_ind=find(ry > ly/2); y_lt_ind=find(ry < - ly/2);
        ry(y_gt_ind) = ry(y_gt_ind) - ly;
        ry(y_lt_ind) = ry(y_lt_ind) + ly;
        rx(y_gt_ind) = rx(y_gt_ind) - xy;
        rx(y_lt_ind) = rx(y_lt_ind) + xy;
        
        rz = XYZ1(i,3) - XYZ2(:,3);
        z_gt_ind=find(rz > lz/2); z_lt_ind=find(rz < - lz/2);
        rz(z_gt_ind) = rz(z_gt_ind) - lz;
        rz(z_lt_ind) = rz(z_lt_ind) + lz;
        rx(z_gt_ind) = rx(z_gt_ind) - xz;
        rx(z_lt_ind) = rx(z_lt_ind) + xz;
        ry(z_gt_ind) = ry(z_gt_ind) - yz;
        ry(z_lt_ind) = ry(z_lt_ind) + yz;
        
        r = sqrt( rx(:,1).^2 + ry(:,1).^2 + rz(:,1).^2 ); % distance calc.
        dist_matrix(:,i)=r;
        X_dist(:,i)=rx;
        Y_dist(:,i)=ry;
        Z_dist(:,i)=rz;
        
        if mod(i,1000)==1
            if i > 1
                i-1
            end
        end
    end
end


% New transposed output
dist_matrix=dist_matrix';
assignin('caller','X_dist',(X_dist)');
assignin('caller','Y_dist',(Y_dist)');
assignin('caller','Z_dist',(Z_dist)');
assignin('caller','analyzed_Box_dim',Box_dim);

