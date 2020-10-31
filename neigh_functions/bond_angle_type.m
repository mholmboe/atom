%% bond_angle_type.m
% * This function tries to find all bonds and angles between the  atom1 and
% atom2 structures. One optional argument like 'reverse' reverses the 
% angle_limit from max to min angle
%
%% Version
% 2.08
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom1 = bond_angle_type(atom1,atom2,Box_dim,1.25,2.25,150)
% # atom1 = bond_angle_type(atom1,atom2,Box_dim,1.25,2.25,150,1)
% # atom1 = bond_angle_type(atom1,atom2,Box_dim,1.25,2.25,150,1,'reverse')

function atom1 = bond_angle_type(atom1,atom2,Box_dim,rmaxshort,rmaxlong,angle_limit,varargin)
%%

% rmin=0;
% rmax=3.5;
% min_angle=120;
% Box_dim=Box_dim(1,:);

if nargin > 6 % was 5, why?
    skip_internal=cell2mat(varargin(1));
else
    skip_internal=0;
end

lx=Box_dim(1);ly=Box_dim(2);lz=Box_dim(3);
if numel(Box_dim)==3
    xy=0;xz=0;yz=0;
elseif numel(Box_dim)==9
    % Box_dim=[lx ly lz 0 0 xy 0 xz yz];
    xy=Box_dim(6);xz=Box_dim(8);yz=Box_dim(9);
end

for i=1:size(atom1,2)
    XYZ_data=[[atom2.x]' [atom2.y]' [atom2.z]'];
    
    if sum(strncmpi([atom1.type],'OW',2))>0 && sum(strncmpi([atom2.type],'HW',2))>0 && skip_internal == 1
        ind_sel=~ismember([atom2.index],[atom1(i).index+1 atom1(i).index+2]);
        XYZ_data=XYZ_data(ind_sel,:);
    end
    
    solute_index=atom1(i).index;
    XYZ_solute=[[atom1(i).x]' [atom1(i).y]' [atom1(i).z]'];
    
    rx=zeros(size(XYZ_data,1),1);ry=zeros(size(XYZ_data,1),1);rz=zeros(size(XYZ_data,1),1);
    
    if size(Box_dim,2)>3
        % Calculate Distance Components for triclic cell with pbc (should work?)
        rx = XYZ_data(:,1) - XYZ_solute(1);
        ry = XYZ_data(:,2) - XYZ_solute(2);
        rz = XYZ_data(:,3) - XYZ_solute(3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        z_gt_ind=find(rz > lz/2);z_lt_ind=find(rz < - lz/2);
        rz(z_gt_ind) = rz(z_gt_ind) - lz;
        rz(z_lt_ind) = rz(z_lt_ind) + lz;
        rx(z_gt_ind) = rx(z_gt_ind) - xz;
        rx(z_lt_ind) = rx(z_lt_ind) + xz;
        ry(z_gt_ind) = ry(z_gt_ind) - yz;
        ry(z_lt_ind) = ry(z_lt_ind) + yz;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        y_gt_ind=find(ry > ly/2);y_lt_ind=find(ry < - ly/2);
        ry(y_gt_ind) = ry(y_gt_ind) - ly;
        ry(y_lt_ind) = ry(y_lt_ind) + ly;
        rx(y_gt_ind) = rx(y_gt_ind) - xy;
        rx(y_lt_ind) = rx(y_lt_ind) + xy;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x_gt_ind=find(rx > lx/2);x_lt_ind=find(rx < - lx/2);
        rx(x_gt_ind) = rx(x_gt_ind) - lx;
        rx(x_lt_ind) = rx(x_lt_ind) + lx;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        %Calculate Distance Components for orthogonal cell with pbc
        rx =  XYZ_data(:,1) - XYZ_solute(1);
        x_gt_ind=find(rx > lx/2);x_lt_ind=find(rx < - lx/2);
        rx(x_gt_ind) = rx(x_gt_ind) - lx;
        rx(x_lt_ind) = rx(x_lt_ind) + lx;
        
        ry = XYZ_data(:,2) - XYZ_solute(2);
        y_gt_ind=find(ry > ly/2);y_lt_ind=find(ry < - ly/2);
        ry(y_gt_ind) = ry(y_gt_ind) - ly;
        ry(y_lt_ind) = ry(y_lt_ind) + ly;
        
        rz = XYZ_data(:,3) - XYZ_solute(3);
        z_gt_ind=find(rz > lz/2);z_lt_ind=find(rz < - lz/2);
        rz(z_gt_ind) = rz(z_gt_ind) - lz;
        rz(z_lt_ind) = rz(z_lt_ind) + lz;
    end
    
    dist = sqrt( rx(:,1).^2 + ry(:,1).^2 + rz(:,1).^2 ); % distance calc.
    
    % Find points inside circle
    in=intersect(find(dist>rmaxshort),find(dist<rmaxlong));
    
    %in=sort(unique(in));
    
    
    %     in = in(find(in~=i));
    neigh.in = in;
    neigh.dist = [dist(in,1)];
    neigh.coords = [XYZ_data(in,1) XYZ_data(in,2) XYZ_data(in,3)];
    neigh.r_vec = [rx(in) ry(in) rz(in)];
    
    neigh_dist=[neigh.dist];
    neigh_ind=[neigh.in];
    neigh_vec=[neigh.r_vec];
    
    neigh_ind(~any(neigh_ind,2),:) = [];
    neigh_vec(~any(neigh_vec,2),:) = [];
    
    a=1;angle_index=zeros(1,12);
    for v=1:size(neigh_ind,1)
        for w=1:size(neigh_ind,1) % From v or from 1?
            angle=rad2deg(atan2(norm(cross(neigh_vec(v,:),neigh_vec(w,:))),dot(neigh_vec(v,:),neigh_vec(w,:))));
            if nargin > 7
                angle_lo=angle_limit;
                angle_hi=angle;
            else
                angle_hi=angle_limit;
                angle_lo=angle;
            end
            if angle > 0 && angle_lo < angle_hi
                if neigh_ind(v,1) < neigh_ind(w,1)
                    angle_index(a,1)= atom2(neigh_ind(v,1)).index;
                    angle_index(a,2)= solute_index;
                    angle_index(a,3)= atom2(neigh_ind(w,1)).index;
                    angle_index(a,4)= angle;
                    angle_index(a,5)= neigh_dist(v);
                    angle_index(a,6)= neigh_dist(w);
                    angle_index(a,7:9)= neigh_vec(v,:);
                    angle_index(a,10:12)= neigh_vec(w,:);
                    a=a+1;
                else
                    angle_index(a,1)= atom2(neigh_ind(w,1)).index;
                    angle_index(a,2)= solute_index;
                    angle_index(a,3)= atom2(neigh_ind(v,1)).index;
                    angle_index(a,4)= angle;
                    angle_index(a,5)= neigh_dist(w);
                    angle_index(a,6)= neigh_dist(v);
                    angle_index(a,7:9)= neigh_vec(w,:);
                    angle_index(a,10:12)= neigh_vec(v,:);
                    a=a+1;
                end
            end
        end
    end
    
    [Y,I]=sort(angle_index(:,2));
    angle_index=angle_index(I,:);
    angle_index = unique(angle_index,'rows','stable');
    angle_index(~any(angle_index,2),:) = [];
    
    neigh.angle=angle_index;
    
    atom1(i).neigh.type=[atom2(neigh.in).type]';
    atom1(i).neigh.index=[atom2(neigh.in).index]';
    atom1(i).neigh.dist=neigh.dist;
    atom1(i).neigh.coords=neigh.coords;
    atom1(i).neigh.vec=neigh.r_vec;
    atom1(i).neigh.angle=neigh.angle;
    
end

% assignin('caller','Neigh',neigh)
% assignin('caller','Angle_index',angle_index)

end

