%% adjust_H_atom.m
% * This function tries to set the angles of structural H2O and bond
% * distances of X-H to reasonable values.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=adjust_H_atom(atom,Box_dim) % Basic input
% # atom=adjust_H_atom(atom,Box_dim,ideal_H_distance) % Sets the X-H bond distance
% # atom=adjust_H_atom(atom,Box_dim,ideal_H_distance,Bond_index) % So that we do not need to recalculate the Bond_index


function atom=adjust_H_atom(atom,Box_dim,varargin)
%%

if nargin < 3
    ideal_dist=0.957;

else
    ideal_dist=varargin{1};
end

ideal_angle_dist=1.514;
if ideal_dist==1
    ideal_angle_dist=1.581;
end

if nargin > 3
    niter=varargin{2};
else
    niter=2;
end

Bond_index=[];
Angle_index=[];
% To get the Bond_index variable
if nargin < 4
    atom_temp=bond_atom(atom,Box_dim,1.25);
else
    Bond_index=varargin{4};
end

% Bond_index(:,3)=[];

% To find all H's
H_ind=find(strncmpi([atom.type],'H',1));

for j=1:niter
    % Correct angles for H-X-H
    i=1;
    while i<size(atom,2)+1
        if strncmpi(atom(i).type,'H',1)
            if ~strncmpi(atom(i).type,'Hw',2)
                molid_ind=find([atom.molid]==[atom(i).molid]);
                %                 Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,1.7);
                Neigh = neighbor_atom([[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,1.25);
                H_Neigh_index=Neigh.in([ismember(Neigh.in,H_ind)]);
                ind_H=find(ismember(Neigh.in,H_ind));
                if [atom(H_Neigh_index).molid]==[atom(i).molid]
                    NewIdealCoords=num2cell([atom(H_Neigh_index).x atom(H_Neigh_index).y atom(H_Neigh_index).z] - 1/2*(ideal_angle_dist/Neigh.dist(ind_H)-1)*[-Neigh.r_vec(ind_H,1) - Neigh.r_vec(ind_H,2) - Neigh.r_vec(ind_H,3)]/norm([-Neigh.r_vec(ind_H,1) - Neigh.r_vec(ind_H,2) -Neigh.r_vec(ind_H,3)]));
                    [atom(H_Neigh_index).x atom(H_Neigh_index).y atom(H_Neigh_index).z]=deal(NewIdealCoords{:});
                end
            end
        end
        i=i+1;
        if mod(i,100)==1
            i-1
        end
    end

    % Correct bonds for X-H
    i=1;All_X_ind=[];
    while i<size(atom,2)+1
        if strncmpi(atom(i).type,'H',1)
            if ~strncmpi(atom(i).type,'Hw',2)
                [Bond_Index_row,Bond_Index_col]=find(Bond_index==i);
                if Bond_Index_col==1
                    X_ind=Bond_index(Bond_Index_row,2);
                else
                    X_ind=Bond_index(Bond_Index_row,1);
                end
                %                     All_X_ind=[All_X_ind X_ind];
                if length(X_ind)>0
                    Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom(X_ind).x]' [atom(X_ind).y]' [atom(X_ind).z]'],Box_dim,2.0);
                    [dist,dist_ind]=min([Neigh.dist]);
                    if numel(dist_ind)>2
                        disp('H bonded to more than one atom')
                        pause(1)
                    end
                    NewIdealCoords=num2cell([atom(i).x atom(i).y (atom(i).z)] - (ideal_dist/dist-1)*[Neigh.r_vec(1,1) Neigh.r_vec(1,2) Neigh.r_vec(1,3)]);
                    [atom(i).x atom(i).y atom(i).z]=deal(NewIdealCoords{:});
                else
                    disp('Found lonely H')
                end
            end
        end
        i=i+1;
        if mod(i,100)==1
            i-1
        end
    end
end

end



