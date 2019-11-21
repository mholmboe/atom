%% bond_angle_atom.m
% * This function tries to find all bonds and angles of the atom struct
% * One optional argument like 'more' will give more bond info
% * atom is the atom struct
% * Box_dim is the box dimension vector
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=bond_angle_atom(atom,Box_dim,1.25,2.25)
% # atom=bond_angle_atom(atom,Box_dim,1.25,2.25,'more')
%

function atom=bond_angle_atom(atom,Box_dim,rmaxshort,rmaxlong,varargin)
%%
% tic

nAtoms=size(atom,2);
% max_short_dist=1.2;
% max_long_dist=2.3;

XYZ_labels=[atom.type]';
XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];

lx=Box_dim(1);ly=Box_dim(2);lz=Box_dim(3);
if length(Box_dim) > 3
    xy=Box_dim(6);xz=Box_dim(8);yz=Box_dim(9);
else
    xy=0;xz=0;yz=0;
end

close_count=1;
Bond_index=zeros(4*size(XYZ_data,1),3);
Angle_index=zeros(4*size(XYZ_data,1),4);
dist_matrix=zeros(nAtoms,nAtoms);
% dist_matrix = dist_matrix_atom(atom,Box_dim); % Do this if you want to calc dist_matrix_atom first
b=1;a=1; overlap_index=[];
for i = 1:size(XYZ_data,1)
    
    %     disp('Setting max_neigh_distance to... ')
    %     max_neigh_distance = max_long_dist %2.25; % Set large max_neigh_distance so that we can find Si-O-H and Al-OH2 bonds %max_long_dist;%2.3
    %     disp('If you have problems, change back to 2.25 on line 43 in bond_angle_atom()')
    %     max_neigh_distance = max_long_dist; % Change to 2.25 if adding edge bonds
    
    if strcmpi(strtrim(XYZ_labels(i)),'Oalhh') || strcmp(strtrim(XYZ_labels(i)),'Osih')
        if rmaxlong < 2.25
            max_neigh_distance = 2.25;
        else
            max_neigh_distance = rmaxlong;
        end
    else
        max_neigh_distance = rmaxlong;
    end
    %                 disp('Looking for M-O-H bonds')
    
    if length(Box_dim)>3
        %Calculate Distance Components for triclic cell
        rz = XYZ_data(i,3) - XYZ_data(:,3);
        ry = XYZ_data(i,2) - XYZ_data(:,2);
        rx = XYZ_data(i,1) - XYZ_data(:,1);
        
        z_gt_ind=find(rz > lz/2); z_lt_ind=find(rz < - lz/2);
        rz(z_gt_ind) = rz(z_gt_ind) - lz;
        rz(z_lt_ind) = rz(z_lt_ind) + lz;
        rx(z_gt_ind) = rx(z_gt_ind) - xz;
        rx(z_lt_ind) = rx(z_lt_ind) + xz;
        ry(z_gt_ind) = ry(z_gt_ind) - yz;
        ry(z_lt_ind) = ry(z_lt_ind) + yz;
        
        y_gt_ind=find(ry > ly/2); y_lt_ind=find(ry < - ly/2);
        ry(y_gt_ind) = ry(y_gt_ind) - ly;
        ry(y_lt_ind) = ry(y_lt_ind) + ly;
        rx(y_gt_ind) = rx(y_gt_ind) - xy;
        rx(y_lt_ind) = rx(y_lt_ind) + xy;
        
        x_gt_ind=find(rx > lx/2); x_lt_ind=find(rx < - lx/2);
        rx(x_gt_ind) = rx(x_gt_ind) - lx;
        rx(x_lt_ind) = rx(x_lt_ind) + lx;
    else
        %Calculate Distance Components for ortogonal cell
        
        rz = XYZ_data(i,3) - XYZ_data(:,3);
        z_gt_ind=find(rz > lz/2); z_lt_ind=find(rz < - lz/2);
        rz(z_gt_ind) = rz(z_gt_ind) - lz;
        rz(z_lt_ind) = rz(z_lt_ind) + lz;
        
        ry = XYZ_data(i,2) - XYZ_data(:,2);
        y_gt_ind=find(ry > ly/2); y_lt_ind=find(ry < - ly/2);
        ry(y_gt_ind) = ry(y_gt_ind) - ly;
        ry(y_lt_ind) = ry(y_lt_ind) + ly;
        
        rx = XYZ_data(i,1) - XYZ_data(:,1);
        x_gt_ind=find(rx > lx/2); x_lt_ind=find(rx < - lx/2);
        rx(x_gt_ind) = rx(x_gt_ind) - lx;
        rx(x_lt_ind) = rx(x_lt_ind) + lx;
        
    end
    
    r = sqrt( rx(:,1).^2 + ry(:,1).^2 + rz(:,1).^2 ); % distance calc.
    bond_in=intersect(find(r > 0), find(r < max_neigh_distance));
    dist_matrix(:,i)=r;
    
    
    % Else if we already called dist_matrix_atom...
    %     r = dist_matrix(:,i);
    %     rx= X_dist(:,i);
    %     ry= Y_dist(:,i);
    %     rz= Z_dist(:,i);
    %     bond_in=intersect(find(r > 0), find(r < max_neigh_distance));
    %
    
    
    n=1;
    Neigh_ind=zeros(12,1);Neigh_vec=zeros(12,3);
    for j=1:length(bond_in)
        if [atom(i).molid]==[atom(bond_in(j)).molid]
            % Original, uncomment this section
            % max_distance = max_short_dist;
            % Only test
            if strncmpi(strtrim(XYZ_labels(i)),'H',1) || strncmpi(strtrim(XYZ_labels(bond_in(j))),'H',1)
                max_distance = rmaxshort;
            elseif strcmpi(strtrim(XYZ_labels(i)),'Oalhh') || strcmp(strtrim(XYZ_labels(i)),'Osih')
                if rmaxlong < 2.25
                    max_distance = 2.25;
                else
                    max_distance = rmaxlong;
                end
                %                 disp('Looking for M-O-H bonds')
            else
                max_distance = rmaxlong;
            end
            
            if r(bond_in(j)) < 0.6
                disp('Atoms too close!!!')
                r(bond_in(j))
                [i,bond_in(j)]
                XYZ_labels(i)
                XYZ_labels(bond_in(j))
                XYZ_data(i,:)
                XYZ_data(bond_in(j),:)
                overlap_index=[overlap_index; {i bond_in(j) r(bond_in(j)) XYZ_labels(i) XYZ_labels(bond_in(j)) XYZ_data(i,:) XYZ_data(bond_in(j),:)}];
            elseif r(bond_in(j)) > max_distance && r(bond_in(j)) < 1.25
                disp('Atoms pretty close...')
                r(bond_in(j))
                [i,bond_in(j)]
                XYZ_labels(i)
                XYZ_labels(bond_in(j))
                XYZ_data(i,:)
                XYZ_data(bond_in(j),:)
                close_count=close_count+1;
                overlap_index=[overlap_index; {i bond_in(j) r(bond_in(j)) XYZ_labels(i) XYZ_labels(bond_in(j)) XYZ_data(i,:) XYZ_data(bond_in(j),:)}];
            end
            if r(bond_in(j)) > max_distance/3 && r(bond_in(j)) < max_distance %strncmpi(XYZ_labels(i),strtrim(XYZ_labels(j)),1) == 0;
                if strncmpi(strtrim(XYZ_labels(i)),{'OW'},2)% || strncmpi(strtrim(XYZ_labels(bond_in(j))),{'OW'},2); % < bond_in(j) && ismember(XYZ_labels(i),{'Ow','Hw','OW','HW','HW1','HW2'}) > 0;
                    if bond_in(j) > i && bond_in(j) < i+3 %strncmpi(strtrim(XYZ_labels(i)),{'Ow'},2) && bond_in(j) > i && bond_in(j) < i+3;
                        Bond_index(b,1)= min([i bond_in(j)]);
                        Bond_index(b,2)= max([i bond_in(j)]);
                        Bond_index(b,3)= r(bond_in(j));
                        b=b+1;
                        if r(bond_in(j)) < rmaxshort % This should always be true right?
                            Neigh_ind(n,1)= bond_in(j);
                            Neigh_vec(n,1:3) = [rx(bond_in(j)) ry(bond_in(j)) rz(bond_in(j))];
                            n=n+1;
                        end
                    end
                elseif strncmpi(strtrim(XYZ_labels(i)),{'HW'},2)
                    % disp('HW')
                elseif strncmpi(strtrim(XYZ_labels(i)),{'H'},1) && strncmpi(strtrim(XYZ_labels(bond_in(j))),{'H'},1)
                    %                     disp('disallowed a H-H bond')
                elseif ismember({'H'},[XYZ_labels(i) XYZ_labels(bond_in(j))]) > 0 || ismember({'Oalhh'},[XYZ_labels(i) XYZ_labels(bond_in(j))]) > 0 || ismember({'Osih'},[XYZ_labels(i) XYZ_labels(j)]) > 0 && ismember(XYZ_labels(bond_in(j)),{'Ow','Hw','OW','HW','HW1','HW2'}) == 0
                    Bond_index(b,1)= min([i bond_in(j)]);
                    Bond_index(b,2)= max([i bond_in(j)]);
                    Bond_index(b,3)= r(bond_in(j));
                    b=b+1;
                    if r(bond_in(j)) < max_distance % This should always be true right?
                        Neigh_ind(n,1)= bond_in(j);
                        Neigh_vec(n,1:3) = [rx(bond_in(j)) ry(bond_in(j)) rz(bond_in(j))];
                        n=n+1;
                    end
                elseif ~strncmpi(strtrim(XYZ_labels(i)),{'OW'},2) || ~strncmpi(strtrim(XYZ_labels(i)),{'HW'},2)
                    if bond_in(j) < i
                        Bond_index(b,1)= min([i bond_in(j)]);
                        Bond_index(b,2)= max([i bond_in(j)]);
                        Bond_index(b,3)= r(bond_in(j));
                        b=b+1;
                    end
                    if r(bond_in(j)) < max_distance % This should always be true right?
                        Neigh_ind(n,1)= bond_in(j);
                        Neigh_vec(n,1:3) = [rx(bond_in(j)) ry(bond_in(j)) rz(bond_in(j))];
                        n=n+1;
                    end
                end
            end
            
        end % test
        
        Neigh_ind(~any(Neigh_ind,2),:) = [];
        Neigh_vec(~any(Neigh_vec,2),:) = [];
        
        for v=1:size(Neigh_ind,1)
            for w=1:size(Neigh_ind,1) % From v or from 1?
                angle=rad2deg(atan2(norm(cross(Neigh_vec(v,:),Neigh_vec(w,:))),dot(Neigh_vec(v,:),Neigh_vec(w,:))));
                if angle > 0 && angle <= 150 % Do we need this??
                    if v < w
                        Angle_index(a,1)= Neigh_ind(v,1);
                        Angle_index(a,2)= i;
                        Angle_index(a,3)= Neigh_ind(w,1);
                        Angle_index(a,4)= angle;
                        Angle_index(a,5:7)= Neigh_vec(v,:);
                        Angle_index(a,8:10)= Neigh_vec(w,:);
                        a=a+1;
                    else
                        Angle_index(a,1)= Neigh_ind(w,1);
                        Angle_index(a,2)= i;
                        Angle_index(a,3)= Neigh_ind(v,1);
                        Angle_index(a,4)= angle;
                        Angle_index(a,5:7)= Neigh_vec(w,:);
                        Angle_index(a,8:10)= Neigh_vec(v,:);
                        a=a+1;
                    end
                end
            end
        end
    end
    if mod(i,100)==1
        i-1
    end
end


[Y,I]=sort(Bond_index(:,1));
Bond_index=Bond_index(I,:);
Bond_index = unique(Bond_index,'rows','stable');

[Y,I]=sort(Angle_index(:,2));
Angle_index=Angle_index(I,:);
Angle_index = unique(Angle_index,'rows','stable');

overlap_index=overlap_index(1:size(overlap_index,1)/2,:);
% [Y,I]=sort(cell2mat(overlap_index(:,2)));
% overlap_index=overlap_index(I,:);
% overlap_index = unique(overlap_index,'rows','stable');

Bond_index(~any(Bond_index,2),:) = [];
Angle_index(~any(Angle_index,2),:) = [];
% overlap_index(~any(overlap_index,2),:) = [];


nBonds=size(Bond_index,1);
nAngles=size(Angle_index,1);

if nargin > 4 %% This will print a whole lot more info to the calling workspace
    % Set a cutoff vector
    rmax=repmat(max_neigh_distance,nAtoms,1);
    % Reduce the cutoff vector for H - X
    rmax(strncmpi(strtrim(XYZ_labels),'H',1))=1.25;
    for i=1:size(atom,2)
        Neigh_ind=intersect(find(dist_matrix(:,i)>0),find(dist_matrix(:,i)<rmax(i)));%radius_ion([atom.type])));
        Neigh_dist=dist_matrix(Neigh_ind,i);
        atom(i).neigh.type = [atom(Neigh_ind).type]';
        atom(i).neigh.index = Neigh_ind;
        atom(i).neigh.dist = Neigh_dist;
        atom(i).neigh.coords = [[atom(Neigh_ind).x]' [atom(Neigh_ind).y]' [atom(Neigh_ind).z]']; % PBC NOT TAKEN INTO ACCOUNT!!!
        atom(i).bond.type = [];
        atom(i).bond.index = [];
        atom(i).bond.dist = [];
        atom(i).angle.type = [];
        atom(i).angle.index = [];
        atom(i).angle.angle = [];
        atom(i).angle.vec1 = [];
        atom(i).angle.vec2 = [];
        
        if ismember(i,Bond_index(:,1:2))
            [A,B]=find(Bond_index(:,1:2)==i);
            atom(i).bond.type = 1;
            atom(i).bond.index = Bond_index(A,1:2);
            atom(i).bond.dist = Bond_index(A,3);
            if ismember(i,Angle_index(:,1:3))
                %                 [C,D]=find(Angle_index(:,1:3)==i);
                [C,D]=find(Angle_index(:,2)==i);
                atom(i).angle.type = 1;
                atom(i).angle.index = Angle_index(C,1:3);
                atom(i).angle.angle = Angle_index(C,4);
                atom(i).angle.vec1 = Angle_index(C,5:7);
                atom(i).angle.vec2 = Angle_index(C,8:10);
            end
        end
        if mod(i,100)==1
            i-1
        end
    end
    
    %%%%%%%%%%
    
    Atom_labels=unique([atom.type]);
    for i=1:length(Atom_labels)
        label_ind=find(strcmpi([atom.type],Atom_labels(i)));
        Tot_dist=[];Tot_type=[];Tot_index=[];Tot_angleindex=[];Tot_bondindex=[];Tot_neighindex=[];Tot_coords=[];Tot_bonds=[];Tot_angles=[];
        for j=label_ind
            if numel([atom(j).neigh])>0
                Tot_index=[Tot_index; repmat(j,numel([atom(j).neigh.index]),1)];
                Tot_dist=[Tot_dist; [atom(j).neigh.dist]];
                Tot_type=[Tot_type; [atom(j).neigh.type]];
                Tot_neighindex=[Tot_neighindex; [atom(j).neigh.index]];
                Tot_coords=[Tot_coords; [atom(j).neigh.coords]];
            end
            
            if numel([atom(j).bond])>0
                Tot_bondindex=[Tot_bondindex; [atom(j).bond.index]];
                Tot_bonds=[Tot_bonds; [atom(j).bond.dist]];
            end
            if numel([atom(j).angle])>0
                Tot_angleindex=[Tot_angleindex; [atom(j).angle.index]];
                Tot_angles=[Tot_angles; [atom(j).angle.angle]];
            end
        end
        try
            assignin('caller',strcat(char(Atom_labels(i)),'_dist')',[num2cell(Tot_index) num2cell(Tot_neighindex) Tot_type num2cell(Tot_dist)]);
            assignin('caller',strcat(char(Atom_labels(i)),'_coords')',[[atom(Tot_neighindex).x]' [atom(Tot_neighindex).y]' [atom(Tot_neighindex).z]']);
            assignin('caller',strcat(char(Atom_labels(i)),'_bonds')',[Tot_bondindex Tot_bonds]);
            assignin('caller',strcat(char(Atom_labels(i)),'_angles')',[Tot_angleindex Tot_angles]);
            assignin('caller',strcat(char(Atom_labels(i)),'_atom')',atom(ismember([atom.type],Atom_labels(i))));
        catch
            assignin('base',strcat(char(Atom_labels(i)),'_dist')',[num2cell(Tot_index) num2cell(Tot_neighindex) Tot_type num2cell(Tot_dist)]);
            assignin('base',strcat(char(Atom_labels(i)),'_coords')',[[atom(Tot_neighindex).x]' [atom(Tot_neighindex).y]' [atom(Tot_neighindex).z]']);
            assignin('base',strcat(char(Atom_labels(i)),'_bonds')',[Tot_bondindex Tot_bonds]);
            assignin('base',strcat(char(Atom_labels(i)),'_angles')',[Tot_angleindex Tot_angles]);
            assignin('base',strcat(char(Atom_labels(i)),'_atom')',atom(ismember([atom.type],Atom_labels(i))));
        end
    end
end

try
    assignin('caller','dist_matrix',dist_matrix);
    assignin('caller','overlap_index',overlap_index);
    assignin('caller','Bond_index',Bond_index);
    assignin('caller','Angle_index',Angle_index);
    assignin('caller','nBonds',nBonds);
    assignin('caller','nAngles',nAngles);
catch
    assignin('base','dist_matrix',dist_matrix);
    assignin('base','overlap_index',overlap_index);
    assignin('base','Bond_index',Bond_index);
    assignin('base','Angle_index',Angle_index);
    assignin('base','nBonds',nBonds);
    assignin('base','nAngles',nAngles);
end

% toc
