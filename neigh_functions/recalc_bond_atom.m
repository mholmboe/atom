%% recalc_bond_atom.m
% * This function recalculates the bonds, angles and dihedrals (latter two
% * optional), based on already existing Bond | Angle | Dihedral_index
% * variables.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=recalc_bond_atom(atom,Box_dim,Bond_index) % Basic input arguments
% # atom=recalc_bond_atom(atom,Box_dim,Bond_index,Angle_index) % Also recalculates the angles
% # atom=recalc_bond_atom(atom,Box_dim,Bond_index,Angle_index,Dihedral_index) % Also recalculates the dihedrals
% # atom=recalc_bond_atom(atom,Box_dim,Bond_index,Angle_index,Dihedral_index,rmaxlong) % Allows setting the max cutoff

function atom = recalc_bond_atom(atom,Box_dim,varargin)

rmaxlong=2.25; % May be reset below

XYZ_labels=[atom.type]';
XYZ_data=single([[atom.x]' [atom.y]' [atom.z]']);

if nargin==2
    disp('You need to atleast supply the Bond_index variable')
elseif nargin==3
    Bond_index=varargin{1};
    Angle_index=[];
    Dihedral_index=[];
elseif nargin==4
    Bond_index=varargin{1};
    Angle_index=varargin{2};
    Dihedral_index=[];
elseif nargin==5
    Bond_index=varargin{1};
    Angle_index=varargin{2};
    Dihedral_index=varargin{3};
elseif nargin>5
    Bond_index=varargin{1};
    Angle_index=varargin{2};
    Dihedral_index=varargin{3};
    rmaxlong=varargin{4};
end

disp('Calculating the distance matrix')


if size(atom,2)>50000 && numel(Box_dim)<9
    dist_matrix = cell_list_dist_matrix_atom(atom,Box_dim);%,1.25,4);
elseif size(atom,2) < 50000
    dist_matrix = dist_matrix_atom(atom,Box_dim);%,1.25,4);
else
    disp('The nAtoms in the atoms truct may be too large, try reducing the number of atoms?')
end

assignin('base','dist_matrix',dist_matrix);

if nargin>2
    disp('Re-calculating the Bond_index')
    Bond_index=Bond_index(:,1:3);
    for i=1:size(Bond_index,1)
        Bond_index(i,3)=dist_matrix(Bond_index(i,1),Bond_index(i,2));
    end
    nBonds=size(Bond_index,1);
    assignin('caller','Bond_index',Bond_index);
    assignin('caller','nBonds',nBonds);

    if nargin>3
        if size(Angle_index,2)>1
            disp('Re-calculating the Angle_index')
            Angle_index=Angle_index(:,1:4);
            for i=1:size(Angle_index,1)
                vec1=[X_dist(Angle_index(i,2),Angle_index(i,1)) Y_dist(Angle_index(i,2),Angle_index(i,1)) Z_dist(Angle_index(i,2),Angle_index(i,1))];
                vec2=[X_dist(Angle_index(i,2),Angle_index(i,3)) Y_dist(Angle_index(i,2),Angle_index(i,3)) Z_dist(Angle_index(i,2),Angle_index(i,3))];
                Angle_index(i,4)=rad2deg(atan2(norm(cross(vec1,vec2)),dot(vec1,vec2)));
                Angle_index(i,5:7)=vec1;
                Angle_index(i,8:10)=vec2;
            end
            nAngles=size(Angle_index,1);
            assignin('caller','nAngles',nAngles);
            assignin('caller','Angle_index',Angle_index);
        end
    else
        Angle_index=[0 0 0 0];
    end

    if nargin>4
        if size(Dihedral_index,2)>1
            disp('Re-calculating the Dihedral_index')
            Ax2=[[Angle_index(:,3) Angle_index(:,2) Angle_index(:,1) Angle_index(:,4) Angle_index(:,8:10) Angle_index(:,5:7)]; Angle_index];
            d=1;
            for i=1:size(Ax2,1)
                for j=i:size(Ax2,1)
                    if isequal([Ax2(i,2) Ax2(i,3)],[Ax2(j,1) Ax2(j,2)])
                        A=cross([Ax2(i,5) Ax2(i,6) Ax2(i,7)],[Ax2(i,8) Ax2(i,9) Ax2(i,10)]);
                        B=cross([Ax2(j,5) Ax2(j,6) Ax2(j,7)],[Ax2(j,8) Ax2(j,9) Ax2(j,10)]);
                        normA=sqrt(sum(A.*A,2));
                        normB=sqrt(sum(B.*B,2));
                        theta=rad2deg(acos(dot(A,B)./(normA*normB)));
                        if Ax2(i,2)<Ax2(i,3)
                            Dihedral_index(d,1:5)=[Ax2(i,1) Ax2(i,2) Ax2(i,3) Ax2(j,3) round2dec(theta,2)];
                        else
                            Dihedral_index(d,1:5)=[Ax2(j,3) Ax2(i,3) Ax2(i,2) Ax2(i,1) round2dec(theta,2)];
                        end
                        d=d+1;
                    end
                end
                if mod(i,1000)==1
                    if i-1>0
                        i-1
                    end
                end
            end

        end

        nDihedrals=size(Dihedral_index,2);

        if nDihedrals>0
            [Y,I] = sort(Dihedral_index(:,2));
            Dihedral_index = Dihedral_index(I,:);
            Dihedral_index = unique(Dihedral_index,'rows','stable');
            Dihedral_index(~any(Dihedral_index,2),:) = [];
        else
            Dihedral_index =[];
        end
        nDihedrals=size(Dihedral_index,1);
        assignin('caller','nDihedrals',nDihedrals);
        assignin('caller','Dihedral_index',Dihedral_index);
    else
        Dihedral_index=[0 0 0 0 0];
    end
end
i

%%%%%%%%%%%%%%%%%%%%%%%

Bx2=[Bond_index;[Bond_index(:,2) Bond_index(:,1) Bond_index(:,3)]];
Bx2=sortrows(Bx2);
disp('Looking for neighbours/bonds')
i=1;
while i<size(atom,2)+1
    k=1;j=1;
    ind=Bx2(:,1)==i;
    neigh_ind=Bx2(ind,2);

    [atom(i).neigh.dist] = [];
    [atom(i).neigh.index] = [];
    [atom(i).neigh.type] = {};
    [atom(i).neigh.coords] = [];
    [atom(i).neigh.r_vec] = [];
    [atom(i).bond.dist] = [];
    [atom(i).bond.index] = [];
    [atom(i).bond.type] = {};
    [atom(i).angle.type] = [];
    [atom(i).angle.index] = [];
    [atom(i).angle.angle] = [];
    [atom(i).angle.vec1] = [];
    [atom(i).angle.vec2] = [];

    while j <= numel(neigh_ind) && k <= numel(neigh_ind) %<= neigh %atom(i).neigh=[];
        if atom(i).molid==atom(neigh_ind(j)).molid
            [atom(i).neigh.dist(k,1)]=dist_matrix(neigh_ind(j),i);
            [atom(i).neigh.index(k,1)]=neigh_ind(j);
            [atom(i).neigh.type(k,1)]=XYZ_labels(neigh_ind(j));
            [atom(i).neigh.coords(k,:)]=[XYZ_data(neigh_ind(j),1) XYZ_data(neigh_ind(j),2) XYZ_data(neigh_ind(j),3)];
            [atom(i).neigh.r_vec(k,:)]=[X_dist(neigh_ind(j),i) Y_dist(neigh_ind(j),i) Z_dist(neigh_ind(j),i)];
            k=k+1;
            j=j+1;
        end
    end
    if ismember(i,Bond_index(:,1:2))
        [A,B]=find(Bond_index(:,1:2)==i);
        atom(i).bond.type = 1;
        atom(i).bond.index = Bond_index(A,1:2);
        atom(i).bond.dist = Bond_index(A,3);
        if size(Angle_index,2)>0 && ismember(i,Angle_index(:,1:3))
            [C,D]=find(Angle_index(:,2)==i);
            atom(i).angle.type = 1;
            atom(i).angle.index = Angle_index(C,1:3);
            atom(i).angle.angle = Angle_index(C,4);
            atom(i).angle.vec1 = Angle_index(C,5:7);
            atom(i).angle.vec2 = Angle_index(C,8:10);
        end
    end
    if mod(i,1000)==1
        if i-1>0
            i-1
        end
    end
    i=i+1;
end
i-1

%% New addition v3.00

CoordNumber=zeros(1,size(atom,2));Remove_ind=0;
if length(Bond_index)>0
    [Y,i] = sort(Bond_index(:,1));
    Bond_index = Bond_index(i,:);
    Bond_index = unique(Bond_index,'rows','stable');
    try
        CoordNumber=arrayfun(@(x) numel(x.neigh.index),atom);
    catch
        CoordNumber=zeros(1,size(atom,2));
        for i=1:size(atom,2)
            i
            [atom(i).neigh.index]
            CoordNumber(i)=numel([atom(i).neigh.index]);
        end
    end
    Remove_ind=find(CoordNumber==0);
    % assignin('caller','Remove_ind',Remove_ind);
    % assignin('caller','CoordNumber',CoordNumber);
end
assignin('caller','Remove_ind',Remove_ind);
assignin('caller','CoordNumber',CoordNumber);

ind=find(tril(dist_matrix)>0);
r=dist_matrix(ind);
[i,j] = ind2sub(size(dist_matrix),ind);

Neigh_index = [j i r];
[Y,i] = sort(Neigh_index(:,1));
Neigh_index = Neigh_index(i,:);
Neigh_index = unique(Neigh_index,'rows','stable');
rm_ind=find(Neigh_index(:,3)>rmaxlong);

Neigh_index(unique([rm_ind]),:) = [];

% Old since v2.13, but too slow in Octave
% rm_ind=[rm_ind];
% for i=1:size(Neigh_index,1)
%     if [atom(Neigh_index(i,1)).molid]~=[atom(Neigh_index(i,2)).molid]
%         rm_ind=[rm_ind; i];
%     end
% end
% Neigh_index(rm_ind,:)=[];

% Extract molid for the first and second columns of Neigh_index
molid1 = [atom(Neigh_index(:,1)).molid];
molid2 = [atom(Neigh_index(:,2)).molid];

% Find rows where molid values differ
rm_ind2 = molid1 ~= molid2;

% Remove rows in one step
Neigh_index(unique([rm_ind2]),:) = [];

[Y,I]=sort(Bond_index(:,1));
Bond_index=Bond_index(I,:);
Bond_index = unique(Bond_index,'rows','stable');
Bond_index(~any(Bond_index,2),:) = [];
nBonds=size(Bond_index,1);

% [Y,I]=sort(Angle_index(:,2));
% Angle_index=Angle_index(I,:);
% Angle_index = unique(Angle_index,'rows','stable');
% Angle_index(~any(Angle_index,2),:) = [];
% nAngles=size(Angle_index,1);
% assignin('caller','Angle_index',Angle_index);
% assignin('caller','Neigh_ind',Neigh_ind);
% assignin('caller','Neigh_vec',Neigh_vec);

[atom.type]=atom.fftype;
atom=order_attributes(atom);

assignin('caller','nBonds',nBonds);
% assignin('caller','radius_limit',radius_limit);
assignin('caller','Bond_index',Bond_index);
assignin('caller','Neigh_index',Neigh_index);
% assignin('caller','bond_matrix',dist_matrix);
assignin('caller','dist_matrix',dist_matrix);
assignin('caller','X_dist',X_dist);
assignin('caller','Y_dist',Y_dist);
assignin('caller','Z_dist',Z_dist);

%% End New addition v3.00

if nargin>6
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

            assignin('base',strcat(char(Atom_labels(i)),'_dist')',[num2cell(Tot_index) num2cell(Tot_neighindex) Tot_type num2cell(Tot_dist)]);
            assignin('base',strcat(char(Atom_labels(i)),'_coords')',[[atom(Tot_neighindex).x]' [atom(Tot_neighindex).y]' [atom(Tot_neighindex).z]']);
            assignin('base',strcat(char(Atom_labels(i)),'_bonds')',[Tot_bondindex Tot_bonds]);
            assignin('base',strcat(char(Atom_labels(i)),'_angles')',[Tot_angleindex Tot_angles]);
            assignin('base',strcat(char(Atom_labels(i)),'_atom')',atom(ismember([atom.type],Atom_labels(i))));
        catch
            assignin('base',strcat(char(Atom_labels(i)),'_dist')',[num2cell(Tot_index) num2cell(Tot_neighindex) Tot_type num2cell(Tot_dist)]);
            assignin('base',strcat(char(Atom_labels(i)),'_coords')',[[atom(Tot_neighindex).x]' [atom(Tot_neighindex).y]' [atom(Tot_neighindex).z]']);
            assignin('base',strcat(char(Atom_labels(i)),'_bonds')',[Tot_bondindex Tot_bonds]);
            assignin('base',strcat(char(Atom_labels(i)),'_angles')',[Tot_angleindex Tot_angles]);
            assignin('base',strcat(char(Atom_labels(i)),'_atom')',atom(ismember([atom.type],Atom_labels(i))));
        end
    end
end

end

