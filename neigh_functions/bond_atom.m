%% bond_atom.m
% * This function tries to assign all bonds to a Bond_matrix and a
% Bond_index variable
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=bond_atom(atom,Box_dim) % Basic input arguments
% # atom=bond_atom(atom,Box_dim,2.25) % Allows setting the max cutoff
% # atom=bond_atom(atom,Box_dim,Bond_index) % Uses an exisitng Bond_index
%
function atom = bond_atom(atom,Box_dim,varargin)

rmaxshort=1.18;

if nargin>2
    rmaxlong=varargin{1}; % Dummy value
else
    rmaxlong=2.25;
end

if nargin>3
    distance_factor=varargin{2};
else
    distance_factor=0.65; % 1.3
end
% distance_factor=1.25; % due to ionic radii and not vdw..

XYZ_labels=[atom.type]';
XYZ_data=single([[atom.x]' [atom.y]' [atom.z]']); % use of single instead of double

[atom.fftype]=atom.type;

if ~isfield(atom,'element')
    atom = element_atom(atom);
end

[atom.type]=atom.element;

Radiiproperties=load('general_functions/Revised_Shannon_radii.mat');
% atom=bond_valence_atom(atom,Box_dim,rmaxshort,2.25);
if size(atom,2)<20000
    disp('Calculating the full distance matrix')
    dist_matrix = dist_matrix_atom(atom,Box_dim); % To calculate a full distance matrix
else
    disp('Calculating the distance matrix with cell lists')
    [dist_matrix,bond_list, dist_list,X_dist,Y_dist,Z_dist] = cell_list_dist_matrix_atom(atom, Box_dim,rmaxshort,rmaxlong);
end

XYZ_radii=zeros(length(XYZ_labels),1);
XYZ_formalcharge=zeros(length(XYZ_labels),1);
Atom_label=sort(unique([atom.type]));
for i=1:length(Atom_label)
    try
        ind=find(strncmpi([Radiiproperties.Ion],Atom_label(i),2));
    catch
        ind=find(strncmpi([Radiiproperties.Ion],Atom_label(i),1));
    end
    % XYZ_radii(ismember([atom.type],Atom_label(i)))=median(Radiiproperties.CrysRadii(ind))';
    % Atom_label(i)
    temp_radii=radius_vdw(Atom_label(i));
    %temp_radii=radius_ion(Atom_label(i));
    XYZ_radii(ismember([atom.type],Atom_label(i)))=temp_radii(1);
    XYZ_formalcharge(ismember([atom.type],Atom_label(i)))=median(Radiiproperties.OxState(ind))';
end

assignin('caller','XYZ_radii',XYZ_radii);
assignin('caller','XYZ_formalcharge',XYZ_formalcharge);

XYZ_radii(XYZ_radii==0)=distance_factor;
XYZ_radii(ismember([atom.type],'H'))=0.3; % Special H radii
radius_matrix=repmat(XYZ_radii,1,length(XYZ_radii));
radius_limit=(radius_matrix+radius_matrix')*distance_factor;
indH=strncmp([atom.type],'H',1);
% radius_limit(indH,:)=rmaxshort;radius_limit(:,indH)=rmaxshort;
radius_limit(radius_limit>rmaxlong)=rmaxlong;
dist_matrix(dist_matrix>radius_limit)=0;

if isfield(atom,'neigh')
    atom=rmfield(atom,'neigh');
end
if isfield(atom,'bond')
    atom=rmfield(atom,'bond');
end

% disp('Looking for neighbours/bonds')
Bond_index=single(zeros(1,3));
Angle_index=single(zeros(1,3));
a=1;b=1;i=1;
while i<size(atom,2)+1
    k=0;j=1;
    Neigh_ind=zeros(12,1);Neigh_vec=zeros(12,3); %%
    bond_ind=find(dist_matrix(:,i)>0);

    [atom(i).neigh.dist]=[];
    [atom(i).neigh.index]=[];
    [atom(i).neigh.type]={};
    [atom(i).neigh.coords]=[];
    [atom(i).neigh.r_vec]=[];
    [atom(i).bond.dist]=[];
    [atom(i).bond.index]=[];
    [atom(i).bond.type]={};

    while j <= numel(bond_ind) && k <= numel(bond_ind) %<= neigh %atom(i).neigh=[];
        if dist_matrix(bond_ind(j),i)>0
            %             if XYZ_formalcharge(i)*XYZ_formalcharge(bond_ind(j))<=0 && atom(i).molid==atom(bond_ind(j)).molid
            if atom(i).molid==atom(bond_ind(j)).molid
                k=k+1;
                [atom(i).neigh.dist(k,1)]=dist_matrix(bond_ind(j),i);
                [atom(i).neigh.index(k,1)]=bond_ind(j);
                [atom(i).neigh.type(k,1)]=XYZ_labels(bond_ind(j));
                [atom(i).neigh.coords(k,:)]=[XYZ_data(bond_ind(j),1) XYZ_data(bond_ind(j),2) XYZ_data(bond_ind(j),3)];
                [atom(i).neigh.r_vec(k,:)]=[X_dist(bond_ind(j),i) Y_dist(bond_ind(j),i) Z_dist(bond_ind(j),i)];
                if [atom(i).molid]==[atom(bond_ind(j)).molid] && dist_matrix(bond_ind(j),i)<rmaxlong
                    [atom(i).bond.dist(k,1)]=dist_matrix(bond_ind(j),i);
                    [atom(i).bond.index(k,:)]=[i bond_ind(j)];
                    [atom(i).bond.type]=1;
                    Bond_index(b,1)=min([i bond_ind(j)]);
                    Bond_index(b,2)=max([i bond_ind(j)]);
                    Bond_index(b,3)=[atom(i).bond.dist(k,1)];
                    %                     b=b+1;
                    %                 end
                    Neigh_ind(b,1) = bond_ind(j);%%
                    Neigh_vec(b,1:3) = -[atom(i).neigh.r_vec(k,:)]; %%
                    b=b+1; %%

                end

            end
        end
        j=j+1;
    end

    %%

    Neigh_ind(~any(Neigh_ind,2),:) = [];
    Neigh_vec(~any(Neigh_vec,2),:) = [];
    for v=1:size(Neigh_ind,1)
        for w=1:size(Neigh_ind,1) % From v or from 1?
            angle=rad2deg(atan2(norm(cross(Neigh_vec(v,:),Neigh_vec(w,:))),dot(Neigh_vec(v,:),Neigh_vec(w,:))));
            if angle > 0 && angle <= 180 % Do we need this??
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

    % if ~isempty(Angle_index)
    %     if ismember(i,Angle_index(:,1:3))
    %         %                 [C,D]=find(Angle_index(:,1:3)==i);
    %         [C,D]=find(Angle_index(:,2)==i);
    %         atom(i).angle.type = 1;
    %         atom(i).angle.index = Angle_index(C,1:3);
    %         atom(i).angle.angle = Angle_index(C,4);
    %         atom(i).angle.vec1 = Angle_index(C,5:7);
    %         atom(i).angle.vec2 = Angle_index(C,8:10);
    %     end
    % end


    %%

    % if mod(i,1000)==1
    %     if i > 1
    %         i-1
    %     end
    % end
    i=i+1;
end
% i-1

[Y,I]=sort(Bond_index(:,1));
Bond_index=Bond_index(I,:);
Bond_index = unique(Bond_index,'rows','stable');
Bond_index(~any(Bond_index,2),:) = [];
nBonds=size(Bond_index,1);

nAngles=0;
if ~isempty(Angle_index)
    [Y,I]=sort(Angle_index(:,2));
    Angle_index=Angle_index(I,:);
    Angle_index = unique(Angle_index,'rows','stable');
    Angle_index(~any(Angle_index,2),:) = [];
    nAngles=size(Angle_index,1);
    assignin('caller','Angle_index',Angle_index);
end

i=1;
while i<size(atom,2)+1
    if ~isempty(Angle_index)
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
    i=i+1;
end

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

rm_ind=[rm_ind];
for i=1:size(Neigh_index,1)
    if [atom(Neigh_index(i,1)).molid]~=[atom(Neigh_index(i,2)).molid]
        rm_ind=[rm_ind; i];
    end
end
Neigh_index(rm_ind,:)=[];

assignin('caller','Neigh_ind',Neigh_ind);
assignin('caller','Neigh_vec',Neigh_vec);

[atom.type]=atom.fftype;
atom=order_attributes(atom);

assignin('caller','nBonds',nBonds);
assignin('caller','nAngles',nAngles);
assignin('caller','radius_matrix',radius_matrix);
assignin('caller','radius_limit',radius_limit);
assignin('caller','Bond_index',Bond_index);
assignin('caller','Neigh_index',Neigh_index);
% assignin('caller','bond_matrix',dist_matrix);
assignin('caller','dist_matrix',dist_matrix);
assignin('caller','X_dist',X_dist);
assignin('caller','Y_dist',Y_dist);
assignin('caller','Z_dist',Z_dist);

