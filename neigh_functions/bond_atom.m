%% bond_atom.m
% * This function tries to assign all bonds to a Bond_matrix and a
% Bond_index variable
%
%% Version
% 2.09
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=bond_atom(atom,Box_dim)
% # atom=bond_atom(atom,Box_dim,2.25) % Last argument is the max bond distance

function atom = bond_atom(atom,Box_dim,varargin)

if nargin>2
    rmaxlong=varargin{1}; % Dummy value
else
    rmaxlong=2.25;
end

if nargin>3
    distance_factor=varargin{2};
else
    distance_factor=0.6; % 1.3
end

XYZ_labels=[atom.type]';
XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];

[atom.fftype]=atom.type;

if ~isfield(atom,'element')
    atom = element_atom(atom);
end

[atom.type]=atom.element;

Radiiproperties=load('Revised_Shannon_radii.mat');
% atom=bond_valence_atom(atom,Box_dim,1.25,2.25);
% clayff_param(sort(unique([atom.type])),'SPC/E');

disp('Calculating the distance matrix')
if size(atom,2)>20000
     dist_matrix = cell_list_dist_matrix_atom(atom,Box_dim,1.25,2.25);
else
     dist_matrix = dist_matrix_atom(atom,Box_dim,1.25,2.25);
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
    %     XYZ_radii(ismember([atom.type],Atom_label(i)))=median(Radiiproperties.CrysRadii(ind))';
    Atom_label(i)
    temp_radii=radius_vdw(Atom_label(i));
    XYZ_radii(ismember([atom.type],Atom_label(i)))=temp_radii(1);
    XYZ_formalcharge(ismember([atom.type],Atom_label(i)))=median(Radiiproperties.OxState(ind))';
end

assignin('caller','XYZ_radii',XYZ_radii);
assignin('caller','XYZ_formalcharge',XYZ_formalcharge);

XYZ_radii(XYZ_radii==0)=distance_factor;
radius_matrix=repmat(XYZ_radii,1,length(XYZ_radii));
radius_limit=(radius_matrix+radius_matrix')*distance_factor;
radius_limit(radius_limit>rmaxlong)=rmaxlong;
dist_matrix(dist_matrix>radius_limit)=0;

if isfield(atom,'neigh')
    atom=rmfield(atom,'neigh');
end
if isfield(atom,'bond')
    atom=rmfield(atom,'bond');
end

disp('Looking for neighbours/bonds')
Bond_index=[];b=1;
for i=1:length(XYZ_labels)
    k=0;j=1;
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
                    b=b+1;
                end
                
            end
        end
        j=j+1;
    end
    if mod(i,1000)==1
        if i > 1
            i-1
        end
    end
end


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
    assignin('caller','CoordNumber',CoordNumber);
    assignin('caller','Remove_ind',Remove_ind);
end

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

[atom.type]=atom.fftype;
nBonds=size(Bond_index,1);

assignin('caller','nBonds',nBonds);
assignin('caller','radius_limit',radius_limit);
assignin('caller','Bond_index',Bond_index);
assignin('caller','Neigh_index',Neigh_index);
% assignin('caller','bond_matrix',dist_matrix);
assignin('caller','dist_matrix',dist_matrix);

