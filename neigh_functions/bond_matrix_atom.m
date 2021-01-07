%% bond_matrix_atom.m
% * This function tries to assign all bonds to a bond_matrix
%
%% Version
% 2.082
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=bond_matrix_atom(atom,Box_dim)

function atom = bond_matrix_atom(atom,Box_dim,varargin)

if ~isfield(atom,'element')
    atom = element_atom(atom);
end

[atom.type]=atom.element;
[atom.fftype]=atom.element;

XYZ_labels=[atom.type]';
XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];

Radiiproperties=load('Revised_Shannon_radii.mat');
% atom=bond_valence_atom(atom,Box_dim,1.25,2.25);
% clayff_param(sort(unique([atom.type])),'SPC/E');

if size(atom,2)>5000
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
    XYZ_radii(ismember([atom.type],Atom_label(i)))=median(Radiiproperties.CrysRadii(ind))';
    XYZ_formalcharge(ismember([atom.type],Atom_label(i)))=median(Radiiproperties.OxState(ind))';
end

assignin('caller','XYZ_radii',XYZ_radii);
assignin('caller','XYZ_formalcharge',XYZ_formalcharge);

distance=1.15;
XYZ_radii(XYZ_radii==0)=distance;
radius_matrix=repmat(XYZ_radii,1,length(XYZ_radii));
radius_limit=(radius_matrix+radius_matrix')*distance;
dist_matrix(dist_matrix==0)=100;
bond_matrix=dist_matrix-radius_limit;
dist_matrix(dist_matrix==100)=0;
bond_matrix(bond_matrix>0)=0;
bond_matrix(bond_matrix<0)=1;
disp('Radii+Radii limits')
% unique(radius_limit)

atom=rmfield(atom,'neigh');
for i=1:length(XYZ_labels)
    k=0;j=1;
    bond_ind=find(bond_matrix(:,i));
%     XYZ_labels(i)
%     numel(bond_ind)
    while j <= numel(bond_ind) && k <= numel(bond_ind) %<= neigh %atom(i).neigh=[];
        if bond_matrix(bond_ind(j),i)==1
            if XYZ_formalcharge(i)*XYZ_formalcharge(bond_ind(j))<0
                k=k+1;
                [atom(i).neigh.dist(k)]=dist_matrix(bond_ind(j),i);
                [atom(i).neigh.index(k)]=bond_ind(j);
                [atom(i).neigh.type(k)]=XYZ_labels(bond_ind(j));
                [atom(i).neigh.coords(k,:)]=[XYZ_data(bond_ind(j),1) XYZ_data(bond_ind(j),2) XYZ_data(bond_ind(j),3)];
                [atom(i).neigh.r_vec(k,:)]=[X_dist(bond_ind(j),i) Y_dist(bond_ind(j),i) Z_dist(bond_ind(j),i)];
            end
        end
        j=j+1;
    end
end

ind=find(tril(bond_matrix));
r=dist_matrix(ind);
[i,j] = ind2sub(size(bond_matrix),ind);

Bond_index = [j i r];
[Y,i] = sort(Bond_index(:,1));
Bond_index = Bond_index(i,:);
Bond_index = unique(Bond_index,'rows','stable');

rm_ind=find(Bond_index(:,3)>2.25);
Bond_index(rm_ind,:)=[];

assignin('caller','Radius_limit',radius_limit);
assignin('caller','Bond_matrix',bond_matrix);
assignin('caller','Bond_index2',Bond_index);
assignin('caller','Dist_matrix',dist_matrix);

