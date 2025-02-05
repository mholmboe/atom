%% remove_sametype_bond.m
% * This function removes covalent bonds, angles and dihedrals (latter two
% * optional) between atoms of the same type element, based on already
% * existing Bond | Angle | Dihedral_index variables.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=remove_sametype_bond(atom,Box_dim,Bond_index) % Basic input arguments
% # atom=remove_sametype_bond(atom,Box_dim,Bond_index,Angle_index) % Will also remove any angles
% # atom=remove_sametype_bond(atom,Box_dim,Bond_index,Angle_index,Dihedral_index) % Will also remove any dihedrals
% # atom=remove_sametype_bond(atom,Box_dim,Bond_index,Angle_index,Dihedral_index,rmaxlong) % Allows setting of the max cutoff

function atom = remove_sametype_bond(atom,Box_dim,varargin)

if nargin==2
    disp('You need to atleast supply the Bond_index variable')
elseif nargin==3
    count=1;
    Bond_index=varargin{1};
    Angle_index=[];
    Dihedral_index=[];
    rmaxlong=0.00001+max(Bond_index(:,3));
elseif nargin==4
    count=2;
    Bond_index=varargin{1};
    Angle_index=varargin{2};
    Dihedral_index=[];
    rmaxlong=0.00001+max(Bond_index(:,3));
elseif nargin==5
    count=3;
    Bond_index=varargin{1};
    Angle_index=varargin{2};
    Dihedral_index=varargin{3};
    rmaxlong=0.00001+max(Bond_index(:,3));
elseif nargin>5
    count=3;
    Bond_index=varargin{1};
    Angle_index=varargin{2};
    Dihedral_index=varargin{3};
    rmaxlong=varargin{4};
end

Atom_labels=unique([atom.type]);

for c=1:count
    index=[];
    if c==1
        index=Bond_index;
    elseif c==2
        index=Angle_index;
    elseif c==3
        index=Dihedral_index;
    end
    if length(index)>0
        all_rm_ind=[];rm_ind=[];
        for i=1:length(Atom_labels)
            ind_type=find(strncmpi([atom.type],Atom_labels(i),1));
            for j=1:c
                [type_row,type_col]=ind2sub(size(index(:,j:j+1)),find(ismember(index(:,j:j+1),ind_type)));
                unique_ind = unique(type_row);
                rm_ind=find(histc(type_row,unique_ind)-1);
                all_rm_ind=[all_rm_ind; rm_ind];
            end
        end
        index(unique(all_rm_ind),:)=[];
        if c==1
            Bond_index=index;
            assignin('caller','Bond_index',Bond_index);
        elseif c==2
            Angle_index=index;
            assignin('caller','Angle_index',Angle_index);
        elseif c==3
            Dihedral_index=index;
            assignin('caller','Dihedral_index',Dihedral_index);
        end
    end

end

atom=recalc_bond_atom(atom,Box_dim,Bond_index,Angle_index,Dihedral_index,rmaxlong);

end

