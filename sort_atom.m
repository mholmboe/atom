%% sort_atom.m - This section orders to atoms with respect to z
function atom = sort_atom(atom,varargin)
%% sort_atom.m
% * This section orders to atoms with respect to z
% * atom is the atom struct
% * Tested when?
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * atom = sort_atom(atom)
% * atom = sort_atom(atom,'x') % Not working yet
% * atom = sort_atom(atom,'y') % Not working yet
% * atom = sort_atom(atom,'z')
% * atom = sort_atom(atom,'molid') % Not working yet
% * atom = sort_atom(atom,'resname') % Not working yet

if nargin>1
    mode=char(varargin{1});
    if strcmpi(mode,'resname')
        Resnames=unique([atom.resname]);
        for i=1:length(Resnames)
           eval(strcat('ind_',num2str(i),'=ismember([atom.resname],Resnames(i));'));
        end
        
        sorted_atom=[];
        for i=1:length(Resnames)
           eval(strcat('sorted_atom=[sorted_atom atom(ind_',num2str(i),')];')');
        end
        atom=update_atom(sorted_atom);
    end
else

atom = median_atom_func(atom)
sorted_coords=sortrows([[atom.index];[atom.molid];[atom.x];[atom.y];[atom.z];[atom.med_z]]',6);
nAtoms=size(atom,2);
Index=sorted_coords(:,1);
MolID=sorted_coords(:,2);
X_coord=sorted_coords(:,3);
Y_coord=sorted_coords(:,4);
Z_coord=sorted_coords(:,5);

nmol=1;first_in=[1];last_in=[];
for i=1:nAtoms;
    if i > 1 && MolID(i) ~= MolID(i-1)
        nmol=nmol+1;
        atom(i).molid=nmol;
        first_in(atom(i).molid,1)=i;last_in(atom(i).molid-1,1)=i-1;
    elseif i > 1 && MolID(i) == MolID(i-1)
        atom(i).molid=atom(i-1).molid;
    elseif i == 1
        atom(i).molid=1;
    end

    %atom(i).resname=Resname(i);
    %atom(i).type=Atomtype(i);
    %atom(i).fftype=Atomtype(i);
    %atom(i).index=mod(i,100000);
    %atom(i).neigh.type  = {};
    atom(i).neigh.index  = zeros(6,1);
    atom(i).neigh.dist  = zeros(6,1);
    atom(i).bond.type  = zeros(6,1);
    atom(i).bond.index  = zeros(6,1);
    atom(i).x=X_coord(i);
    atom(i).y=Y_coord(i);
    atom(i).z=Z_coord(i);
end
last_in(atom(end).molid,1)=nAtoms;

% for i=1:max([atom(:).molid])
%     ind=find([atom.molid]==i);
%     [atom(ind).ave_x]=deal(mean([atom(ind).x]));
%     [atom(ind).ave_y]=deal(mean([atom(ind).y]));
%     [atom(ind).ave_z]=deal(mean([atom(ind).z]));
% end

atom = update_atom(atom);

atom = median_atom(atom);

end