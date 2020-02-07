%% heal_atom.m
% * This function heals sites in the atom struct given by the index
% vector ind, by adding a certain atomtype to a new atom struct called
% healed_atom. It does so by placing the new atom type opposite to the
% mean position of all neighbours within rcut [Å] of the healed site.
% * Note that you first need to find which sites that should be healed,
% using for instance the bond_valende_atom function, and then decide with
% what atomtypes the unsaturated sites should be healed with.
%
%% Similar
% fuse_atom
% protonate_atom
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # healed_atom = heal_atom(atom,Box_dim,[6 16 26 36])
% # healed_atom = heal_atom(atom,Box_dim,[6:10:960],3)
% # healed_atom = heal_atom(atom,Box_dim,[3 202 493],3,'He')
% # healed_atom = heal_atom(atom,Box_dim,[3 202 493],2.5,'He')
% # healed_atom = heal_atom(atom,Box_dim,[3 202 493],2.5,'He',1.87)
% # healed_atom = heal_atom(atom,Box_dim,[3 202 493],2.5,'He','Shannon')

function healed_atom = heal_atom(atom,Box_dim,ind,varargin) % optional varargs == rmaxlong,heal_type,r|'Shannon');

if nargin==2
    disp('You must supply an index vector indicating which sites that should be healed');
    pause
end
if nargin==3
    rmaxlong=2.25;
else
    rmaxlong=varargin{1};
end
if nargin<4
    heal_type={'H'};
else
    heal_type=varargin{2};
    if ~iscell(heal_type)
        heal_type={heal_type};
    end
end

element = element_atom(atom);
element = neigh_atom(element,Box_dim,1.25,rmaxlong);

bonddist=ones(1,numel(ind));
if nargin==6
    rtype=varargin{3};
    if ~isnumeric(rtype)
        load('Revised_Shannon_radii.mat') % To heal metal sites with Oxygen atoms
        n=1;
        for i=ind
            Ion_M_ind=find(strcmp([element(i).type],Ion));
            if numel(Ion_M_ind)==0
                Ion_M_ind=find(strncmpi([element(i).type],Ion,1));
            end
            CN_M=numel(element(i).neigh.index)+1; % +1 since we have not added the extra ligand yet
            CN_M_ind=find(CN_M==CN);
            ind_M=intersect(Ion_M_ind,CN_M_ind);
            crysradii_M=CrysRadii(ind_M);
            
            Ion_O_ind=find(strcmp(Ion,'O'));
            if CN_M==1
                CN_M=2;
            elseif CN_M==5
                CN_M=6;
            elseif CN_M==7 || CN_M>8
                CN_M=8;
            end
            CN_O_ind=find(CN_M==CN);
            ind_O=intersect(Ion_O_ind,CN_O_ind);
            crysradii_O=CrysRadii(ind_O);
            bonddist(n)=crysradii_M+crysradii_O
            n=n+1;
        end
    else
        bonddist=rtype*ones(numel(ind)); % Set the bond-bond distance manually
    end
else
    danglingtype_radii = radius_ion([element(ind(1)).type]);
    if ~strncmpi(heal_type,'H',1)
        healtype_radii = radius_ion(heal_type);
        bonddist=0.8*(danglingtype_radii+healtype_radii)*ones(numel(ind));
        if bonddist(1) < 1.5
            bonddist=1.5*ones(numel(ind));
        end
    else
        bonddist=0.95*ones(numel(ind)); % like in O-H
    end
end

healed_atom=[];n=1;
for i=ind
    Neigh_cell = sort([element(i).neigh.type]);
    if isempty(Neigh_cell) > 0 && iscell(Neigh_cell)
        Neighbours=strcat(Neigh_cell{:});
    else
        Neighbours={'Nan'};
    end
    if numel(healed_atom)==0
        healed_atom=element(1);
    else
        healed_atom(size(healed_atom,2)+1)=healed_atom(end);
    end
    [healed_atom(end).type]=heal_type;[healed_atom(end).fftype]=heal_type;
    healed_atom(end).index=size(healed_atom,2);
    r_vec=element(i).neigh.r_vec;
    XYZ_data=num2cell([element(i).x element(i).y (element(i).z)]-bonddist(n)*mean([r_vec(:,1) r_vec(:,2) r_vec(:,3)],1)/norm(mean([r_vec(:,1) r_vec(:,2) r_vec(:,3)],1)));
    [healed_atom(end).x,healed_atom(end).y,healed_atom(end).z]=deal(XYZ_data{:});
    n=n+1;
end

dist_matrix=dist_matrix_atom(healed_atom,Box_dim);

i=1;rmind_tot=[];
while i < size(healed_atom,2)
    rmind=find(dist_matrix(:,i)<rmax);
    if numel(rmind)>1
        x1=[healed_atom(i).x];
        y1=[healed_atom(i).y];
        z1=[healed_atom(i).z];

        healed_atom(rmind) = translate_atom(healed_atom(rmind),[Box_dim(1)/2-x1 Box_dim(2)/2-y1 Box_dim(3)/2-z1]);
        healed_atom(rmind) = wrap_atom(healed_atom(rmind),Box_dim);
        
        [healed_atom(i).x]=mean([healed_atom(rmind).x]);
        [healed_atom(i).y]=mean([healed_atom(rmind).y]);
        [healed_atom(i).z]=mean([healed_atom(rmind).z]);
        
        healed_atom(rmind) = translate_atom(healed_atom(rmind),[-Box_dim(1)/2+x1 -Box_dim(2)/2+y1 -Box_dim(3)/2+z1]);
        rmind_tot=[rmind_tot rmind(rmind>i)];
    end
    i=i+1;
end
healed_atom(rmind_tot)=[];

if isstruct(healed_atom)
    try
        if ~isfield(atom,'element')
            healed_atom=rmfield(healed_atom,'element');
        end
    catch
    end
    
%     try
%         if isfield(healed_atom,'xfrac')
%             healed_atom=rmfield(healed_atom,'xfrac');
%             healed_atom=rmfield(healed_atom,'yfrac');
%             healed_atom=rmfield(healed_atom,'zfrac');
%         end
%     catch
%     end
end
