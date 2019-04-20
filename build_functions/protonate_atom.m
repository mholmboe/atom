%% protonate_atom.m
% * This function protonates the sites in the atom struct given by the 
% index vector ind by adding a H's to a new H atom struct. It does so by 
% placing the H opposite to the mean position of all neughbours within 2.5
% Ångström of the site to be protonated
%
%% Version
% 2.03
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # Hatom = protonate_atom(atom,Box_dim) % Protonating all O's that are only single bonded
% # Hatom = protonate_atom(atom,Box_dim,ind) % Protonates all sites with index ind
% # Hatom = protonate_atom(atom,Box_dim,ind,rmaxlong) % rcut can be used to change the default cutoff 2.5 Ångström
% # Hatom = protonate_atom(atom,Box_dim,ind,rmaxlong,{'He'}) % {'He'} can be used to change the default atomtype H to He
% # Hatom = protonate_atom(atom,Box_dim,ind,rmaxlong,{'He'},'minus') % 'minus' or default 'plus' denotes the tilt direction of the added H in the Z-direction
%
function H_atom = protonate_atom(atom,Box_dim,varargin)

if numel(Box_dim)==1
    Box_dim(1)=Box_dim(1);
    Box_dim(2)=Box_dim(1);
    Box_dim(3)=Box_dim(1);
    xy=0;xz=0;yz=0;
elseif numel(Box_dim)==3
    xy=0;xz=0;yz=0;
else
    xy=Box_dim(6);xz=Box_dim(8);yz=Box_dim(9);
end

if nargin==2
    disp('Assuming all oxygen atoms should have 2 neighbours...');
    disp('else also supply an ind vector for sites to protonate!');
end

atom = element_atom(atom);

if nargin > 2
    ind=varargin{1};
else
    ind=[];
end

if nargin > 3
    rmaxlong=varargin{2};
else
    rmaxlong=2.25;
end

if nargin < 4
    heal_type={'H'};
else
    heal_type=varargin{3};
    if ~iscell(heal_type)
        heal_type={heal_type};
    end
end

atom = neigh_atom(atom,Box_dim,1.25,rmaxlong);
    
if numel(ind)==0
    i=1; H_atom=[];
    while i<=size(atom,2)
        Neigh_cell = sort([atom(i).neigh.type]);
        if isempty(Neigh_cell) > 0 && iscell(Neigh_cell)
            Neighbours=strcat(Neigh_cell{:});
        else
            Neighbours={'Nan'};
        end
        % numel(atom(i).neigh.index)
        if strncmpi([atom(i).type],'O',1) && numel(atom(i).neigh.index) < 2
            disp('Adding H to O')
            if numel(H_atom)==0
                H_atom=atom(1);
            else
                H_atom(size(H_atom,2)+1)=H_atom(end);
            end
            [H_atom(end).type]=heal_type;[H_atom(end).fftype]=heal_type;
            H_atom(end).index=size(H_atom,2);
            r_vec=atom(i).neigh.r_vec;
            H_coords=num2cell([atom(i).x atom(i).y (atom(i).z)]-0.9572*mean([r_vec(:,1) r_vec(:,2) r_vec(:,3)],1)/norm(mean([r_vec(:,1) r_vec(:,2) r_vec(:,3)],1)));
            [H_atom(end).x H_atom(end).y H_atom(end).z]=deal(H_coords{:});
        end
        i=i+1;
    end
elseif nargin > 2
    ind=varargin{1};
    H_atom=[];
    for i=ind
        Neigh_cell = sort([atom(i).neigh.type]);
        if isempty(Neigh_cell) > 0 && iscell(Neigh_cell)
            Neighbours=strcat(Neigh_cell{:});
        else
            Neighbours={'Nan'};
        end
        if numel(H_atom)==0
            H_atom=atom(1);
        else
            H_atom(size(H_atom,2)+1)=H_atom(end);
        end
        [H_atom(end).type]=heal_type;[H_atom(end).fftype]=heal_type;
        H_atom(end).index=size(H_atom,2);
        r_vec=atom(i).neigh.r_vec;
        H_coords=num2cell([atom(i).x atom(i).y (atom(i).z)]-0.9572*mean([r_vec(:,1) r_vec(:,2) r_vec(:,3)],1)/norm(mean([r_vec(:,1) r_vec(:,2) r_vec(:,3)],1)));
        [H_atom(end).x H_atom(end).y H_atom(end).z]=deal(H_coords{:});
    end
end

dist_matrix=dist_matrix_atom(H_atom,Box_dim);

i=1;rmind_tot=[];
while i < size(H_atom,2)
    rmind=find(dist_matrix(:,i)<0.85);
    if numel(rmind)>1
        x1=[H_atom(i).x];
        y1=[H_atom(i).y];
        z1=[H_atom(i).z];

        H_atom(rmind) = translate_atom(H_atom(rmind),[Box_dim(1)/2-x1 Box_dim(2)/2-y1 Box_dim(3)/2-z1]);
        H_atom(rmind) = wrap_atom(H_atom(rmind),Box_dim);
        
        [H_atom(i).x]=mean([H_atom(rmind).x]);
        [H_atom(i).y]=mean([H_atom(rmind).y]);
        [H_atom(i).z]=mean([H_atom(rmind).z]);
        
        H_atom(rmind) = translate_atom(H_atom(rmind),[-Box_dim(1)/2+x1 -Box_dim(2)/2+y1 -Box_dim(3)/2+z1]);
        rmind_tot=[rmind_tot rmind(rmind>i)];
    end
    i=i+1;
end
H_atom(rmind_tot)=[];

if isstruct(H_atom)
    H_atom=rmfield(H_atom,'element');
end
