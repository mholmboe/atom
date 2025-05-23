%% protonate_atom.m
% * This function protonates the sites in the atom struct given by the
% index vector ind by adding a H's to a new H atom struct. It does so by
% placing the H opposite to the mean position of all neughbours within 2.5
% �ngstr�m of the site to be protonated
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # Hatom = protonate_atom(atom,Box_dim) % Protonating all O's that are only single bonded
% # Hatom = protonate_atom(atom,Box_dim,ind) % Protonates all sites with index ind
% # Hatom = protonate_atom(atom,Box_dim,ind,rmaxlong) % rcut can be used to change the default cutoff 2.25 �ngstr�m
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

% atom = element_atom(atom);

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

if nargin < 5
    heal_type={'H'};
else
    heal_type=varargin{3};
    if ~iscell(heal_type)
        heal_type={heal_type};
    end
end

atom = neigh_atom(atom,Box_dim,1.25,rmaxlong);

if numel(ind)<1
    i=1;
    while i<=size(atom,2)
        if strncmpi([atom(i).type],'O',1) && numel(atom(i).neigh.index) < 2
            ind=[ind i];
        end
        i=i+1;
    end
end

if numel(ind)<1
    prop=analyze_atom(atom,Box_dim);
    ind=heal_ind;
end

disp('Guessing this many H�s!')

if numel(ind) > 0
    ind=unique(ind);
    H_atom=[];
    for i=ind
        i
        atom(i).neigh.type
        Neigh_cell = sort([atom(i).neigh.type]);
        if isempty(Neigh_cell) > 0 && iscell(Neigh_cell)
            Neighbours=strcat(Neigh_cell{:});
        else
            Neighbours={'Nan'};
        end
        if numel(H_atom)==0
            H_atom = xyz2atom({'H'},[0 0 0],Box_dim,[atom(1).resname],[]);
            H_atom.molid=atom(end).molid;
            H_atom.index=1;
        else
            H_atom(size(H_atom,2)+1)=H_atom(end);
        end
        [H_atom(end).type]=heal_type;
        [H_atom(end).fftype]=heal_type;
        H_atom(end).index=size(H_atom,2);
        r_vec=atom(i).neigh.r_vec;
        H_coords=num2cell([atom(i).x atom(i).y (atom(i).z)]-0.9572*mean([r_vec(:,1) r_vec(:,2) r_vec(:,3)],1)/norm(mean([r_vec(:,1) r_vec(:,2) r_vec(:,3)],1)));
        [H_atom(end).x H_atom(end).y H_atom(end).z]=deal(H_coords{:});
    end

    if size(H_atom,2)>1
        dist_matrix=dist_matrix_atom(H_atom,Box_dim);
        i=1;rmind_tot=[];
        while i < size(H_atom,2)
            rmind=find(dist_matrix(:,i)<+.85);
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

                try
                    rmind_tot=[rmind_tot; rmind(rmind>i)];
                catch
                    rmind_tot=[rmind_tot; rmind(rmind>i)];
                end
            end
            i=i+1;
        end
        H_atom(rmind_tot)=[];
    end
end
disp('Created this many H�s!')
size(H_atom,2)

if isfield(H_atom,'element')
    H_atom=rmfield(H_atom,'element');
end

end