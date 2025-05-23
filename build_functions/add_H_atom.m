%% add_H_atom.m
% * This function protonates the sites in the atom struct given by the
% index vector ind by adding a H's to a new H atom struct. It does so by
% placing the H opposite to the mean position of all neughbours within rmaxlong
% Ångström of the site to be protonated. H_atom can be used with the function call
% [SOL,atom]=find_H2O(atom,Box_dim,1.05) to add the Hw's to the Ow's... 
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # Hatom = add_H_atom(atom,Box_dim,ind) % Protonates all sites with index ind
% # Hatom = add_H_atom(atom,Box_dim,ind,nH) % nH = 1 or 2
% # Hatom = add_H_atom(atom,Box_dim,ind,nH,rmaxlong) % rcut can be used to change the default cutoff 2.25 Ångström
% # Hatom = add_H_atom(atom,Box_dim,ind,nH,rmaxlong,{'He'}) % {'He'} can be used to change the default atomtype H to He
%

function H_atom = add_H_atom(atom,Box_dim,varargin)

% Default parameters
rOH    = 0.957;
HOHangle  = 104.5 * pi / 180;

if nargin==2
    disp('Assuming all oxygen atoms should have 2 neighbours...');
    disp('else also supply an ind vector for sites to protonate!');
end

if nargin > 2
    ind=varargin{1};
else
    ind=[];
end

if nargin > 3
    nH=varargin{2};
else
    nH=2;
end

if nargin > 4
    rmaxlong=varargin{3};
else
    rmaxlong=2.25;
end

if nargin < 6
    H_type={'H'};
else
    H_type=varargin{4};
    if ~iscell(H_type)
        H_type={H_type};
    end
end

% Build neighbor lists
atom = neigh_atom(atom,Box_dim,1.25,rmaxlong);

H_atom = [];index=atom(end).index;
for i = reshape(ind,1,[])

    % Compute bisector direction
    r_vecs = atom(i).neigh.r_vec;                     % Ni×3
    mean_dir = mean(r_vecs,1);
    u = -mean_dir / norm(mean_dir);                    % unit direction

    % Choose placement directions
    if nH==1
        direction = u;
    elseif nH==2
        % construct perpendicular basis
        tmp = [1 0 0];
        if abs(dot(tmp,u))>0.9, tmp = [0 1 0]; end
        v1 = cross(u, tmp); v1 = v1/norm(v1);
        theta = HOHangle/2;
        direction = [cos(theta)*u + sin(theta)*v1;
                     cos(theta)*u - sin(theta)*v1];
    else
        error('Only 1 or 2 H supported.');
    end

    % Create H atoms
    basePos = [atom(i).x, atom(i).y, atom(i).z];
    for k = 1:size(direction,1)
        pos = basePos + rOH * direction(k,:)
        Hn = xyz2atom(H_type, pos,Box_dim,[atom(1).resname],[]);
        Hn.molid = atom(i).molid;
        Hn.index = index + 1;
        H_atom = [H_atom, Hn];
        index = index + 1;
    end
end

if isfield(H_atom,'element')
    H_atom=rmfield(H_atom,'element');
end

end
