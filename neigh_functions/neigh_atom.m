%% neigh_atom.m
% * This function checks which neighbors each atom has and ouputs their info
% * Todo... check if support for triclinic Box_dim works, because its untested...
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = neigh_atom(atom,Box_dim,rmax) % Basic input arguments
% # atom = neigh_atom(atom,Box_dim,rmax) % Allows settng the max cutoff
% # atom = neigh_atom(atom,Box_dim,rmax,101) % To analyze selected indexes
%
function atom = neigh_atom(atom,Box_dim,varargin)

nAtoms_ind=1:size(atom,2);

orig_index=[atom.index];

% skip_ind=[];
if nargin==2
    % Simple way to set the rmax for each atom
    rdist = 2.25*ones(size(atom,2),1);
    % Special thing for H's and water oxygens
    rdist(strncmpi([atom.type],'H',1))=1.25;
    rdist(strncmpi([atom.type],'Ow',2))=1.25;
elseif nargin==3
    % Simple way to set the rmax for each atom
    rdist = varargin{1}*ones(size(atom,2),1);
    % Special thing for H's and water oxygens
    rdist(strncmpi([atom.type],'H',1))=1.25;
    rdist(strncmpi([atom.type],'Ow',2))=1.25;
elseif nargin>3
    % Simple way to set the rmax for each atom
    rdist = varargin{2}*ones(size(atom,2),1);
    % Special thing for H's and water oxygens
    rdist(strncmpi([atom.type],'H',1))=varargin{1};
    rdist(strncmpi([atom.type],'Ow',2))=varargin{1};
end

if nargin>4
    nAtoms_ind=varargin{3};
end

XYZ_data=single([[atom.x]' [atom.y]' [atom.z]']);

for i=nAtoms_ind
    
    XYZ_solute=XYZ_data(i,:);
    
    rx=zeros(size(XYZ_data,1),1);ry=zeros(size(XYZ_data,1),1);rz=zeros(size(XYZ_data,1),1);
    
    lx=Box_dim(1);ly=Box_dim(2);lz=Box_dim(3);
    if numel(Box_dim)==3
        xy=0;xz=0;yz=0;
    else
        xy=Box_dim(6);xz=Box_dim(8);yz=Box_dim(9);
    end
    
    % Calculate Distance Components
    rx = XYZ_data(:,1) - XYZ_solute(1);
    ry = XYZ_data(:,2) - XYZ_solute(2);
    rz = XYZ_data(:,3) - XYZ_solute(3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z_gt_ind=find(rz > lz/2);z_lt_ind=find(rz < - lz/2);
    rz(z_gt_ind) = rz(z_gt_ind) - lz;
    rz(z_lt_ind) = rz(z_lt_ind) + lz;
    rx(z_gt_ind) = rx(z_gt_ind) - xz;
    rx(z_lt_ind) = rx(z_lt_ind) + xz;
    ry(z_gt_ind) = ry(z_gt_ind) - yz;
    ry(z_lt_ind) = ry(z_lt_ind) + yz;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_gt_ind=find(ry > ly/2);y_lt_ind=find(ry < - ly/2);
    ry(y_gt_ind) = ry(y_gt_ind) - ly;
    ry(y_lt_ind) = ry(y_lt_ind) + ly;
    rx(y_gt_ind) = rx(y_gt_ind) - xy;
    rx(y_lt_ind) = rx(y_lt_ind) + xy;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_gt_ind=find(rx > lx/2);x_lt_ind=find(rx < - lx/2);
    rx(x_gt_ind) = rx(x_gt_ind) - lx;
    rx(x_lt_ind) = rx(x_lt_ind) + lx;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dist = sqrt( rx(:,1).^2 + ry(:,1).^2 + rz(:,1).^2 ); % distance calc.
    
    % Find points inside circle
    in=intersect(find(dist>0),find(dist<rdist(i)));
    % in=sort(unique(in));
    % in = in(find(in~=i));
    % in = in(in~=i&~ismember(in,skip_ind));
    % neigh(i).in = [in(find(in~=i))];
    % atom(i).neigh.index = in;
    atom(i).neigh.index = orig_index(in)'; % New in version 2.04
    atom(i).neigh.type = deal([atom(in).type])';
    atom(i).neigh.dist = [dist(in,1)];
    atom(i).neigh.coords = [XYZ_data(in,1) XYZ_data(in,2) XYZ_data(in,3)];
    atom(i).neigh.r_vec = [rx(in) ry(in) rz(in)];
    
    if mod(i,1000)== 0
        i
    end
    
end

assignin('caller','dist',[atom(1).neigh.dist]);

end