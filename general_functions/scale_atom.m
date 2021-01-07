%% scale_atom.m
% * This function scales the coordinates in the atom struct
%
%% Version
% 2.09
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = scale_atom(atom,Box_dim,[1 1 1.5],'all')
% # atom = scale_atom(atom,Box_dim,0.9 0.9 0.9,'SOL')
%
function atom = scale_atom(atom,Box_dim,scale_vec,varargin)

if nargin>3
    Resname=varargin{1};
else
    Resname='ALL';
end

% disp('Scaling the coordinates')

nAtoms=size([atom.x],2);

if strncmpi(Resname,'ALL',3)
    ind_resname=1:nAtoms;
elseif strncmpi(Resname,'SYSTEM',3)
    ind_resname=1:nAtoms;
else
    ind_resname=find(strcmp([atom.resname],Resname));
end

if size(Box_dim,2) > 3
    lx = Box_dim(1);
    ly = Box_dim(2);
    lz = Box_dim(3);
    xy = Box_dim(6);
    xz = Box_dim(8);
    yz = Box_dim(9);
else
    lx = Box_dim(1);
    ly = Box_dim(2);
    lz = Box_dim(3);
    xy = 0;
    xz = 0;
    yz = 0;
end

a=lx;
b=(ly^2+xy^2)^.5;
c=(lz^2+xz^2+yz^2)^.5;
alfa=rad2deg(acos((ly*yz+xy*xz)/(b*c)));
beta=rad2deg(acos(xz/c));
gamma=rad2deg(acos(xy/b));

orto=orto_atom(atom,Box_dim);

Box_dim=[orto_Box_dim(1)*scale_vec(1) orto_Box_dim(2)*scale_vec(2) orto_Box_dim(3)*scale_vec(3)];

x_shift=num2cell([orto(ind_resname).x]*scale_vec(1)); [orto(ind_resname).x]=deal(x_shift{:});
y_shift=num2cell([orto(ind_resname).y]*scale_vec(2)); [orto(ind_resname).y]=deal(y_shift{:});
z_shift=num2cell([orto(ind_resname).z]*scale_vec(3)); [orto(ind_resname).z]=deal(z_shift{:});

atom = triclinic_atom(orto,Box_dim,[alfa beta gamma],'angles');

Box_dim=triclinic_Box_dim;

Box_dim(Box_dim<0.0001&Box_dim>-0.0001)=0;
if sum(Box_dim(4:end))== 0
    Box_dim=Box_dim(1:3);
end

atom=update_atom(atom);

assignin('caller','Box_dim',Box_dim);
