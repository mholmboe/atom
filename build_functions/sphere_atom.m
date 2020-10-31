%% sphere_atom.m
% * This function slices a spherical particle (like a colloid) of the atom
% * struct. If the spherical particle radius is larger then the atom structs
% * Box_dim dimensions in x/y/z-directions, the atom struct will be replicated
%
%
%% Version
% 2.08
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = sphere_atom(atom,Box_dim,30)
% # atom = sphere_atom(atom,Box_dim,30,[50 50 50]) % In order to set a new Box_dim
%
function atom = sphere_atom(atom,Box_dim,radius,varargin)

if min(Box_dim(1:3))<radius
    rep_factors=ceil((2*radius)./[max([atom.x]) max([atom.y]) max([atom.z])]);
    rep_factors+(~rem(rep_factors, 2) == 0); % To make sure its an even number
    atom = replicate_atom(atom,Box_dim,rep_factors);
end

if numel(Box_dim)==1
    Box_dim(1)=Box_dim(1);
    Box_dim(2)=Box_dim(1);
    Box_dim(3)=Box_dim(1);
    xy=0;
    xz=0;
    yz=0;
elseif numel(Box_dim)==3
    lx=Box_dim(1);
    ly=Box_dim(2);
    lz=Box_dim(3);
    xy=0;
    xz=0;
    yz=0;
elseif size(Box_dim,2)==9
    lx=Box_dim(1);
    ly=Box_dim(2);
    lz=Box_dim(3);
    xy=Box_dim(6);
    xz=Box_dim(8);
    yz=Box_dim(9);
end

% Triclininc Box parameters
a=lx;
b=(ly^2+xy^2)^.5;
c=(lz^2+xz^2+yz^2)^.5;
alfa=rad2deg(acos((ly*yz+xy*xz)/(b*c)));
beta=rad2deg(acos(xz/c));
gamma=rad2deg(acos(xy/b));

% Volume
v=(1 - cos(deg2rad(alfa))^2 - cos(deg2rad(beta))^2 - cos(deg2rad(gamma))^2 + 2*cos(deg2rad(alfa))*cos(deg2rad(beta))*cos(deg2rad(gamma)))^.5;

% Transformation matrix
FromFrac=[a b*cos(deg2rad(gamma)) c*cos(deg2rad(beta));...
    0 b*sin(deg2rad(gamma))  c*(cos(deg2rad(alfa))-cos(deg2rad(beta))*cos(deg2rad(gamma)))/sin(deg2rad(gamma));...
    0 0 c*v/sin(deg2rad(gamma))];
% The center of the Box_dim
Box_center=FromFrac*[0.5 0.5 0.5]'

distance = zeros(1,size(atom,2));
for i=1:size(atom,2)
    distance(i) = (([atom(i).x]-Box_center(1))^2 + ([atom(i).y]-Box_center(2))^2 + ([atom(i).z]-Box_center(3))^2)^0.5;
end
distance(distance<=radius)=1;
distance(distance>=radius)=0;
atom=update_atom(atom(logical(distance)));

if nargin>3
    Box_dim=varargin{1};
    atom=center_atom(atom,Box_dim)
end

assignin('caller','Box_dim',Box_dim)

