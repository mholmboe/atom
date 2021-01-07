%% unreplicate_atom.m
% * This function 'condenses' a atom struct in the xyz-directions with the
% condensation factors given by the 1x3 vector condense_factors.
%
%
%% Version
% 2.09
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = unreplicate_atom(atom,Box_dim,condense_factors)
% # atom = unreplicate_atom(atom,Box_dim,condense_factors,10)
% # atom = unreplicate_atom(atom,Box_dim,condense_factors,0,'occupancy')

function atom = unreplicate_atom(atom,Box_dim,condense_factors,varargin)

if nargin > 3
    ind=varargin{1};
    if ind>0
        atom=translate_atom(atom,-[[atom(ind(1)).x] [atom(ind(1)).y] [atom(ind(1)).z]]);
    end
end

if size(Box_dim(1,:),2) > 3
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

orto=wrap_atom(orto,orto_Box_dim);

Box_dim=[orto_Box_dim(1)/condense_factors(1) orto_Box_dim(2)/condense_factors(2) orto_Box_dim(3)/condense_factors(3)];

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
XYZ_data(:,1) = XYZ_data(:,1) - Box_dim(1)*floor(XYZ_data(:,1)./Box_dim(1));
XYZ_data(:,2) = XYZ_data(:,2) - Box_dim(2)*floor(XYZ_data(:,2)./Box_dim(2));
XYZ_data(:,3) = XYZ_data(:,3) - Box_dim(3)*floor(XYZ_data(:,3)./Box_dim(3));

X=num2cell(XYZ_data(:,1));
Y=num2cell(XYZ_data(:,2));
Z=num2cell(XYZ_data(:,3));

[orto.x]=X{:};
[orto.y]=Y{:};
[orto.z]=Z{:};

atom = triclinic_atom(orto,Box_dim,[alfa beta gamma],'angles');

Box_dim=triclinic_Box_dim;

Box_dim(Box_dim<0.0001&Box_dim>-0.0001)=0;
if sum(Box_dim(4:end))== 0
    Box_dim=Box_dim(1:3);
    atom = wrap_atom(atom,Box_dim);
end

if nargin > 4
    [atom.occupancy]=deal(1/(condense_factors(1)*condense_factors(2)*condense_factors(3)));
end

assignin('caller','Box_dim',Box_dim);

