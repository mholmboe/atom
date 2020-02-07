%% unreplicate_atom.m
% * This function 'condenses' a atom struct in the xyz-directions with the
% condensation factors given by the 1x3 vector condense_factors.
%
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = unreplicate_atom(atom,Box_dim,condense_factors)
% # atom = unreplicate_atom(atom,Box_dim,condense_factors,10)

function atom = unreplicate_atom(atom,Box_dim,condense_factors,varargin)

if nargin > 3
    ind=varargin{1};
    if ind>0
        atom=translate_atom(atom,-[[atom(ind(1)).x] [atom(ind(1)).y] [atom(ind(1)).z]]);
    end
end

atom = wrap_atom(atom,Box_dim);
% atom = slice_triclinic_atom(atom,Box_dim);

if size(Box_dim(1,:),2) > 3
    Lx = Box_dim(1);
    Ly = Box_dim(2);
    Lz = Box_dim(3);
    xy = Box_dim(6);
    xz = Box_dim(8);
    yz = Box_dim(9);
else
    Lx = Box_dim(1);
    Ly = Box_dim(2);
    Lz = Box_dim(3);
    xy = 0;
    xz = 0;
    yz = 0;
end

Box_dim=[Lx/condense_factors(1) Ly/condense_factors(2) Lz/condense_factors(3)...
    0 0 xy/condense_factors(1) 0 xz/condense_factors(2) yz/condense_factors(3)];

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
if size(Box_dim(1,:),2)==3
    XYZ_data(:,1) = XYZ_data(:,1) - Box_dim(1)*floor(XYZ_data(:,1)./Box_dim(1));
    XYZ_data(:,2) = XYZ_data(:,2) - Box_dim(2)*floor(XYZ_data(:,2)./Box_dim(2));
    XYZ_data(:,3) = XYZ_data(:,3) - Box_dim(3)*floor(XYZ_data(:,3)./Box_dim(3));
elseif size(Box_dim(1,:),2)==9
    XYZ_data(:,1) = XYZ_data(:,1) - xy*floor(XYZ_data(:,2)./Box_dim(2));
    XYZ_data(:,1) = XYZ_data(:,1) - xz*floor(XYZ_data(:,3)./Box_dim(3));
    XYZ_data(:,2) = XYZ_data(:,2) - yz*floor(XYZ_data(:,3)./Box_dim(3));
    XYZ_data(:,1) = XYZ_data(:,1) - Box_dim(1)*floor(XYZ_data(:,1)./Box_dim(1));
    XYZ_data(:,2) = XYZ_data(:,2) - Box_dim(2)*floor(XYZ_data(:,2)./Box_dim(2));
    XYZ_data(:,3) = XYZ_data(:,3) - Box_dim(3)*floor(XYZ_data(:,3)./Box_dim(3));
end

X=num2cell(XYZ_data(:,1));
Y=num2cell(XYZ_data(:,2));
Z=num2cell(XYZ_data(:,3));

[atom.x]=X{:};
[atom.y]=Y{:};
[atom.z]=Z{:};

if nargin > 4
    [atom.occupancy]=deal(1/(condense_factors(1)*condense_factors(2)*condense_factors(3)));
end

assignin('caller','Box_dim',Box_dim);

% XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
% if size(Box_dim(1,:),2)==3
%     XYZ_data(:,1) = floor(XYZ_data(:,1)./Box_dim(1));
%     XYZ_data(:,2) = floor(XYZ_data(:,2)./Box_dim(2));
%     XYZ_data(:,3) = floor(XYZ_data(:,3)./Box_dim(3));
% elseif size(Box_dim(1,:),2)==9
%     XYZ_data(:,1) = floor(XYZ_data(:,1)./Box_dim(1))+xy*floor(XYZ_data(:,2)./Box_dim(2));
%     XYZ_data(:,2) = floor(XYZ_data(:,2)./Box_dim(2))+xz*floor(XYZ_data(:,3)./Box_dim(3));
%     XYZ_data(:,3) = floor(XYZ_data(:,3)./Box_dim(3))+yz*floor(XYZ_data(:,3)./Box_dim(3));
% end
% 
% [Y,i] = sort(XYZ_data(:,1),1);
% XYZ_data = XYZ_data(i,:);
% XYZ = unique(XYZ_data,'rows','stable');
% 
% assignin('caller','rep_Y',Y);
% assignin('caller','rep_factors',XYZ_data);
% assignin('caller','rep_XYZ',XYZ);
% assignin('caller','rep_ind_i',i);
