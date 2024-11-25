%% frac2atom.m
% * This function transforms fractional coordinates to cartesian
% * Its assumed the fractional coordinates are stored in the atom.x/atom.y/atom.z
% * Box_dim(1:3) should contain the real a,b,c unit cell lengths
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = frac2atom(atom,Box_dim(1:3),[82 93.1212 104.3232],'angle')
% # atom = frac2atom(atom,Box_dim(1:3),[0.0232 0.002131 0.0],'tiltfactors')
%
function atom = frac2atom(atom,Box_dim,angleparam,angletype)

a=Box_dim(1);
b=Box_dim(2);
c=Box_dim(3);

alfa=angleparam(1);
beta=angleparam(2);
gamma=angleparam(3);

lx = a;
xy = b * cos(deg2rad(gamma));
ly = (b^2-xy^2)^.5;
xz = c*cos(deg2rad(beta));
yz = (b*c*cos(deg2rad(alfa))-xy*xz)/ly;
lz = (c^2 - xz^2 - yz^2)^0.5;

Box_dim(Box_dim<0.00001&Box_dim>-0.00001)=0;
if sum(Box_dim(4:end))== 0
    Box_dim=Box_dim(1:3);
end

% atom_orto = orto_atom(atom,Box_dim,angleparam,angletype);

X=num2cell([atom.x]*lx);[atom(:).x]=deal(X{:});
Y=num2cell([atom.y]*ly);[atom(:).y]=deal(Y{:});
Z=num2cell([atom.z]*lz);[atom(:).z]=deal(Z{:});
%
atom = triclinic_atom(atom,[a b c],[alfa beta gamma],'angle');

% Box_dim(1:3)=[lx ly lz];

assignin('caller','Box_dim',Box_dim);

