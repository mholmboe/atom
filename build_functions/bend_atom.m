%% ben_atom.m
% * This simple function tries to bend the coordinates in  an atom struct... 
% downwards. It works best if the original atom struct initally is centered.
%
%% Version
% 2.09
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = bend_atom(atom,Box_dim,Radii)

function atom = bend_atom(atom,Box_dim,radius,varargin)

trans=translate_atom(atom,-Box_dim./2);
  
x=[trans.x];
z=[trans.z];
y=[trans.y];
 
zr=(z-(radius-x.^2)./radius);
yr=y.*((radius-(z-min(zr)))/radius).^.5;
xr=x.*((radius-(z-min(zr)))/radius).^.5;

X=num2cell(xr);
Y=num2cell(yr);
Z=num2cell(zr);
[trans.x]=X{:};
[trans.y]=Y{:};
[trans.z]=Z{:};

atom=translate_atom(trans,Box_dim./2);

%plot3(xr,yr,zr,'.k')