%% orto_atom.m
% * This function transforms a triclinic atom struct to an orthogonal one
%
%% Version
% 2.11
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = orto_atom(atom,Box_dim)
% # atom = orto_atom(atom,Box_dim,[alfa beta gamma],'angle')
% # atom = orto_atom(atom,Box_dim,[xy xz yz],'tilt')

function atom = orto_atom(atom,Box_dim,varargin)

if nargin > 2
    angleparam=cell2mat(varargin(1));
    angletype=varargin(2);
    if strncmpi(angletype,'angle',5)
        a=Box_dim(1);
        b=Box_dim(2);
        c=Box_dim(3);
        % Angle values
        alfa=angleparam(1);
        beta=angleparam(2);
        gamma=angleparam(3);
        lx = a;
        xy = b * cos(deg2rad(gamma));
        ly = (b^2-xy^2)^.5;
        xz = c*cos(deg2rad(beta));
        yz = (b*c*cos(deg2rad(alfa))-xy*xz)/ly;
        lz = (c^2 - xz^2 - yz^2)^0.5;
    else
        % Tilt values
        a=Box_dim(1);
        b=Box_dim(2);
        c=Box_dim(3);
        xy=angleparam(1);
        xz=angleparam(2);
        yz=angleparam(3);
        lx = a;
        ly = (b^2-xy^2)^.5;
        lz = (c^2 - xz^2 - yz^2)^0.5;
        alfa=rad2deg(acos((ly*yz+xy*xz)/(b*c)))
        beta=rad2deg(acos(xz/c));
        gamma=rad2deg(acos(xy/b));
    end
elseif size(Box_dim,2) == 9
    lx=Box_dim(1);
    ly=Box_dim(2);
    lz=Box_dim(3);
    xy=Box_dim(6);
    xz=Box_dim(8);
    yz=Box_dim(9);
    
    a=lx;
    b=(ly^2+xy^2)^.5;
    c=(lz^2+xz^2+yz^2)^.5;
    alfa=rad2deg(acos((ly*yz+xy*xz)/(b*c)));
    beta=rad2deg(acos(xz/c));
    gamma=rad2deg(acos(xy/b));
elseif size(Box_dim,2) == 3
    lx=Box_dim(1);
    ly=Box_dim(2);
    lz=Box_dim(3);
    xy=0;
    xz=0;
    yz=0;
    
    a=lx;
    b=(ly^2+xy^2)^.5;
    c=(lz^2+xz^2+yz^2)^.5;
    alfa=rad2deg(acos((ly*yz+xy*xz)/(b*c)));
    beta=rad2deg(acos(xz/c));
    gamma=rad2deg(acos(xy/b));
else
    disp('Something is wrong with Box_dim, what are you trying to do?')
end

Box_dim(Box_dim<0.00001&Box_dim>-0.00001)=0;
if sum(Box_dim(4:end))== 0
    Box_dim=Box_dim(1:3);
end

% Straight from wikipedia
% From fractional coordinates
volume=a*b*c*(1 - cos(deg2rad(alfa))^2 - cos(deg2rad(beta))^2 - cos(deg2rad(gamma))^2 + 2*cos(deg2rad(alfa))*cos(deg2rad(beta))*cos(deg2rad(gamma)))^.5;
v=(1 - cos(deg2rad(alfa))^2 - cos(deg2rad(beta))^2 - cos(deg2rad(gamma))^2 + 2*cos(deg2rad(alfa))*cos(deg2rad(beta))*cos(deg2rad(gamma)))^.5;

FromFrac=[a b*cos(deg2rad(gamma)) c*cos(deg2rad(beta));...
    0 b*sin(deg2rad(gamma))  c*(cos(deg2rad(alfa))-cos(deg2rad(beta))*cos(deg2rad(gamma)))/sin(deg2rad(gamma));...
    0 0 c*v/sin(deg2rad(gamma))];

% To fractional coordinates
ToFrac=[1/a -cos(deg2rad(gamma))/(a*sin(deg2rad(gamma))) (cos(deg2rad(alfa))*cos(deg2rad(gamma))-cos(deg2rad(beta)))/(a*v*sin(deg2rad(gamma)));...
    0 1/(b*sin(deg2rad(gamma)))  (cos(deg2rad(beta))*cos(deg2rad(gamma))-cos(deg2rad(alfa)))/(b*v*sin(deg2rad(gamma)));...
    0 0 sin(deg2rad(gamma))/(c*v)];

% FromFrac*ToFrac

if size(atom,2)>0
    
    XYZ_labels=[atom.type]';
    XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
    XYZ_data_frac=XYZ_data;XYZ_data_orto=XYZ_data;
    for i=1:size(atom,2)
        XYZ_data_frac(i,:)=ToFrac*[XYZ_data(i,1) XYZ_data(i,2) XYZ_data(i,3)]';
        XYZ_data_orto(i,:)=[lx ly lz].*XYZ_data_frac(i,:);
        atom(i).x=round(XYZ_data_orto(i,1),4);
        atom(i).y=round(XYZ_data_orto(i,2),4);
        atom(i).z=round(XYZ_data_orto(i,3),4);
        atom(i).xfrac=round(XYZ_data_frac(i,1),4);
        atom(i).yfrac=round(XYZ_data_frac(i,2),4);
        atom(i).zfrac=round(XYZ_data_frac(i,3),4);
    end
    
end

Box_dim=[lx ly lz];

% lx
% ly
% lz
% a
% b
% c
% volume


assignin('caller','orto_Box_dim',Box_dim);
assignin('caller','Box_volume',a*b*c*v);


