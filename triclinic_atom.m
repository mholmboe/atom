%% triclinic_atom.m
% * This function transforms an orthogonal atom struct to a triclinic with 
% * the angles alfa, beta, gamma or tilt factors xy, xz, yz
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se


%% Examples
% * atom = triclinic_atom(atom,Box_dim,[alfa beta gamma],'angle')
% * atom = triclinic_atom(atom,Box_dim,[xy xz yz],'tilt')


function atom = triclinic_atom(atom,Box_dim,angleparam,angletype)

if strncmpi(angletype,'angle',5)
    % Angle values
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
else
    % Tilt values
    
    lx=Box_dim(1);
    ly=Box_dim(2);
    lz=Box_dim(3);
%     if size(Box_dim,2)==9
%         xy=Box_dim(6);
%         xz=Box_dim(8);
%         yz=Box_dim(9);
%     end
    xy=angleparam(1);
    xz=angleparam(2);
    yz=angleparam(3);
    a=lx;
    b=(ly^2+xy^2)^.5;
    c=(lz^2+xz^2+yz^2)^.5;
    alfa=rad2deg(acos((ly*yz+xy*xz)/(b*c)));
    beta=rad2deg(acos(xz/c));
    gamma=rad2deg(acos(xy/b));
    %     lx = a;
    %     ly = (b^2-xy^2)^.5;
    %     lz = (c^2 - xz^2 - yz^2)^0.5;
    %     alfa=rad2deg(acos((ly*yz+xy*xz)/(b*c)))
    %     beta=rad2deg(acos(xz/c));
    %     gamma=rad2deg(acos(xy/b));
end

% Straight from wikipedia
% From fractional coordinates
% FromFrac*ToFrac

% Straight from wikipedia
% From fractional coordinates
% v=a*b*c*(1 - cos(deg2rad(alfa))^2 - cos(deg2rad(beta))^2 - cos(deg2rad(gamma))^2 + 2*cos(deg2rad(alfa))*cos(deg2rad(beta))*cos(deg2rad(gamma)))^.5;
v=(1 - cos(deg2rad(alfa))^2 - cos(deg2rad(beta))^2 - cos(deg2rad(gamma))^2 + 2*cos(deg2rad(alfa))*cos(deg2rad(beta))*cos(deg2rad(gamma)))^.5;

FromFrac=[a b*cos(deg2rad(gamma)) c*cos(deg2rad(beta));...
    0 b*sin(deg2rad(gamma))  c*(cos(deg2rad(alfa))-cos(deg2rad(beta))*cos(deg2rad(gamma)))/sin(deg2rad(gamma));...
    0 0 c*v/sin(deg2rad(gamma))];

% To fractional coordinates
ToFrac=[1/a -cos(deg2rad(gamma))/(a*sin(deg2rad(gamma))) (cos(deg2rad(alfa))*cos(deg2rad(gamma))-cos(deg2rad(beta)))/(a*v*sin(deg2rad(gamma)));...
    0 1/(b*sin(deg2rad(gamma)))  (cos(deg2rad(beta))*cos(deg2rad(gamma))-cos(deg2rad(alfa)))/(b*v*sin(deg2rad(gamma)));...
    0 0 sin(deg2rad(gamma))/(c*v)];

% FromFrac*ToFrac

XYZ_labels=[atom.type]';
XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
XYZ_data_frac=XYZ_data;XYZ_data_tric=XYZ_data;
for i=1:size(atom,2)
    XYZ_data_frac(i,:)=[XYZ_data(i,1)/lx XYZ_data(i,2)/ly XYZ_data(i,3)/lz]';
    XYZ_data_tric(i,:)=FromFrac*[XYZ_data_frac(i,1) XYZ_data_frac(i,2) XYZ_data_frac(i,3)]';
    atom(i).x=XYZ_data_tric(i,1);
    atom(i).y=XYZ_data_tric(i,2);
    atom(i).z=XYZ_data_tric(i,3);
end

Box_dim=[lx ly lz 0 0 xy 0 xz yz];

Box_dim(Box_dim<0.00001&Box_dim>-0.00001)=0;
if sum(Box_dim(4:end))== 0
    Box_dim=Box_dim(1:3);
end

assignin('caller','triclinic_Box_dim',Box_dim);


