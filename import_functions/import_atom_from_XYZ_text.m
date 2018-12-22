%// Read lines from input file
clear all
close all
filename='test1.txt';
outfilename='Na12H2O-2L-Mt-D-vdW-DF.gro';
resname='MMT';
numColumns=4;

%// Edit these lines with the cif data
% a = 13.7100959
% b = 14.0372098
% c = 50;
% alfa = 90; % Angles
% beta = 90;
% gamma = 90.049927

in_Box_dim=[5.297602060 0.001107473 -0.075521525;...
    0.003337687 9.181726505 0.010987329;...
-0.187064721 -0.000079513 15.912446655;...
    ];

lx=in_Box_dim(1);
ly=in_Box_dim(5);
lz=in_Box_dim(9);
xy=in_Box_dim(4);
xz=in_Box_dim(7);
yx=in_Box_dim(2);
yz=in_Box_dim(8);
zx=in_Box_dim(3);
zy=in_Box_dim(6);

% lx = a;
% xy = b * cos(deg2rad(gamma));
% ly = (b^2-xy^2)^.5;
% xz = c*cos(deg2rad(beta));
% yz = (b*c*cos(deg2rad(alfa))-xy*xz)/ly;
% lz = (c^2 - xz^2 - yz^2)^0.5;

Box_dim=[lx ly lz xy xz yx yz zx zy];
% Box_dim=[lx ly lz 0 0 xy 0 xy xz];

% Box_dim(Box_dim<0.00001&Box_dim>-0.00001)=0;
% if sum(Box_dim(4:end))== 0
%     Box_dim=Box_dim(1:3);
% end

File_as_char = fileread(filename);
File_cell=strsplit(File_as_char);

XYZ=reshape(File_cell(1:end),numColumns,[])';

XYZ_labels=XYZ(:,1);
XYZ_data=zeros(size(XYZ(:,numColumns-2:end),1),size(XYZ(:,numColumns-2:end),2));
XYZ_data=str2double(XYZ(:,numColumns-2:end));

atom = xyz2atom(XYZ_labels,XYZ_data,Box_dim,resname,'molid');
[atom.molid]=deal(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Straight from wikipedia
% % % From fractional coordinates
% % % v=a*b*c*(1 - cos(deg2rad(alpha))^2 - cos(deg2rad(beta))^2 - cos(deg2rad(gamma))^2 + 2*cos(deg2rad(alpha))*cos(deg2rad(beta))*cos(deg2rad(gamma)))^.5;
% % v=(1 - cos(deg2rad(alfa))^2 - cos(deg2rad(beta))^2 - cos(deg2rad(gamma))^2 + 2*cos(deg2rad(alfa))*cos(deg2rad(beta))*cos(deg2rad(gamma)))^.5;
% %
% % FromFrac=[a b*cos(deg2rad(gamma)) c*cos(deg2rad(beta));...
% %     0 b*sin(deg2rad(gamma))  c*(cos(deg2rad(alfa))-cos(deg2rad(beta))*cos(deg2rad(gamma)))/sin(deg2rad(gamma));...
% %     0 0 c*v/sin(deg2rad(gamma))];
% %
% % % % To fractional coordinates
% % % ToFrac=[1/a -cos(deg2rad(gamma))/(a*sin(deg2rad(gamma))) (cos(deg2rad(alpha))*cos(deg2rad(gamma))-cos(deg2rad(beta)))/(a*v*sin(deg2rad(gamma)));...
% % %     0 1/(b*sin(deg2rad(gamma)))  (cos(deg2rad(beta))*cos(deg2rad(gamma))-cos(deg2rad(alpha)))/(b*v*sin(deg2rad(gamma)));...
% % %     0 0 sin(deg2rad(gamma))/(c*v)];
% %
% % % XYZ_labels=[atom.type]';
% % % XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
% % XYZ_data_frac=XYZ_data;XYZ_data_tric=XYZ_data;
% % for i=1:size(atom,2)
% %     %     XYZ_data_frac(i,:)=[XYZ_data(i,1)/lx XYZ_data(i,2)/ly XYZ_data(i,3)/lz]';
% %     XYZ_data_tric(i,:)=FromFrac*[XYZ_data_frac(i,1) XYZ_data_frac(i,2) XYZ_data_frac(i,3)]';
% %     atom(i).x=XYZ_data_tric(i,1);
% %     atom(i).y=XYZ_data_tric(i,2);
% %     atom(i).z=XYZ_data_tric(i,3);
% % end

atom = rotate_atom(atom,Box_dim,[-rad2deg(asin(zx/lx)) 0 0]);
atom = rotate_atom(atom,Box_dim,[0 -rad2deg(asin(yx/lx)) 0]);
atom = rotate_atom(atom,Box_dim,[0 0 -rad2deg(asin(zy/ly))]);
atom = translate_atom(atom,-[atom(15).x atom(15).y atom(15).z]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename(regexp(filename,'.txt'))=[];
filename(regexp(filename,'.cif'))=[];
filename(regexp(filename,'.dat'))=[];
write_atom_gro(atom,Box_dim,outfilename)
% write_atom_pdb(atom,Box_dim,strcat('out_',filename,'.pdb'))
vmd(atom,Box_dim)
