%% import_mclf_C8dispersion.m
% * This function imports the C8 dispersion terms from the
% mclf code. Not tested alot..
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_atom_ddec('XYZ_even_tempered_net_atomic_charges.xyz')
%
function [atom,Box_dim] = import_mclf_C8(varargin)

if nargin>0
    filename=varargin{1};
else
    filename='MCLF_C8_dispersion_coefficients.xyz';
end

fileID = fopen(filename,'r');
line1 = {fgets(fileID)};
line2 = {fgets(fileID)};
nAtoms=str2double(line1);
Box_string=strsplit(char(line2));

if strcmp(Box_string(1),'Nonperiodic')
    disp('Using dummy Box_dim info')
    Box_string(2)={'10'};
    Box_string(3)={'10'};
    Box_string(4)={'10'};
end

%     lx=Box_dim(1);
%     ly=Box_dim(2);
%     lz=Box_dim(3);
%     xy=Box_dim(6);
%     xz=Box_dim(8);
%     yz=Box_dim(9);
%
%     a=lx;
%     b=(ly^2+xy^2)^.5;
%     c=(lz^2+xz^2+yz^2)^.5;
%     alfa=rad2deg(acos((ly*yz+xy*xz)/(b*c)));
%     beta=rad2deg(acos(xz/c));
%     gamma=rad2deg(acos(xy/b));
%
%     Cell=[a b c alfa beta gamma];
%
% UNITCELL [ax, ay, az, bx, by, bz, cx, cy, cz] ?
%  CELL| Volume [angstrom^3]:                                            40.037000
%  CELL| Vector a [angstrom]:       3.126     0.000     0.000   |a| =     3.126000
%  CELL| Vector b [angstrom]:      -1.563     2.707     0.000   |b| =     3.126000
%  CELL| Vector c [angstrom]:       0.000     0.000     4.731   |c| =     4.731000
%  CELL| Angle (b,c), alpha [degree]:                                    90.000000
%  CELL| Angle (a,c), beta  [degree]:                                    90.000000
%  CELL| Angle (a,b), gamma [degree]:                                   120.000000
%  CELL| Requested initial symmetry:                                     TRICLINIC

if length(Box_string)<=5
    Box_dim=zeros(1,3);
    Box_dim(1)=str2double(Box_string(2)); % Lx,a
    Box_dim(2)=str2double(Box_string(3)); % Ly,b
    Box_dim(3)=str2double(Box_string(4)); % Lz,c
else
    Box_dim=zeros(1,9);
    Box_dim(1)=str2double(Box_string(11)); % Lx,a
    Box_dim(2)=str2double(Box_string(17)); % Ly,b
    Box_dim(3)=str2double(Box_string(23)); % Lz,c

    Box_dim(6)=str2double(Box_string(16)); % xy
    Box_dim(8)=str2double(Box_string(21)); % xz
    Box_dim(9)=str2double(Box_string(22)); % yz
end

%% Close the text file.
fclose(fileID);

filetempID = fopen(filename,'r');
line1 = {fgets(filetempID)};
line2 = {fgets(filetempID)};
for i=1:nAtoms
    line = fgetl(filetempID);
    XYZ_string=strsplit((line));
    XYZ_labels(i,1) = XYZ_string(1);
    X(i) = XYZ_string(2);
    Y(i) = XYZ_string(3);
    Z(i) = XYZ_string(4);
    mclf_C8(i) = XYZ_string(5);

end
fclose(filetempID);
mclf_C8=str2double(mclf_C8);
XYZ_data=[str2double(X)' str2double(Y)' str2double(Z)'];

for i=1:nAtoms
    atom(i).resname={'MOL'};
    atom(i).molid=1;
    atom(i).type        = XYZ_labels(i,1);
    atom(i).fftype  = XYZ_labels(i,1);
    atom(i).index=mod(i,100000);
    atom(i).neigh.type  = {};
    atom(i).neigh.index  = zeros(6,1);
    atom(i).neigh.dist  = zeros(6,1);
    atom(i).bond.type  = zeros(6,1);
    atom(i).bond.index  = zeros(6,1);
    atom(i).angle.type  = zeros(6,1);
    atom(i).angle.index  = zeros(6,1);
    atom(i).x = XYZ_data(i,1);
    atom(i).y = XYZ_data(i,2);
    atom(i).z = XYZ_data(i,3);
    atom(i).vx=0;
    atom(i).vy=0;
    atom(i).vz=0;
    atom(i).C8_aubhr8=mclf_C8(i);
    atom(i).C8_kJmolnm8=mclf_C8(i)*4.35974472220718E-18*6.022E+23/1000*(0.052917721^8);
end

% mclf_C8

Atom_labels=unique([atom.type]);

All_C8=[];
for i=1:size(Atom_labels,2)
    ind=ismember([atom.type],Atom_labels(i));
    Ave_C8_aubhr8=mean([atom(ind).C8_aubhr8]);
    Ave_C8_kJmolnm8=mean([atom(ind).C8_kJmolnm8]);
    All_C8=[All_C8; [Atom_labels(i) {Ave_C8_aubhr8} {Ave_C8_kJmolnm8}]];
end

writecell(All_C8,'All_Ave_C8.dat');

assignin('caller','All_C8',All_C8);

end



