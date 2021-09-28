%% import_atom_ddec.m
% * This function imports the partial net atomic charges (NAC's) from the
% chargemol/DDEC6 code. Not tested alot..
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_atom_ddec('XYZ6_even_tempered_net_atomic_charges.xyz')
%
function atom = import_ddec_charges(varargin)

if nargin>0
    filename=varargin{1};
else
    filename='DDEC6_even_tempered_net_atomic_charges.xyz';
end

fileID = fopen(filename,'r');
line1 = {fgets(fileID)};
line2 = {fgets(fileID)};
nAtoms=str2double(line1);
Box_string=strsplit(char(line2));
Box_dim=zeros(1,9);
Box_dim(1)=str2double(Box_string(11)); % Lx,a
Box_dim(2)=str2double(Box_string(17)); % Ly,b
Box_dim(3)=str2double(Box_string(23)); % Lz,c

Box_dim(6)=str2double(Box_string(16)); % xy
Box_dim(8)=str2double(Box_string(21)); % xz
Box_dim(9)=str2double(Box_string(22)); % yz

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
    ddec_charge(i) = XYZ_string(5);
    
end
fclose(filetempID);
ddec_charge=str2double(ddec_charge);
XYZ_data=[str2double(X)' str2double(Y)' str2double(Z)'];

% size(XYZ_labels)
% size(XYZ_data)
% XYZ_labels = dataArray{:, 1};
% XYZ_data = [dataArray{:,2} dataArray{:,3} dataArray{:,4}];

for i=1:nAtoms
    atom(i).resname={'MOL'};
    atom(i).molid=1;
    atom(i).type        = XYZ_labels(i,1);
    atom(i).fftype  = XYZ_labels(i,1);
    atom(i).charge      = 0;
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
    atom(i).charge=ddec_charge(i);
end

% atom = resname_atom(atom);

write_atom_pqr(atom,Box_dim,'ddec.pqr')

assignin('caller','Box_dim',Box_dim)
assignin('caller','XYZ_labels',XYZ_labels)
assignin('caller','XYZ_data',XYZ_data)

