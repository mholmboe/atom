%% import_cp2k.m
% * This function imports the type of .restart files that the package CP2K
% uses.
%
%% Version
% 2.11
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% #  atom = import_cp2k(filename)
%
function atom = import_cp2k(varargin)

if nargin>0
    filename=varargin{1};
else
    filename='MIN-1.restart';
end

% Open and read the input file
inputfile = fopen(filename, 'r');
C = textscan(inputfile, '%s', 'Delimiter', '\n');
fclose(inputfile);

nRows = size(C{1,1},1);
nColumns=size(strsplit(char(C{1,1}(end-1,:))),2);

% Search a specific string and find all rows containing matches
Cell = strfind(C{1}, '&SUBSYS');
Cell_rows = find(~cellfun('isempty', Cell));
Cell_data=C{1,1}(Cell_rows+2:Cell_rows+4);
% Cell_data = regexprep(Cell_data,'point','');
% Cell_data = regexprep(Cell_data,'DG','');
% Cell_data = regexprep(Cell_data,',','');
% Cell_data = regexprep(Cell_data,'+/-','');

CellData=[];
for i=1:size(Cell_data,1)
    temp=regexp(Cell_data{i,:}, '\s+', 'split');
    temp_data = cellfun(@str2double,temp);
    CellData=[CellData;temp_data(1,2:end)];
end

a=norm(CellData(1,:));
b=norm(CellData(2,:));
c=norm(CellData(3,:));
alpha=rad2deg(atan2(norm(cross(CellData(2,:),CellData(3,:))),dot(CellData(2,:),CellData(3,:))));
beta=rad2deg(atan2(norm(cross(CellData(1,:),CellData(3,:))),dot(CellData(1,:),CellData(3,:))));
gamma=rad2deg(atan2(norm(cross(CellData(1,:),CellData(2,:))),dot(CellData(1,:),CellData(2,:))));
Cell=[a b c alpha beta gamma];
Box_dim = Cell2Box_dim(Cell);

Coord = strfind(C{1}, '&COORD');
Coord_row = find(~cellfun('isempty', Coord));
CoordEnd = strfind(C{1}, '&END COORD');
CoordEnd_row = find(~cellfun('isempty', CoordEnd));
Coord_data=C{1,1}(Coord_row+1:CoordEnd_row-1);

nAtoms=size(Coord_data,1);
for i=1:nAtoms
    XYZ_string=strsplit(Coord_data{i,1});
    XYZ_labels(i,1) = XYZ_string(1);
    X(i) = XYZ_string(2);
    Y(i) = XYZ_string(3);
    Z(i) = XYZ_string(4);
end

XYZ_data=[str2double(X)' str2double(Y)' str2double(Z)'];

for i=1:nAtoms
    atom(i).resname = {'MOL'};
    atom(i).molid = 1;
    atom(i).type = XYZ_labels(i,1);
    atom(i).fftype = XYZ_labels(i,1);
    atom(i).charge = 0;
    atom(i).index = mod(i,100000);
    atom(i).neigh.type = {};
    atom(i).neigh.index = zeros(6,1);
    atom(i).neigh.dist = zeros(6,1);
    atom(i).bond.type = zeros(6,1);
    atom(i).bond.index = zeros(6,1);
    atom(i).angle.type = zeros(6,1);
    atom(i).angle.index = zeros(6,1);
    atom(i).x = XYZ_data(i,1);
    atom(i).y = XYZ_data(i,2);
    atom(i).z = XYZ_data(i,3);
    atom(i).vx = 0;
    atom(i).vy = 0;
    atom(i).vz = 0;
end

if min(Box_dim(1:3))<2.5
    nAtoms_init=size(atom,2);
  
    Box_dim_init=Box_dim;
    rep_atom=replicate_atom(atom,Box_dim,[2 2 2]);
    prop=analyze_atom(rep_atom,Box_dim);
    prop=prop(1:nAtoms_init);
    diff_valence=diff_valence(1:nAtoms_init);
    Box_dim=Box_dim_init;

else
    prop=analyze_atom(atom,Box_dim);
end

atom=bond_angle_atom(atom,Box_dim);

assignin('caller','CellMatrix',CellData);
assignin('caller','Cell',Cell);
assignin('caller','Box_dim',Box_dim);
assignin('caller','XYZ_labels',XYZ_labels);
assignin('caller','XYZ_data',XYZ_data);
assignin('caller','Angle_index',Angle_index);
assignin('caller','Bond_index',Bond_index);

assignin('caller','Tot_valence',Tot_valence);
assignin('caller','Tot_valence_oxstate',Tot_valence_oxstate);
assignin('caller','GII',GII);
assignin('caller','GII_noH',GII_noH);
assignin('caller','BondSummary',BondSummary);

assignin('caller','diff_valence',diff_valence);
assignin('caller','prop_atom',element);

try
    assignin('caller','dist_matrix',dist_matrix);
    assignin('caller','diff_bond_bv',[properties.rdiffvalence]');
catch
    
end

