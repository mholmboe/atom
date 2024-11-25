%% import_atom_poscar.m
% * This function imports an .poscar file into the atom struct
% * It tries to guess the Box_dim, so watch out!
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_atom_poscar('POSCAR')
%
function atom = import_atom_poscar(filename)

if regexp(filename,'POSCAR') ~= false
    filename = filename;
else

    filename = strcat(filename,'.poscar');
end


%% Open the text file.
FID = fopen(filename,'r');

All_data = textscan(FID, '%s', 'delimiter', '\n', 'whitespace', '');

%% Close the text file.
fclose(FID);

DataRow = strfind(All_data{1}, 'Cartesian');
startRow = find(~cellfun('isempty', DataRow), 1,'first')+1;
numElements=startRow-7;

fileID = fopen(filename,'r');
% fileID = fopen(filename,'r');
filename_in = {fgets(fileID)};
Repfactor = {fgets(fileID)};
Boxline1 = {fgets(fileID)};
Boxline2 = {fgets(fileID)};
Boxline3 = {fgets(fileID)};
Box_string1=strsplit(char(Boxline1));
Box_string2=strsplit(char(Boxline2));
Box_string3=strsplit(char(Boxline3));
Box_vec1 = str2double(Box_string1);
Box_vec2 = str2double(Box_string2);
Box_vec3 = str2double(Box_string3);
Box_vec1=Box_vec1(~isnan(Box_vec1));
Box_vec2=Box_vec2(~isnan(Box_vec2));
Box_vec3=Box_vec3(~isnan(Box_vec3));

Box_matrix=[Box_vec1;Box_vec2;Box_vec3]';

%% Box vectors for the .gro format is (free format, space separated reals), values:
% v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
% the last 6 values may be omitted (they will be set to zero) when all angles are 90
% GROMACS only supports boxes with v1(y)=v1(z)=v2(z)=0.

%% Box matrix
% v1(x) v2(x) v3(x)    v1(x) v2(x) v3(x)
% v1(y) v2(y) v3(y) == 0     v2(y) v3(y)
% v1(z) v2(z) v3(z)    0     0     v3(z)

Box_dim=[Box_matrix(1,1) Box_matrix(2,2) Box_matrix(3,3) 0 0 Box_matrix(1,2) 0 Box_matrix(1,3) Box_matrix(2,3)]
Cell=Box_dim2Cell(Box_dim)

Elements=[];
for i=1:numElements
    Elements{i} = {fgets(fileID)};
    if i==1
        ElementTypes=strsplit(char(strtrim(Elements{i})));
    else
        ElementNumbers=str2double(strsplit(char(strtrim(Elements{i}))));
    end

end

%% Close the text file.
fclose(fileID);

filetempID = fopen(filename,'r');
for i=1:startRow-1
    temp = {fgets(filetempID)};
end

for i=1:sum(ElementNumbers)
    line = fgetl(filetempID);
    XYZ_string=str2double(strsplit((line)));
    XYZ_string=XYZ_string(~isnan(XYZ_string));
    X(i) = XYZ_string(1);
    Y(i) = XYZ_string(2);
    Z(i) = XYZ_string(3);
end
fclose(filetempID);

XYZ_data=[X' Y' Z'];

XYZ_labels=[];
for i=1:numel(ElementTypes)
    XYZ_labels=[XYZ_labels; repmat(ElementTypes(i),ElementNumbers(i),1)];
end

for i=1:size(XYZ_data,1)
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
end

atom=element_atom(atom);

assignin('caller','Box_dim',Box_dim)
assignin('caller','XYZ_labels',XYZ_labels)
assignin('caller','XYZ_data',XYZ_data)

