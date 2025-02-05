%% import_atom_xyz.m
% * This function imports an .xyz file into the atom struct
% * It tries to guess the Box_dim, so watch out!
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_atom_xyz('molecule.xyz')
%
function [atom,Box_dim] = import_atom_xyz(filename)

if regexp(filename,'.xyz') ~= false
    filename = filename;
else
    filename = strcat(filename,'.xyz');
end

fileID = fopen(filename,'r');
line1 = {fgets(fileID)};
line2 = {fgets(fileID)};
nAtoms=str2double(line1);
Box_string=strsplit(char(line2));

% startRow = 3; endRow = 2+nAtoms;
% formatSpec = '%s%f%f%f%[^\n\r]';
% delimiter = {'\t',' '};
% dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', 0,'ReturnOnError', false);
%
% for block=2:length(startRow)
%     frewind(fileID);
%     textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
%     dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
%     for col=1:length(dataArray)
%         dataArray{col} = [dataArray{col};dataArrayBlock{col}];
%     end
% end

%% Close the text file.
fclose(fileID);

XYZ_labels = cell(nAtoms,1);   % pre-allocate as a cell array
filetempID = fopen(filename,'r');
line1 = {fgets(filetempID)};
line2 = {fgets(filetempID)};
for i=1:nAtoms
    line = fgetl(filetempID);
    XYZ_string = strsplit(line);
    
    % Use curly braces to access the contents of the first token
    XYZ_labels{i,1} = XYZ_string{1};
    
    X(i) = XYZ_string(2);
    Y(i) = XYZ_string(3);
    Z(i) = XYZ_string(4);
end
fclose(filetempID);

XYZ_data=[str2double(X)' str2double(Y)' str2double(Z)'];

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
end

idx = strfind(char(Box_string(1)),'#');

Box_dim=[];
if numel(idx) >0
    for i=1:length(Box_string)
        [num, status] = str2num(char(Box_string(i)));
        j=1;
        if status==1
            Box_dim=[Box_dim num];
            j=j+1;
        end
    end
end

if length(Box_dim)==6
    Box_dim=Box_dim(1:6);
    if sum(Box_dim(4:6)) > 100 % If anglees... and not tilt-factors
        a=Box_dim(1);
        b=Box_dim(2);
        c=Box_dim(3);
        alfa=Box_dim(4);
        beta=Box_dim(5);
        gamma=Box_dim(6);
        lx = a;
        xy = b * cos(deg2rad(gamma));
        ly = (b^2-xy^2)^.5;
        xz = c*cos(deg2rad(beta));
        yz = (b*c*cos(deg2rad(alfa))-xy*xz)/ly;
        lz = (c^2 - xz^2 - yz^2)^0.5;
        Box_dim=[lx ly lz 0 0 xy 0 xz yz];
        Box_dim(Box_dim<0.00001&Box_dim>-0.00001)=0;
    elseif sum(Box_dim(4:6)) < 100
        Box_dim=[Box_dim(1:3) 0 0 Box_dim(4) 0 Box_dim(5) Box_dim(6)];
    end
    
    if sum(find(Box_dim(4:end)))<0.0001
        Box_dim=Box_dim(1:3);
    end
end

% if length(Box_dim)==3
%     Box_dim=Box_dim;
% elseif length(Box_dim)==6
%     Box_dim=[Box_dim(1:3) 0 0 Box_dim(4) 0 Box_dim(5) Box_dim(6)];
% elseif length(Box_dim)==9
%     Box_dim=Box_dim(1:9);
% end

if numel(Box_dim)==0
    disp('Guessing the box dimensions to be .1% larger than max coordinates')
    Box_dim = [(max([atom.x])-min([atom.x]))*1.001    (max([atom.y])-min([atom.y]))*1.001  (max([atom.z])-min([atom.z]))*1.001 0 0 0 0 0 0];
%     pause(2)
end

% atom = resname_atom(atom);

assignin('caller','Box_dim',Box_dim)
assignin('caller','XYZ_labels',XYZ_labels)
assignin('caller','XYZ_data',XYZ_data)

