%% import_ave_gro.m
% * This function import .gro files into an atom struct variable
% * This function is faster than the usual import_atom_gro but is not
% * recommended.
% 
%
%% Version
% 2.10
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_ave_gro('molecule.gro')

function atom = import_ave_gro(filename,varargin)
tic

if regexp(filename,'.gro') ~= false
    filename = filename;
else
    filename = strcat(filename,'.gro');
end

% Get the number of atoms and Box_dim
fileID = fopen(filename,'r');
Line1 = {fgets(fileID)};
Line2 = {fgets(fileID)};
%Title=strsplit(char(Line1));
nAtoms=str2double(Line2);
Box_string = textscan(fileID, '%s',1,'delimiter', '\n','HeaderLines', nAtoms);
Box_dim=str2double(strsplit(char(Box_string{1,1})))*10;
fclose(fileID);

% Read columns of data as strings:
formatSpec = '%5s%5s%5s%5.0f%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, nAtoms, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines', 2, 'ReturnOnError', false);
fclose(fileID);

% AtomID ResName
nAtoms=size(dataArray{:,5}(:),1);
MolID = str2double((dataArray{:,1})); % Converts to double
ind=find(dataArray{:,4}(:)>99999);
dataArray{1,4}(ind)=dataArray{1,4}(ind)-100000;
X_coord = single(dataArray{:,5})*10;
Y_coord = single(dataArray{:,6})*10;
Z_coord = single(dataArray{:,7})*10;

nmol=1;first_in=[1];last_in=[];
for i=1:nAtoms
    if i > 1 && MolID(i) ~= MolID(i-1)
        nmol=nmol+1;
        atom(i).molid=nmol;
        first_in(atom(i).molid,1)=i;last_in(atom(i).molid-1,1)=i-1;
    elseif i > 1
        atom(i).molid=atom(i-1).molid;
    elseif i == 1
        atom(i).molid=1;
    end
    atom(i).x=X_coord(i);
    atom(i).y=Y_coord(i);
    atom(i).z=Z_coord(i);
end
last_in(atom(end).molid,1)=nAtoms;

assignin('caller','Box_dim',Box_dim)

disp('.gro file imported')
toc

