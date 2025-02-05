%% import_atom_gro.m
% * This function import .gro files into an atom struct variable
% * varargin can be used to translate, alt. center+translate the molecule
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_atom_gro('molecule.gro')
% # atom = import_atom_gro('molecule.gro',[10 5 2])
% # atom = import_atom_gro('molecule.gro',[10 5 0],[35.24 24.23 52.23])

function [atom,Box_dim] = import_atom_gro_legacy(filename,varargin)
tic

if exist('matlabroot', 'builtin') == 5

    if regexp(filename,'.gro') ~= false
        filename = filename;
    else
        filename = strcat(filename,'.gro');
    end

    % Get the number of atoms and Box_dim
    fileID = fopen(filename,'r');
    Line1 = {fgets(fileID)};
    Line2 = {fgets(fileID)};
    Title=strsplit(char(Line1));
    nAtoms=str2double(Line2);
    Box_string = textscan(fileID, '%s',1,'delimiter', '\n','HeaderLines', nAtoms);
    Box_dim=str2double(strsplit(char(Box_string{1,1})))*10;
    fclose(fileID);

    %% Box vectors for the .gro format is (free format, space separated reals), values:
    % v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
    % the last 6 values may be omitted (they will be set to zero) when all angles are 90
    % GROMACS only supports boxes with v1(y)=v1(z)=v2(z)=0.

    %% Box matrix
    % v1(x) v2(x) v3(x)    v1(x) v2(x) v3(x)
    % v1(y) v2(y) v3(y) == 0     v2(y) v3(y)
    % v1(z) v2(z) v3(z)    0     0     v3(z)


    % Read columns of data as strings:
    formatSpec = '%5s%5s%5s%5.0f%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f%[^\n\r]';

    % Open the text file.
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, nAtoms, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines', 2, 'ReturnOnError', false);
    fclose(fileID);

    % AtomID ResName
    nAtoms=size(dataArray{:,5}(:),1);
    MolID = str2double((dataArray{:,1})); % Converts to double
    Resname = strtrim(dataArray{:,2});
    Atomtype = strtrim(dataArray{:,3});
    ind=find(dataArray{:,4}(:)>99999);
    dataArray{1,4}(ind)=dataArray{1,4}(ind)-100000;
    AtomID = dataArray{:,4}; % Converts to double
    X_coord = dataArray{:,5}*10;
    Y_coord = dataArray{:,6}*10;
    Z_coord = dataArray{:,7}*10;
    X_velo = dataArray{:,8}*10;
    Y_velo = dataArray{:,9}*10;
    Z_velo = dataArray{:,10}*10;

    % Preallocate the atom struct attributes/fields
    atom(1).molid=[];
    atom(1).resname={};
    atom(1).type={};
    atom(1).fftype={};
    atom(1).index=[];
    atom(1).neigh.type  = {};
    atom(1).neigh.index  = [0;0;0;0;0;0];
    atom(1).neigh.dist  = [0;0;0;0;0;0];
    atom(1).bond.type  = [0;0;0;0;0;0];
    atom(1).bond.index  = [0;0;0;0;0;0];
    atom(1).angle.type  = [0;0;0;0;0;0];
    atom(1).angle.index  = [0;0;0;0;0;0];
    atom(1).x=[];
    atom(1).y=[];
    atom(1).z=[];
    atom(1).vx=[];
    atom(1).vy=[];
    atom(1).vz=[];
    atom=repmat(atom,1,nAtoms);

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
        atom(i).resname=Resname(i);
        atom(i).type=Atomtype(i);
        atom(i).fftype=Atomtype(i);
        atom(i).index=mod(i,100000);
        atom(i).neigh.type  = {};
        atom(i).neigh.index  = [0;0;0;0;0;0];
        atom(i).neigh.dist  = [0;0;0;0;0;0];
        atom(i).bond.type  = [0;0;0;0;0;0];
        atom(i).bond.index  = [0;0;0;0;0;0];
        atom(i).angle.type  = [0;0;0;0;0;0];
        atom(i).angle.index  = [0;0;0;0;0;0];
        atom(i).x=X_coord(i);
        atom(i).y=Y_coord(i);
        atom(i).z=Z_coord(i);
        %     atom(i).fx=X_coord(i)/Box_dim(1);
        %     atom(i).fy=Y_coord(i)/Box_dim(2);
        %     atom(i).fz=Z_coord(i)/Box_dim(3);
        atom(i).vx=X_velo(i);
        atom(i).vy=Y_velo(i);
        atom(i).vz=Z_velo(i);
    end
    last_in(atom(end).molid,1)=nAtoms;

    if nargin==2
        atom = translate_atom(atom,cell2mat(varargin(1))+[0 0 -median([atom.z])],'all');
    end

    if nargin==3
        atom = center_atom(atom,cell2mat(varargin(2)),'all','xyz');
        atom = translate_atom(atom,cell2mat(varargin(1))+[0 0 -median([atom.z])],'all');
    end

    XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
    XYZ_labels=[atom.type]';

    Cell=Box_dim2Cell(Box_dim);

else
    atom = import_atom_gro_octave(filename,varargin);
end

% atom = resname_atom(atom);
assignin('caller','XYZ_labels',XYZ_labels)
assignin('caller','XYZ_data',XYZ_data)
assignin('caller','atom',atom)
assignin('caller','nAtoms',nAtoms)
assignin('caller','Box_dim',Box_dim)
assignin('caller','Cell',Cell)
assignin('caller','MolID',MolID)

disp('.gro file imported')
toc

