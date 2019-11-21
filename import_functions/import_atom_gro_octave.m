%% import_atom_gro.m
% * This function import .gro files into an atom struct variable
% * This function is supposedly compatible with Octave but a little bit
% * slower than import_atom_gro_matlab.m
% * varargin can be used to translate, alt. center+translate the molecule
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_atom('molecule.gro')
% # atom = import_atom('molecule.gro',[10 5 2])
% # atom = import_atom('molecule.gro',[10 5 0],[35.24 24.23 52.23])
%
function  atom = import_atom_gro(filename,varargin)

fileID = fopen(filename, 'r');

Line1 = {fgets(fileID)};
Line2 = {fgets(fileID)};
Title=strsplit(char(Line1));
nAtoms=str2double(Line2);

tic
nmol=1;first_in=[1];last_in=[];
for i=1:nAtoms
    
    line = fgetl(fileID);
    MolID(i) = str2double(line(1:5));
% This section reorders the Molid's, if they are scrambled. Use update_atom() instead.
%     MolID(i) = str2num(line(1:5));
%     if i > 1 && MolID(i) ~= MolID(i-1)
%         nmol=nmol+1;
%         atom(i).molid=nmol;
%         first_in(atom(i).molid,1)=i;last_in(atom(i).molid-1,1)=i-1;
%     elseif i > 1
%         atom(i).molid=atom(i-1).molid;
%     elseif i == 1
%         atom(i).molid=1;
%     end
    atom(i).molid   = MolID(i);
    atom(i).resname = {strtrim(line(6:10))};
    atom(i).type    = {strtrim(line(11:15))};
    atom(i).fftype    = {strtrim(line(11:15))};
    atom(i).index=mod(i,100000);
    atom(i).neigh.type   = {};
    atom(i).neigh.index  = [0;0;0;0;0;0];
    atom(i).neigh.dist   = [0;0;0;0;0;0];
    atom(i).bond.type    = [0;0;0;0;0;0];
    atom(i).bond.index   = [0;0;0;0;0;0];
    atom(i).angle.type   = [0;0;0;0;0;0];
    atom(i).angle.index  = [0;0;0;0;0;0];
    atom(i).x            = 10*str2double(line(21:28));
    atom(i).y            = 10*str2double(line(29:36));
    atom(i).z            = 10*str2double(line(37:44));
    if numel(line)>45
        atom(i).vx           = str2double(line(45:52));
        atom(i).vy           = str2double(line(53:60));
        atom(i).vz           = str2double(line(61:68));
    else
        atom(i).vx           = NaN;
        atom(i).vy           = NaN;
        atom(i).vz           = NaN;
    end
end

toc 
Box_string = fgetl(fileID);

fclose(fileID);

Box_dim=str2double(strsplit(char(Box_string)))*10;
Box_dim(isnan(Box_dim))=[];

if nargin==2
    atom = translate_atom(atom,cell2mat(varargin(1))+[0 0 -median([atom.z])],'all');
end

if nargin==3
    atom = center_atom(atom,cell2mat(varargin(2)),'all','xyz');
    atom = translate_atom(atom,cell2mat(varargin(1))+[0 0 -median([atom.z])],'all');
end

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
XYZ_labels=[atom.type]';

% atom = resname_atom(atom);

assignin('caller','XYZ_labels',XYZ_labels)
assignin('caller','XYZ_data',XYZ_data)
assignin('caller','atom',atom)
assignin('caller','nAtoms',nAtoms)
assignin('caller','Box_dim',Box_dim)
assignin('caller','MolID',MolID)

disp('.gro file imported')

