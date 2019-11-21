%% import_atom_pdb.m
% * This function import .pdb files into the atom struct
% * varargin can be used to translate, alt. center+translate the molecule
%
%% Version
% 2.05
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_atom('molecule.pdb')
% # atom = import_atom('molecule.pdb',[10 5 2])
% # atom = import_atom('molecule.pdb',[10 5 0],[35.24 24.23 52.23])
%
function atom = import_atom_pdb(filename,varargin)

% See http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
% COLUMNS        DATA  TYPE    FIELD        DEFINITION
% -------------------------------------------------------------------------------------
% 1 -  6         Record name   "ATOM  "
% 7 - 11         Integer       Serial       Atom  serial number.
% 13 - 16        Atom          Atom type    Atom name.   ->17 by MH
% 17             Character     AltLoc       Alternate location indicator.
% 18 - 20        Residue name  ResName      Residue name.
% 22             Character     ChainID      Chain identifier.
% 23 - 26        Integer       ResSeq       Residue sequence number.
% 27             AChar         Code         Code for insertion of residues.
% 31 - 38        Real(8.3)     X            Orthogonal coordinates for X in Angstroms.
% 39 - 46        Real(8.3)     Y            Orthogonal coordinates for Y in Angstroms.
% 47 - 54        Real(8.3)     Z            Orthogonal coordinates for Z in Angstroms.
% 55 - 60        Real(6.2)     Occupancy    Occupancy.
% 61 - 66        Real(6.2)     TempFactor   Temperature  factor.
% 73 - 76        LString(4)    Segment identifier, left-justified. % Not used
% 77 - 78        LString(2)    Element      Element symbol, right-justified.
% 79 - 80        LString(2)    Charge       Charge on the atom.
filename
fid = fopen(filename,'r');
% fullText = fread(fid,'char=>char')';
% data = strread(fullText,'%s','delimiter','\n');% use textscan instead?
data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', ''); % New addition
data=data{1}; % New addition
fclose(fid);

IndexCRYS = strfind(data,'CRYS');
Index = find(not(cellfun('isempty',IndexCRYS)));

Box_cell=strsplit(char(data(Index)));Box_dim=[];
for i=1:length(Box_cell)
    [num, status] = str2num(char(Box_cell(i)));
    j=1;
    if status==1
        Box_dim=[Box_dim num];
        j=j+1;
    end
end

if length(Box_dim)>3
    Box_dim=Box_dim(1:6);
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
    if sum(find(Box_dim(4:end)))<0.0001
        Box_dim=Box_dim(1:3);
    end
end

j = 0;atom=[];
for i = 1:length(data)
    line = data{i};
    if ((length(line)>=4) && (strcmp(line(1:4), 'ATOM') || strcmp(line(1:6), 'HETATM')))
        j = j + 1;
        atom(j).molid = str2double(line(23:26));
        atom(j).resname = {strtrim(line(18:20))};
%         atom(j).type = {strtrim(line(13:16))}; % Changed to 17 for better
%         compatiblity
%         atom(j).fftype = {strtrim(line(13:16))}; % Changed to 17 for better
%         compatiblity
        atom(j).type = {strtrim(line(13:17))};
        atom(j).fftype = {strtrim(line(13:17))};
        atom(j).index = str2double(line(7:11));
        atom(j).neigh.type = {};
        atom(j).neigh.index = [0;0;0;0;0;0];
        atom(j).neigh.dist = [0;0;0;0;0;0];
        atom(j).bond.type = [0;0;0;0;0;0];
        atom(j).bond.index = [0;0;0;0;0;0];
        atom(j).angle.type = [0;0;0;0;0;0];
        atom(j).angle.index = [0;0;0;0;0;0];
        atom(j).x = str2double(line(31:38));
        atom(j).y = str2double(line(39:46));
        atom(j).z = str2double(line(47:54));
        atom(j).vx = NaN;
        atom(j).vy = NaN;
        atom(j).vz = NaN;
        
        occupancy(j,1)=str2double(line(55:60));
        tempfactor(j,1)=str2double(line(61:66));
        
        atom(j).occupancy=occupancy(j,1);
        atom(j).B=tempfactor(j,1);
    end
end

nAtoms=size(atom,2);

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

assignin('caller','occupancy',occupancy)
assignin('caller','tempfactor',tempfactor)

assignin('caller','XYZ_labels',XYZ_labels)
assignin('caller','XYZ_data',XYZ_data)
assignin('caller','atom',atom)
assignin('caller','nAtoms',nAtoms)
assignin('caller','Box_dim',Box_dim)
assignin('caller','MolID',[atom.molid])

disp('.pdb file imported')
