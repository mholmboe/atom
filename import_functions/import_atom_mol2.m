%% import_atom_mol2.m
% * This function imports basic .mol2 files and the @<TRIPOS>MOLECULE, 
% @<TRIPOS>CRYSIN, @<TRIPOS>ATOM, @<TRIPOS>BOND sections into the atom 
% struct. 
% * varargin can be used to translate, alt. center+translate the molecule
% if the Box_dim information exists (which is usually not the case with 
% with .mol2 files).
%
%% Version
% 2.06
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_atom('molecule.mol2')
% # atom = import_atom('molecule.mol2',[10 5 2])
% # atom = import_atom('molecule.mol2',[10 5 0],[35.24 24.23 52.23])
%
function atom = import_atom_mol2(filename,varargin)

fid = fopen(filename,'r');
data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', ''); % New addition
data=data{1}; % New addition
fclose(fid);

IndexMOLECULE = strfind(data,'@<TRIPOS>MOLECULE');
IndexMOLECULE = find(not(cellfun('isempty',IndexMOLECULE)));

IndexCRYSIN = strfind(data,'@<TRIPOS>CRYSIN');
IndexCRYSIN = find(not(cellfun('isempty',IndexCRYSIN)));

IndexATOM = strfind(data,'@<TRIPOS>ATOM');
IndexATOM = find(not(cellfun('isempty',IndexATOM)));

IndexBOND = strfind(data,'@<TRIPOS>BOND');
IndexBOND = find(not(cellfun('isempty',IndexBOND)));

IndexSUBSTRUCTURE = strfind(data,'@<TRIPOS>SUBSTRUCTURE');
IndexSUBSTRUCTURE = find(not(cellfun('isempty',IndexSUBSTRUCTURE)));

IndexHEADTAIL = strfind(data,'@<TRIPOS>HEADTAIL');
IndexHEADTAIL = find(not(cellfun('isempty',IndexHEADTAIL)));

IndexRESIDUECONNECT = strfind(data,'@<TRIPOS>RESIDUECONNECT');
IndexRESIDUECONNECT = find(not(cellfun('isempty',IndexRESIDUECONNECT)));

MOLECULE=data{IndexMOLECULE+1};
MOLECULEPROPERTIES=data{IndexMOLECULE+2};
MOLECULEPROPERTIES=strsplit(data{IndexMOLECULE+2});
if numel(MOLECULEPROPERTIES{1})==0
    MOLECULEPROPERTIES(1)=[];
end
MOLECULEPROPERTIES=str2double(MOLECULEPROPERTIES)
nAtoms=MOLECULEPROPERTIES(1);
nBonds=MOLECULEPROPERTIES(2);

if numel(IndexCRYSIN)<0
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
        if length(Box_dim)>6
            Spacegroup=Box_dim(7);
            if length(Box_dim)>7
                Setting=Box_dim(8);
            end
        end
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
else
    disp('No Box information found in the .mol2 file')
    Box_dim=[0 0 0];
end

j = 0;atom=[];
for i = IndexATOM+1:IndexATOM+nAtoms
        line = data{i};
    if length(line)>=60
        j = j + 1;
        atom(j).molid = str2double(line(47:49));
        atom(j).resname = {strtrim(line(53:55))};
        atom(j).type = {strtrim(line(5:7))};
        atom(j).fftype = {strtrim(line(44:46))};
        atom(j).index = str2double(line(1:4));
        atom(j).neigh.type = {};
        atom(j).neigh.index = [0;0;0;0;0;0];
        atom(j).neigh.dist = [0;0;0;0;0;0];
        atom(j).bond.type = [0;0;0;0;0;0];
        atom(j).bond.index = [0;0;0;0;0;0];
        atom(j).angle.type = [0;0;0;0;0;0];
        atom(j).angle.index = [0;0;0;0;0;0];
        atom(j).x = str2double(line(10:20));
        atom(j).y = str2double(line(21:31));
        atom(j).z = str2double(line(32:42));
        atom(j).vx = NaN;
        atom(j).vy = NaN;
        atom(j).vz = NaN;
        atom(j).charge=str2double(line(57:63));
    end
end

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
XYZ_labels=[atom.type]';

if nBonds>0
    disp('No checks for PBC was made when reading the bonds from the mol2 file')
j = 0;Bond_index=[];
for i = IndexBOND+1:IndexBOND+nBonds
    line = data{i};
    j = j + 1;
    Bond_index(j,1) = str2double(line(6:11));
    Bond_index(j,2) = str2double(line(12:17));
    rx=XYZ_data(Bond_index(j,2),1)-XYZ_data(Bond_index(j,1),1);
    ry=XYZ_data(Bond_index(j,2),2)-XYZ_data(Bond_index(j,1),2);
    rz=XYZ_data(Bond_index(j,2),3)-XYZ_data(Bond_index(j,1),3);
    Bond_index(j,3)=sqrt( rx(:,1).^2 + ry(:,1).^2 + rz(:,1).^2 );
end
end

if nargin==2
    atom = translate_atom(atom,cell2mat(varargin(1))+[0 0 -median([atom.z])],'all');
end

if nargin==3
    atom = center_atom(atom,cell2mat(varargin(2)),'all','xyz');
    atom = translate_atom(atom,cell2mat(varargin(1))+[0 0 -median([atom.z])],'all');
end

assignin('caller','Bond_index',Bond_index)
assignin('caller','XYZ_labels',XYZ_labels)
assignin('caller','XYZ_data',XYZ_data)
assignin('caller','atom',atom)
assignin('caller','nAtoms',nAtoms)
assignin('caller','Box_dim',Box_dim)
assignin('caller','MolID',[atom.molid])

disp('.mol2 file imported')
