%% import_atom_car.m
% * This function imports .car files from Hendrik Heinz INTERFACE ff 
% distribution, and then tries to write out a Gromacs molecular topology 
% file (.itp) and a new .pdb file.
% * varargin could be ...,remove_atomtype,[center to this Box_dim],
% [translate_vector])
% * The remove_atomtype char/cell could be used to remove counterions like
% NA+ from the initial structures...
%
%% Version
% 2.0
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_atom_car('molecule.car')
% # atom = import_atom_car('molecule.car','no_counterions')

function atom = import_atom_car(filename,varargin)


fid = fopen(filename,'r');
fullText = fread(fid,'char=>char')';
fclose(fid); 

data = strread(fullText,'%s','delimiter','\n'); % use textscan instead?

IndexPBC = strfind(data,'PBC');
Index = find(not(cellfun('isempty',IndexPBC)));
IndexEND = strfind(data,'end');
IndexEND = find(not(cellfun('isempty',IndexEND)));
IndexEND = IndexEND(end);

Box_cell=strsplit(char(data(Index(end))));Box_dim=[];
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
    if sum(find(Box_dim(4:end)))<0.0001
        Box_dim=Box_dim(1:3);
    end
end

cardata=data(Index(end)+1:IndexEND-1,:);
j = 0;atom=[];
for i = 1:length(cardata)
    line = cardata{i};
    if numel(line) > 5
        j = j + 1;
        atom(j).molid = 1;%str2double(line(57));
        atom(j).resname = {strtrim(line(52:55))};
%         atom(j).type = {strtrim(line(13:16))}; % Changed to 17 for better
%         compatiblity
%         atom(j).fftype = {strtrim(line(13:16))}; % Changed to 17 for better
%         compatiblity
        atom(j).type = {strtrim(line(1:4))};
        atom(j).fftype = {upper(strtrim(line(64:67)))};
        atom(j).index = j;
        atom(j).neigh.type = {};
        atom(j).neigh.index = zeros(6,1);
        atom(j).neigh.dist = zeros(6,1);
        atom(j).bond.type = zeros(6,1);
        atom(j).bond.index = zeros(6,1);
        atom(j).angle.type = zeros(6,1);
        atom(j).angle.index = zeros(6,1);
        atom(j).x = str2double(line(7:20));
        atom(j).y = str2double(line(22:35));
        atom(j).z = str2double(line(37:50));
        atom(j).vx = NaN;
        atom(j).vy = NaN;
        atom(j).vz = NaN;
        atom(j).element = str2double(line(72:73));
        atom(j).charge = str2double(line(74:80));
    end
end
% Box_dim(3)=21;
% atom=replicate_atom(atom,Box_dim,[1 1 2]);

if length([atom(1).resname]) > 3
    temp=[atom(1).resname];
    [atom.resname]=deal({upper(temp(1:3))});
end

if regexp(char([atom(1).resname]),'X') ~= false
    [atom.resname]=deal({upper(filename(1:3))});
end

atom = resname_atom(atom);

if nargin>1
    remove_type=varargin(1)
    remove_type={'NA+' 'K+' 'CA2+' 'LI+'};
    ind_rm=ismember([atom.type],remove_type);
    ind_rm2=ismember([atom.fftype],remove_type);
    ind_rm=unique([find(ind_rm) find(ind_rm2)]);
    ion=atom(ind_rm);
    atom(ind_rm)=[];
%     atomwion=update_atom({atom ion});
%     write_atom_pdb(atomwion,Box_dim,strcat(filename(1:end-4),'_gmx.pdb'));
    atom=update_atom(atom);
end

nAtoms=size(atom,2);

if nargin==3
    atom = translate_atom(atom,cell2mat(varargin(2))+[0 0 -median([atom.z])],'all');
end

if nargin==4
    atom = center_atom(atom,cell2mat(varargin(3)),'all','xyz');
    atom = translate_atom(atom,cell2mat(varargin(2))+[0 0 -median([atom.z])],'all');
end

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
XYZ_labels=[atom.type]';

atom = resname_atom(atom);

assignin('caller','XYZ_labels',XYZ_labels)
assignin('caller','XYZ_data',XYZ_data)
assignin('caller','atom',atom)
assignin('caller','nAtoms',nAtoms)
assignin('caller','Box_dim',Box_dim)
assignin('caller','MolID',[atom.molid])

disp('.car file imported')
disp('and the charge was found to be...')
sum([atom.charge])
 
write_atom_psf(atom,Box_dim,strcat(filename(1:end-4)),1.25,2.25,'interface_car','tip3p')
write_atom_itp(atom,Box_dim,strcat(filename(1:end-4)),1.25,2.25,'interface_car','tip3p')
write_atom_pdb(atom,Box_dim,strcat(filename(1:end-4),'.pdb'));

