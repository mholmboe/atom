%% import_atom.m
% * This function imports a .xyz|.gro|.pdb and basic .mol2 files and puts
% the data in the structure variable called atom
% * varargin can be used to translate, alt. center+translate the molecule
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_atom('molecule.gro') %
% # atom = import_atom('molecule.pdb',[10 5 2]) % Translates the atom struct by a [x y z] Ångströms
% # atom = import_atom('molecule.gro',[10 5 0],[35.24 24.23 52.23]) %
% Translates the atom struct by a [x y z] Ångströms and

function atom = import_atom(filename,varargin)
%%

if iscell(filename)
    filename=char(filename);
end

Box_dim = [0 0 0];

if regexp(filename,'.gro') > 1
    disp('Found .gro file');
    if nargin==2
        trans_vec=cell2mat(varargin(1));
        atom = import_atom_gro(filename,trans_vec);
    elseif nargin==3
        trans_vec=cell2mat(varargin(1));
        NewBox_dim=cell2mat(varargin(2));
        atom = import_atom_gro(filename,trans_vec,NewBox_dim);
    else
        atom=import_atom_gro(filename);
    end
    
elseif regexp(filename,'.pdb') > 1
    disp('Found .pdb file');
    if nargin==2
        trans_vec=cell2mat(varargin(1));
        atom = import_atom_pdb(filename,trans_vec);
    elseif nargin==3
        trans_vec=cell2mat(varargin(1));
        NewBox_dim=cell2mat(varargin(2));
        atom = import_atom_pdb(filename,trans_vec,NewBox_dim);
    else
        atom = import_atom_pdb(filename);
    end
    
    assignin('caller','occupancy',occupancy)
    assignin('caller','tempfactor',tempfactor)
    
elseif regexp(filename,'.pqr') > 1
    disp('Found .pdb file');
    if nargin==2
        trans_vec=cell2mat(varargin(1));
        atom = import_atom_pqr(filename,trans_vec);
    elseif nargin==3
        trans_vec=cell2mat(varargin(1));
        NewBox_dim=cell2mat(varargin(2));
        atom = import_atom_pqr(filename,trans_vec,NewBox_dim);
    else
        atom = import_atom_pqr(filename);
    end
    
    assignin('caller','occupancy',occupancy)
    assignin('caller','tempfactor',tempfactor)
    
elseif regexp(filename,'.mol2') > 1
    disp('Found .mol2 file');
    if nargin==2
        trans_vec=cell2mat(varargin(1));
        atom = import_atom_mol2(filename,trans_vec);
    elseif nargin==3
        trans_vec=cell2mat(varargin(1));
        NewBox_dim=cell2mat(varargin(2));
        atom = import_atom_mol2(filename,trans_vec,NewBox_dim);
    else
        atom = import_atom_mol2(filename);
    end
    
elseif regexp(filename,'.xyz') > 1
    disp('Found .xyz file');
    atom=import_atom_xyz(filename);
    %Box_dim = [(max([atom.x])-min([atom.x]))*1.0001 (max([atom.y])-min([atom.y]))*1.0001 (max([atom.z])-min([atom.z]))*1.0001];
    filenamegro='temp.gro';
    write_atom_gro(atom,Box_dim,filenamegro)
    scale=1;
    if nargin==2
        trans_vec=cell2mat(varargin(1));
        atom = import_atom_gro(filenamegro,trans_vec);
    elseif nargin==3
        trans_vec=cell2mat(varargin(1));
        NewBox_dim=cell2mat(varargin(2));
        atom = import_atom_gro(filenamegro,trans_vec,NewBox_dim);
    else
        atom=import_atom_gro(filenamegro);
    end
end

%%%%%%%%%% To analyze/find bonds
if size(atom,2)<2000
    try
        disp('Looking up bonds')
        atom=bond_atom(atom,Box_dim,2.4);
        atom=update_atom(atom);
    catch
        disp('Could not find bonds')
    end
end
%%%%%%%%%% To analyze composition
disp('Composition info:')
composition_atom(atom);
disp('The box dimensions are:')
Box_dim

disp('The cell dimensions are:')
Cell=Box_dim2Cell(Box_dim)

%atom = resname_atom(atom);
if numel(Box_dim)>0
    atom = mass_atom(atom,Box_dim);
    assignin('caller','Box_volume',Box_volume);
    assignin('caller','Box_density',Box_density);
    assignin('caller','Mw',Mw);
    assignin('caller','Mw_occupancy',Mw_occupancy);
else
    atom = mass_atom(atom);
end

assignin('caller','Atom_types',Atom_types);
assignin('caller','composition',composition);
assignin('caller','XYZ_labels',XYZ_labels)
assignin('caller','XYZ_data',XYZ_data)
assignin('caller','atom',atom)
assignin('caller','nAtoms',nAtoms)
assignin('caller','MolID',MolID)
assignin('caller','Box_dim',Box_dim)
assignin('caller','Cell',Cell)


delete('./#*'); delete('./temp*');
