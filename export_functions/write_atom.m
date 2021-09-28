%% write_atom.m
% * This function imports a .xyz|.gro|.pdb file and puts the data in the structure variable called atom
% * Not finished for all writing functions yet... Tested 15/04/2017
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = write_atom(atom,Box_dim,'molecule.gro')
% # atom = write_atom(atom,Box_dim,'molecule.pdb')
% # atom = write_atom(atom,Box_dim,'molecule.xyz')
% # atom = write_atom(atom,Box_dim,'molecule.pqr')
% # atom = write_atom(atom,Box_dim,'molecule.mol2')
% # atom = write_atom(atom,Box_dim,'molecule.cif')
% # atom = write_atom(atom,Box_dim,'molecule.sdf')
% # atom = write_atom(atom,Box_dim,'molecule.itp')

function atom = write_atom(atom,Box_dim,filename,varargin)

if iscell(filename)
    filename=char(filename);
end

if regexp(filename,'.gro') > 1
    disp('Will try to write a .gro file');
    atom=write_atom_gro(atom,Box_dim,filename);
elseif regexp(filename,'.pdb') > 1
    disp('Will try to write a .pdb file');
    atom=write_atom_pdb(atom,Box_dim,filename);
elseif regexp(filename,'.xyz') > 1
    disp('Will try to write a .xyz file');
    atom=write_atom_xyz(atom,Box_dim,filename);
elseif regexp(filename,'.mol2') > 1
    disp('Will try to write a .mol2 file');
    atom=write_atom_mol2(atom,Box_dim,filename);
elseif regexp(filename,'.pqr') > 1
    disp('Will try to write a .pqr file');
    atom=write_atom_pqr(atom,Box_dim,filename);
elseif regexp(filename,'.cif') > 1
    disp('Will try to write a .cif file');
    atom=write_atom_cif(atom,Box_dim,filename);
elseif regexp(filename,'.sdf') > 1
    disp('Will try to write a .sdf file');
    atom=write_atom_sdf(atom,Box_dim,filename);
elseif regexp(filename,'.itp') > 1
    disp('Will try to write a .itp file');
    atom=write_atom_itp(atom,Box_dim,filename);
end

