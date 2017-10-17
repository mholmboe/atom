%% List of import functions

%% Import of data files
import_xvg(filename) % This function imports a Gromacs .xvg file

%% Import of structure files
import_atom_car(filename,varargin) % This function imports .car files from Hendrik Heinz INTERFACE ff distribution, and then tries to write out a Gromacs molecular topology file (.itp) and a new .pdb file
import_atom_gro(filename) % This function imports .gro files into the atom struct
import_atom_multiple_gro(filename,nFrames) % This function import multiple .gro files
import_atom_pdb(filename) % This function imports .pdb files into the atom struct
import_atom_xyz_matlab(filename) % This function imports an .xyz file into the atom struct
import_atom_xyz(filename) % This imports an .xyz file into the atom struct
import_atom(filename) % import_atom.m - This imports a .xyz/.gro/.pdb file and puts the data in a structure variable called atom
import_xyz(filename) % This function imports an .xyz file. Atom types should be made of letters, not numbers... Try the import_atom_xyz function instead...

%% Import of trajectory files
import_gro_traj(filename,varargin) % This function imports an strcture and an .gro trajectory file
import_pdb_traj(filename,varargin) % This function imports an strcture and an .pdb trajectory file.
import_traj(filenameconf,filenametraj) % This function imports an strcture and an dcd, trr, xtc, xyz or gro  trajectory file.
import_trr(filenameconf,filenametraj) % This function imports an strcture and an trr  trajectory file
import_trrv2(filenameconf,filenametraj) % This function imports an strcture and an trr  trajectory file
import_xtc(filenameconf,filenamextc) % import_atom_xtc.m - This imports a structure file and a xtc file
import_xtcv2(filenameconf,filenamextc) % import_atom_xtc.m - This imports a structure file and a xtc file
import_atom_xtc_v2(filenameconf,filenamextc) % This function imports a structure file and a xtc file
import_atom_xtc(filenameconf,filenamextc) % This function imports a structure file and a xtc file
import_xyz_traj(filenametraj) % import_xyz_traj.m - This imports an strcture and an .xyz trajectory file.
