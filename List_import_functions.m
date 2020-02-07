%% List of import functions
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Import of data files
% # <import_xvg.html import_xvg(filename)> % This function imports a Gromacs .xvg file

%% Import of structure files
% # <import_atom_car.html import_atom_car(filename,varargin)> % This function imports .car files from Hendrik Heinz INTERFACE ff distribution, and then tries to write out a Gromacs molecular topology file (.itp) and a new .pdb file
% # <import_atom_gro.html import_atom_gro(filename)> % This function imports .gro files into the atom struct
% # <import_atom_multiple_gro.html import_atom_multiple_gro(filename,nFrames)> % This function import multiple .gro files
% # <import_atom_pdb.html import_atom_pdb(filename)> % This function imports .pdb files into the atom struct
% # <import_atom_xyz.html import_atom_xyz(filename)> % This imports an .xyz file into the atom struct
% # <import_atom.html import_atom(filename)> % import_atom.m - This imports a .xyz/.gro/.pdb file and puts the data in a structure variable called atom
% # <import_xyz.html import_xyz(filename)> % This function imports an .xyz file. Atom types should be made of letters, not numbers... Try the import_atom_xyz function instead...

%% Import of trajectory files (these functions may depend other third party functions...)
% # <import_gro_traj.html import_gro_traj(filename,varargin)> % This function imports an strcture and an .gro trajectory file
% # <import_atom_multiple_gro.html import_atom_multiple_gro(filename,nFrames)> % This function import multiple .gro files
% # <import_mc_pdb_traj.html import_mc_pdb_traj(filename,varargin)> % This function imports an structure and an .pdb trajectory file, and can handle changing number of particles
% # <import_pdb_traj.html import_pdb_traj(filename,varargin)> % This function imports an strcture and an .pdb trajectory file.
% # <import_traj.html import_traj(filenameconf,filenametraj)> % This function imports an strcture and an dcd, trr, xtc, xyz or gro  trajectory file.
% # <import_trr.html import_trr(filenameconf,filenametraj)> % This function imports an structure and an trr  trajectory file
% # <import_trrv2.html import_trrv2(filenameconf,filenametraj)> % This function imports an structure and an trr  trajectory file
% # <import_xtc.html import_xtc(filenameconf,filenamextc)> % import_atom_xtc.m - This imports a structure file and a xtc file
% # <import_xtcv2.html import_xtcv2(filenameconf,filenamextc)> % import_atom_xtc.m - This imports a structure file and a xtc file
% # <import_xyz_traj.html import_xyz_traj(filenametraj)> % import_xyz_traj.m - This imports an strcture and an .xyz trajectory file.

