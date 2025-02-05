%% List of import functions
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%

%% Import of topology files
% # <import_itp.html import_itp(filename)> % Import a Gromacs .itp topology file

%% Import of data files
% # <import_xvg.html import_xvg(filename)> % Import a Gromacs .xvg file.
% # <import_gmx_energy.html import_gmx_energy(filename)> % Import a Gromacs energy file.
% # <import_bar.html import_bar(filename)> % Import a BAR file.
% # <import_dat.html import_dat(filename)> % Import a .dat file.
% # <import_red_charges.html import_red_charges(filename)> % Import reduced charges from a file.
% # <import_ddec_charges.html import_ddec_charges(filename)> % Import DDEC charges from a file.

%% Import of structure files
% # <import_atom_car.html import_atom_car(filename,varargin)> % Import .car files from Hendrik Heinz INTERFACE force field distribution, then write out a Gromacs molecular topology file (.itp) and a new .pdb file.
% # <import_atom_gro.html import_atom_gro(filename)> % Import .gro files into the atom struct.
% # <import_atom_gro_fscanf.html import_atom_gro_fscanf(filename)> % Import .gro files using fscanf (alternative method).
% # <import_atom_gro_octave.html import_atom_gro_octave(filename)> % Import .gro files using Octave-compatible code.
% # <import_atom_mol2.html import_atom_mol2(filename)> % Import .mol2 files into the atom struct.
% # <import_atom_pdb.html import_atom_pdb(filename)> % Import .pdb files into the atom struct.
% # <import_atom_pqr.html import_atom_pqr(filename)> % Import .pqr files into the atom struct.
% # <import_atom_xyz.html import_atom_xyz(filename)> % Import an .xyz file into the atom struct.
% # <import_atom_poscar.html import_atom_poscar(filename)> % Import a VASP POSCAR file into the atom struct.
% # <import_atom.html import_atom(filename)> % Import a .xyz, .gro, or .pdb file into a structure variable called atom.
% # <import_xyz.html import_xyz(filename)> % Import an .xyz file. Atom types should be made of letters, not numbers. Use import_atom_xyz instead.

%% Import of trajectory files
% # <import_gro_traj.html import_gro_traj(filename,varargin)> % Import a structure and a .gro trajectory file.
% # <import_xyz_traj.html import_xyz_traj(filenametraj)> % Import a structure and an .xyz trajectory file.
% # <import_mc_pdb_traj.html import_mc_pdb_traj(filename,varargin)> % Import a structure and a .pdb trajectory file, handling changing numbers of particles.
% # <import_pdb_traj.html import_pdb_traj(filename,varargin)> % Import a structure and a .pdb trajectory file.
% # <import_traj.html import_traj(filenameconf,filenametraj)> % Import a structure and a .dcd, .trr, .xtc, .xyz, or .gro trajectory file.
% # <import_trr.html import_trr(filenameconf,filenametraj)> % Import a structure and a .trr trajectory file.
% # <import_trrv2.html import_trrv2(filenameconf,filenametraj)> % Import a structure and a .trr trajectory file (version 2).
% # <import_xtc.html import_xtc(filenameconf,filenamextc)> % Import a structure and an .xtc file.
% # <import_xtcv2.html import_xtcv2(filenameconf,filenamextc)> % Import a structure and an .xtc file (version 2).

%% Miscellaneous
% # <import_ave_gro.html import_ave_gro(filename)> % Import an averaged structure from a .gro trajectory.
% # <import_mclf.html import_mclf(filename)> % Import a multi-configurational London force (MCLF) file.
% # <import_mclf_C6.html import_mclf_C6(filename)> % Import MCLF C6 dispersion parameters.
% # <import_mclf_C8.html import_mclf_C8(filename)> % Import MCLF C8 dispersion parameters.
% # <import_mclf_C10.html import_mclf_C10(filename)> % Import MCLF C10 dispersion parameters.
% # <import_mclf_dispersion.html import_mclf_dispersion(filename)> % Import dispersion parameters for MCLF.
% # <import_cp2k.html import_cp2k(filename)> % Import a CP2K output file.
% # <import_cp2k_resp.html import_cp2k_resp(filename)> % Import RESP charges from CP2K.

%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se