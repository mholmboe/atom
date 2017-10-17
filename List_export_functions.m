%% List of export functions

%% Export of structure/toplogy files
write_atom_gro(atom,Box_dim,filename_out) %  write_atom_gro.m - This writes a gro file. Does it also write velocities?
write_atom_itp(atom,Box_dim,filename,varargin) % This script creates and prints a gromacs .itp file. Works best for clayff or interface ff with spc, spce or tip3p
write_atom_lmp(atom,Box_dim,filename,varargin) % This script creates and prints a lammps data file (.lj). Works best for Clayff systems
write_atom_mol2(atom,Bond_index,Box_dim,filename_out) % This function writes a .mol2 file from the atom struct
write_atom_oplsaa_go_itp(atom,Box_dim,filename,varargin) % This custom made script creates and prints a gromacs .itp file for 
write_atom_pdb(atom,Box_dim,filename_out) % This function writes an .pdb file from the atom struct using Gromacs
write_atom_pqr(atom,Box_dim,filename_out,varargin) % This function writes an .pqr file from the atom struct
write_atom_psf(atom,Box_dim,filename_out,varargin) % This function writes an .psf file from the atom struct
write_atom_xyz(atom,Box_dim,filename_out) % This function writes an XYZ file from the atom struct
write_atom(atom,Box_dim,filename_out,varargin) % This function tries to write various files for you. Works best for systems designed for Clayff...

%% Import of trajectory files
write_atom_multiple_gro(atom,traj,filename_out) % This function writes a .gro trajectory
write_gro_traj(atom,traj,Box_dim,filename_out) % This function writes a .gro trajectory 
write_pdb_traj(atom,traj,Box_dim,filename_out) % This function writes a .pdb trajectory 
write_xyz_traj(atom,traj,Box_dim,filename_out) % This function writes a .xyz trajectory 