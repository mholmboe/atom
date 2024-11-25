%% List of export functions
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%

%% Export of structure/topology files
% # <write_atom_all.html write_atom_all(atom,Box_dim,filename_out,varargin)> % Write various file types for the atom struct, best suited for Clayff systems.
% # <write_atom_cif.html write_atom_cif(atom,Box_dim,filename_out)> % Write a basic .cif file from the atom struct.
% # <write_atom_gro.html write_atom_gro(atom,Box_dim,filename_out)> % Write a .gro file from the atom struct, optionally including velocities.
% # <write_atom_itp.html write_atom_itp(atom,Box_dim,filename_out,varargin)> % Create and print a Gromacs .itp file, best for Clayff or interface force fields.
% # <write_atom_lmp.html write_atom_lmp(atom,Box_dim,filename_out,varargin)> % Create and print a LAMMPS data file (.lj), suited for Clayff systems.
% # <write_atom_mol2.html write_atom_mol2(atom,Bond_index,Box_dim,filename_out)> % Write a .mol2 file from the atom struct.
% # <write_atom_oplsaa_go_itp.html write_atom_oplsaa_go_itp(atom,Box_dim,filename_out,varargin)> % Create and print a Gromacs .itp file for OPLS-AA or GO systems.
% # <write_atom_pdb.html write_atom_pdb(atom,Box_dim,filename_out)> % Write a .pdb file from the atom struct using Gromacs.
% # <write_atom_pqr.html write_atom_pqr(atom,Box_dim,filename_out,varargin)> % Write a .pqr file from the atom struct.
% # <write_atom_psf.html write_atom_psf(atom,Box_dim,filename_out,varargin)> % Write a .psf file from the atom struct.
% # <write_atom_xyz.html write_atom_xyz(atom,Box_dim,filename_out)> % Write an XYZ file from the atom struct.
% # <write_atom.html write_atom(atom,Box_dim,filename_out,varargin)> % Write different file types (.gro, .pdb, .xyz, .itp, etc.) based on filename and parameters.
% # <write_atom_sdf.html write_atom_sdf(atom,Box_dim,filename_out)> % Write an SDF file from the atom struct.
% # <write_atom_dodecahedron_gro.html write_atom_dodecahedron_gro(atom,Box_dim,filename_out)> % Write a .gro file using a dodecahedron-shaped simulation box.
% # <write_atom_top.html write_atom_top(atom,Box_dim,filename_out)> % Write a topology file (.top) from the atom struct.
% # <CONECT_atom.html CONECT_atom(atom,Box_dim,filename_out)> % Write CONECT records for a PDB file.

%% Export of trajectory files
% # <write_gro_traj.html write_gro_traj(atom,traj,Box_dim,filename_out)> % Write a .gro trajectory file.
% # <write_atom_multiple_gro.html write_atom_multiple_gro(atom,traj,filename_out)> % Write multiple .gro files for trajectory output.
% # <write_pdb_traj.html write_pdb_traj(atom,traj,Box_dim,filename_out)> % Write a .pdb trajectory file.
% # <write_xyz_traj.html write_xyz_traj(atom,traj,Box_dim,filename_out)> % Write a .xyz trajectory file.
% # <write_ave_gro.html write_ave_gro(atom,traj,Box_dim,filename_out)> % Write an average structure from a .gro trajectory.
% # <write_ave_pdb.html write_ave_pdb(atom,traj,Box_dim,filename_out)> % Write an average structure from a .pdb trajectory.

%% Export of other formats and data
% # <export_ndx.html export_ndx(atom,Box_dim,filename_out)> % Export an index file (.ndx) from the atom struct.
% # <write_xvg.html write_xvg(filename, data)> % Export data in Gromacs .xvg format.
% # <write_tabulated_potentials.html write_tabulated_potentials(filename, data)> % Write tabulated potential files.
% # <write_ff.html write_ff(atom,filename_out)> % Write force field parameters.
% # <write_ffnonbonded.html write_ffnonbonded(atom,filename_out)> % Write non-bonded force field parameters.
% # <write_ffnonbonded_C6C12.html write_ffnonbonded_C6C12(atom,filename_out)> % Write non-bonded parameters (C6, C12) for a force field.

%% Miscellaneous
% # <replace_string.html replace_string(filename_in,filename_out,old_string,new_string)> % Replace strings in files.

%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%