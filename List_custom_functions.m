%% List of custom Gromacs and topology tools
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%

%% Gromacs index (.ndx) file tools
% # <atom_make_ndx.html atom_make_ndx(filename,groupname,atomtypes,molid)> % Print custom Gromacs .ndx files for the given atom types and MolID's.
% # <atom_make_ndx_bonded.html atom_make_ndx_bonded(filename,groupname,atomtypes)> % Print one custom Gromacs .ndx group based on the bonded atom types.
% # <atom_mk_angndx.html atom_mk_angndx(filename,groupname,atomtypes)> % Print one custom Gromacs angle .ndx group based on the atom types.
% # <atom_read_ndx.html atom_read_ndx(filename,varargin)> % Read Gromacs index (.ndx) files into arrays.
% # <gmx_make_ndx.html gmx_make_ndx(id,groupname,varargin)> % Print custom Gromacs .ndx files from a list of atom indexes.
% # <gmx_mk_angndx.html gmx_mk_angndx(atom,Box_dim,atomtype1,atomtype2,atomtype3,groupname,varargin)> % Print custom Gromacs angle .ndx files with triplet indices for atomtype1-3.
% # <import_ndx.html import_ndx(filename,varargin)> % Read Gromacs index (.ndx) files into arrays.
% # <pdb2make_ndx.html pdb2make_ndx(filename,groupname,atomtypes,molid)> % Print custom Gromacs .ndx files directly from a .pdb file.

%% Gromacs .mdp and .itp file tools
% # <import_itp.html import_itp(filename)> % Import a Gromacs .itp topology file
% # <import_mdp.html import_mdp(mdp_filename)> % Import a Gromacs .mdp file into a mdp struct.
% # <modify_itp.html modify_itp(itp,index,varargin)> % Modify an itp struct, currently by adding values to all indexed values or removing single ones.
% # <write_itp.html write_itp(itp,filename)> % Exports a Gromacs .itp topology file from a itp struct
% # <write_mdp.html write_mdp(mdp,varargin)> % Write a Gromacs .mdp file from an imported mdp struct.

%% Trajectory tools
% # <decompose_traj.html decompose_traj(trajname,outtrajname,varargin)> % Template function used to extract trajectory components from a Gromacs .xtc or .trr file.
% # <extract_vxvy_from_trr.html extract_vxvy_from_trr(trajname,outtrajname,varargin)> % Extract the x and y components of the velocities in a Gromacs .trr file, using libxdrfile.
% # <extract_vz_from_trr.html extract_vz_from_trr(trajname,outtrajname,varargin)> % Extract the z component of the velocities in a Gromacs .trr file, using libxdrfile.

%% Miscellaneous
% # <fit2lattice.html fit2lattice> % Special script that imports a model like PO43- and tries to fit it into a crystal lattice holding such sites.

%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
