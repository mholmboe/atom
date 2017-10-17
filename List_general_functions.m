%% List of general functions

%% Convert water functions
spc2tip4p(filename) % This function converts a .gro or .pdb file with spc water to some tip4p water
spc2tip5p(filename) % This function converts a .gro or .pdb file with spc water to some tip5p water
spce2tip4p(filename) % This function converts a .gro or .pdb file with spc water to some tip4p water
tip3p2tip4p(filename) % This function converts a .gro file with spc water to some tip4p water
% 864_spc.gro % equilibrated SPC water box
% 864_spce.gro % equilibrated SPC/E water box
% 864_tip3p.gro % equilibrated TIP3P water box
% 864_tip4p.gro % equilibrated TIP4P water box
% 864_tip5p.gro % equilibrated TIP5P water box
% 96spc_hex_ice_h.gro % equilibrated SPC hex-ice water box
% 96tip4p_hex_ice_h.gro % equilibrated TIP4P hex-ice water box

%% Various distance and bond functions
bond_angle_atom(atom,Box_dim,max_short_dist,max_long_dist,varargin) % This function tries to find all bonds and angles of the atom struct 'more' is an optional varargin argument
bond_angle_type(atom1,atom2,Box_dim,rmin,rmax,angle_limit,varargin) % This tries to find all bonds and angles of the atom types
CONECT_atom(atom,Box_dim,short_r,long_r) % This prints conect records for pdb files
dist_matrix_atom(atom,Box_dim) % This calculates the distance matrix from the atom struct
neigh_atom(atom,Box_dim,rmax,varargin) % This function checks which neighbors each atom has and ouputs their info
new_neigh_atom(atom,Box_dim,rmax,varargin) % Not finished yet...
radius_atom(atom,ffname,watermodel) % This function fetches the ion radius from clayff or interface or interface2015 ff's and
radius_ion(Atom_label) % This function fetches the ionic radius, originally taken from the link below

%% Other general functions
atom2make_ndx(filename,groupname,atomtypes,molid) % This little script can be used to write custom made Gromacs index files
ave_atom(atom) % This function calculates the mean of the atom coordinates
cat_atom.m % This is a special script (and not a function) that imports and appends atom structs into a .gro trajectory file, useful to make a trajectory with varying number of particles
COM_atom(atom,MolID) % This function calculates the COM for certain elements
COM_func(MolID,XYZ_data,Atom_label,XYZ_labels,Box_dim) % This calculates the center of mass for water. Slow due to pbc...
COM_molid(atom,MolID) % This function calculates the COM for certain elements
COM_SOL(MolID,XYZ_data,Atom_label,XYZ_labels,Box_dim) % Computes the COM of SPC water?
composition_atom(atom) % This function looks at the composition of the atom struct
concatenate_atom(atom_1,atom_2) % This function concatenats atom sections.
condense_atom(atom,Box_dim,s) % This function tries to minimize the box size and remove gaps between molids?
dipoles_atom(Elements,Box_dim) % This function calculates the dipole vector of water. Similar to the COM_func
draw_box_atom(Box_dim,LineColor,LineThickness) % Draws a box
element_atom(atom,varargin) % Converts atomtypes to element types. This function replaces the atomtypes names with the element names
element_color(Atom_label) % This function assigns a certain color to each element. Estethic improvements are welcome...
find_bonded_atom(atom,bond_matrix,label1,label2) % This function does a cross check of the bond matrix
frac2atom(atom,Box_dim,angleparam,angletype) % This function transforms fractional coordinates to cartesian
frame2atom(atom,traj,frame,Box_dim,varargin) % This function extracts a frame to the trajectory matrix
gmx_make_ndx(groupname,ind) % This function helps you print custom gromacs .ndx files
hist_atom(atom,s) % This function is used to calculate density profiles in the Z-direction
mass_atom(atom) % This function fetches the mass for each atomtype and put it into atom.mass
median_atom(atom) % This function calculates the median position of the atom struct
orto_atom(atom,Box_dim) % This transforms a triclinic atom struct to an orthogonal atom struct. Box_dim must look like [lx ly lz 0 0 xy 0 xz yz]
PATH2VMD() % The VMD path on your computer
plot_atom(atom,Box_dim,varargin) % This function draws the atom struct in 3D. Its very simplistic with no cool features
rename_type(atom,atomtype,new_atomtype,varargin) % This function renames atoms in the atom
reorder_atom_gro(atom,atomlist,Box_dim,filename_out) % This function reorders the atoms in a .gro file
resname_atom(atom) % This function tries to guess the resname of all atom types
sort_atom(atom) % sort_atom.m - This section orders to atoms with respect to z
sort_molid(Molid) % This function sorts the molecular indexes in an ascending order
update_atom(atom) % This function updates the molid index and the atoms index in the atom struct
vmd(atom,Box_dim) % This function plots the atom struct

%% Keep/remove functions
keep_atom(atom,resname) % keep_atom.m - This removes all but resname
keep_resname(atom,resnames) % keep_resname.m - This removes all but the resnames
remove_molid(atom,MolID) %  remove_molid.m - This removes residue with molid MolID = [1 2 3 .....]
remove_residues(atom,resnames,lo,hi,dim) % This function section is used to remove residues in the simulation box between limits lo and hi
remove_resname(atom,resnames) % This function removes residue with molid MolID, resnames = {'SOL' 'Protein'}
remove_SOL(atom,atomname,lo,hi,dim) %  This section is used to remove residues in the simulation box between limits lo and hi
remove_type(atom,typescell) % This function removes atomtypes with types as in typescell = {'OW' 'HW1' 'HW2'}


