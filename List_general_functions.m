%% List of general functions

%% Convert water functions
% # <spc2tip4p.html spc2tip4p(filename)> % This function converts a .gro or .pdb file with spc water to some tip4p water
% # <spc2tip5p.html spc2tip5p(filename)> % This function converts a .gro or .pdb file with spc water to some tip5p water
% # <spce2tip4p.html spce2tip4p(filename)> % This function converts a .gro or .pdb file with spc water to some tip4p water
% # <tip3p2tip4p.html tip3p2tip4p(filename)> % This function converts a .gro file with spc water to some tip4p water

%% Various distance and bond functions
% # <bond_angle_atom.html bond_angle_atom(atom,Box_dim,max_short_dist,max_long_dist,varargin)> % This function tries to find all bonds and angles of the atom struct 'more' is an optional varargin argument
% # <bond_angle_type.html bond_angle_type(atom1,atom2,Box_dim,rmin,rmax,angle_limit,varargin)> % This tries to find all bonds and angles of the atom types
% # <bond_valence_atom.html bond_valence_atom(atom,Box_dim,varargin)> % This function tries to calculate the bond valence values according to the bond valence method
% # <bond_valence_data.html bond_valence_data(ion1,ion2,R,varargin)> % This function fetches the data and matches it to the passed atom types used to calculate the bond valence values according to http://www.iucr.org/resources/data/datasets/bond-valence-parameters
% # <cell_list_dist_matrix_atom.html cell_list_dist_matrix_atom(atom,Box_dim,varargin)> % This function calculates the distance matrix from the atom struct, using a cell list algorithm adapted from the Matlab MDtoolbox by Yasuhiro Matsunaga
% # <dist_matrix_atom.html dist_matrix_atom(atom,Box_dim)> % This calculates the distance matrix from the atom struct
% # <neigh_atom.html neigh_atom(atom,Box_dim,rmax,varargin)> % This function checks which neighbors each atom has and ouputs their info
% # <neighbor_func.html neighbor_func(solute_index,XYZ_solute,XYZ_data,Box_dim,radius)> % This function scans xyz data and checks who is within a certain radius. It outputs neighbour index,
% # <radius_atom.html radius_atom(atom,ffname,watermodel)> % This function fetches the ion radius from clayff or interface or interface2015 ff's and
% # <radius_ion.html radius_ion(Atom_label)> % This function fetches the ionic radius, originally taken from the link below
% # <radius_vdw.html radius_vdw(Atom_label)> % This function fetches the rdw radius, originally taken from below from 'A cartography of the van der Waals territories' Santiago Alvarez doi:10.1039/c3dt50599e
% # <xrd_atom.html xrd_atom(varargin)> % This function calculates theoretical XRD patterns from a .pdb|.gro file or from an atom struct and Box_dim.
% # <xrd_sc_atom.html xrd_atom(varargin)> % This function calculates theoretical XRD patterns from a .pdb|.gro file or from an atom struct and Box_dim, trying ti use the supercell approach by Ufer implemented in the BGMN XRD code.

%% Other general functions
% # <add2atom.html add2atom(XYZ_labels,XYZ_data,varargin)> % This function appends so-called XYZ atomtype labels and XYZ data to an existing atom struct
% # <analyze_atom.html analyze_atom(atom,Box_dim,max_H_dist,max_dist)> % This function analyzes variofus things of the MMT atom struct
% # <atomic_scattering_factors.html atomic_scattering_factors(Atom_label,lambda,twotheta,DW)> % This function retrieves the atomic scattering factor vs 2theta using the 11 coeff parameters from Waasmaier Kirfel, 1995
% # <ave_atom.html ave_atom(atom)> % This function calculates the mean of the atom coordinates
% # <cat_atom.html cat_atom.m> % This is a special script (and not a function) that imports and appends atom structs into a .gro trajectory file, useful to make a trajectory with varying number of particles
% # <COM_atom.html COM_atom(atom,MolID)> % This function calculates the COM for certain elements
% # <COM_func.html COM_func(MolID,XYZ_data,Atom_label,XYZ_labels,Box_dim)> % This calculates the center of mass for water. Slow due to pbc...
% # <COM_molid.html COM_molid(atom,MolID)> % This function calculates the COM for certain elements
% # <COM_SOL.html COM_SOL(MolID,XYZ_data,Atom_label,XYZ_labels,Box_dim)> % Computes the COM of SPC water?
% # <composition_atom.html composition_atom(atom)> % This function looks at the composition of the atom struct
% # <dipoles_atom.html dipoles_atom(Elements,Box_dim)> % This function calculates the dipole vector of water. Similar to the COM_func
% # <draw_box_atom.html draw_box_atom(Box_dim,LineColor,LineThickness)> % Draws a box
% # <element_atom.html element_atom(atom,varargin)>  % Converts atomtypes to element types. This function replaces the atomtypes names with the element names
% # <element_color.html element_color(Atom_label)> % This function assigns a certain color to each element. Estethic improvements are welcome...
% # <frac2atom.html frac2atom(atom,Box_dim,angleparam,angletype)> % This function transforms fractional coordinates to cartesian
% # <frame2atom.html frame2atom(atom,traj,frame,Box_dim,varargin)> % This function extracts a frame to the trajectory matrix
% # <G2_atom.html G2_atom(atom,Box_dim)> % This function calculates the continuous G2 factor fromthe cos and sin terms and also saves a struct variable for G2_calc_func(). You might wnat to edit the atomtype names below to fit your needs...
% # <hist_atom.html hist_atom(atom,s)> % This function is used to calculate density profiles in the X|Y|Z-direction
% # <histz_atom.html histz_atom(atom,s)> % This function is used to calculate density profiles in the Z-direction
% # <median_atom.html median_atom(atom)> % This function calculates the median position of the atom struct
% # <orto_atom.html orto_atom(atom,Box_dim)> % This transforms a triclinic atom struct to an orthogonal atom struct. Box_dim must look like [lx ly lz 0 0 xy 0 xz yz]
% # <PATH2GMX.html PATH2<gmx()> % The Gromacs path on your computer
% # <PATH2VMD.html PATH2VMD()> % The VMD path on your computer
% # <place_atom.html place_atom(atom,position)> % This function places the atom struct according to the position vector called position, trying to use the COM of the molecule
% # <plot_density_atom.html plot_density_atom(atom,Box_dim,varargin)> % This function draws the atom struct in 3D adjoined by some density profiles
% # <plot_atom.html plot_atom(atom,Box_dim,varargin)> % This function draws the atom struct in 3D. Its very simplistic with no cool features
% # <properties_atom.html properties_atom(atom,Box_dim,varargin)> % This function fetches various preperties of the atoms in the atom struct, using for instance the bond valence method and for instance the radii originally taken from below Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides. R. D. Shannon Acta Cryst. (1976) A32, 751-767.
% # <reduced_mass.html reduced_mass(Atom_label1,varargin)> % This function calculates the reduced mass.
% # <triclinic_atom.html triclinic_atom(atom,Box_dim,angleparam,angletype)> %  triclinic_atom.m - This transforms an orthogonal atom struct to a triclinic with the angles alfa, beta, gamma or tilt factors xy, xz, yz
% # <update_atom.html update_atom(atom)> % This function updates the molid index and the atoms index in the atom struct
% # <vmd.html vmd(atom,Box_dim)> % This function plots the atom struct



%% Keep/remove functions
% # <keep_atom.html keep_atom(atom,resname)> % keep_atom.m - This removes all but resname
% # <keep_resname.html keep_resname(atom,resnames)> % keep_resname.m - This removes all but the resnames
% # <remove_molid.html remove_molid(atom,MolID)> %  remove_molid.m - This removes residue with molid MolID = [1 2 3 .....]
% # <remove_occypancy_atom.html remove_occypancy_atom(atom)> % This function removes all succeding particles in the atom struct that has identical coordinates to a preceding particle
% # <remove_residues.html remove_residues(atom,resnames,lo,hi,dim)> % This function section is used to remove residues in the simulation box between limits lo and hi
% # <remove_resname.html remove_resname(atom,resnames)> % This function removes residue with molid MolID, resnames = {'SOL' 'Protein'}
% # <remove_SOL.html remove_SOL(atom,atomname,lo,hi,dim)> %  This section is used to remove residues in the simulation box between limits lo and hi
% # <remove_type.html remove_type(atom,typescell)> % This function removes atomtypes with types as in typescell = {'OW' 'HW1' 'HW2'}


