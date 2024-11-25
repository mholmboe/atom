%% List of general functions
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%

%% Convert water functions
% # <spc2tip4p.html spc2tip4p(filename)> % Convert a .gro or .pdb file with spc water to tip4p water.
% # <spc2tip5p.html spc2tip5p(filename)> % Convert a .gro or .pdb file with spc water to tip5p water.
% # <spce2tip4p.html spce2tip4p(filename)> % Convert a .gro or .pdb file with spce water to tip4p water.
% # <tip3p2tip4p.html tip3p2tip4p(filename)> % Convert a .gro file with tip3p water to tip4p water.

%% Various distance and bond functions
% # <bond_angle_atom.html bond_angle_atom(atom,Box_dim,max_short_dist,max_long_dist,varargin)> % Find all bonds and angles of the atom struct. 'More' is an optional argument.
% # <bond_angle_dihedral_atom.html bond_angle_dihedral_atom(atom,Box_dim,varargin)> % Find all bonds, angles, and dihedrals of the atom struct. Optional arguments include Rmaxshort and Rmaxlong.
% # <bond_angle_type.html bond_angle_type(atom1,atom2,Box_dim,rmin,rmax,angle_limit,varargin)> % Find all bonds and angles of the atom types.
% # <bond_valence_atom.html bond_valence_atom(atom,Box_dim,varargin)> % Calculate bond valence values using the bond valence method.
% # <bond_valence_data.html bond_valence_data(ion1,ion2,R,varargin)> % Fetch data to calculate bond valence values for specified atom types.
% # <dist_matrix_atom.html dist_matrix_atom(atom,Box_dim)> % Calculate the distance matrix from the atom struct.
% # <neigh_atom.html neigh_atom(atom,Box_dim,rmax,varargin)> % Check neighbors for each atom and output their information.
% # <neighbor_func.html neighbor_func(solute_index,XYZ_solute,XYZ_data,Box_dim,radius)> % Scan xyz data and check neighbors within a specified radius.
% # <radius_atom.html radius_atom(atom,ffname,watermodel)> % Fetch the ion radius from Clayff, Interface, or Interface2015 force fields.
% # <radius_ion.html radius_ion(Atom_label)> % Fetch the ionic radius.
% # <radius_vdw.html radius_vdw(Atom_label)> % Fetch the van der Waals radius.
% # <rdf_atom.html rdf_atom(atom,Box_dim,varargin)> % Calculate the radial distribution function and coordination number.
% # <xrd_atom.html xrd_atom(varargin)> % Calculate theoretical XRD patterns from a .pdb, .gro file, or atom struct.

%% Other general functions
% # <add2atom.html add2atom(XYZ_labels,XYZ_data,varargin)> % Append XYZ atom type labels and XYZ data to an existing atom struct.
% # <analyze_atom.html analyze_atom(atom,Box_dim)> % Fetch properties of atoms using methods like bond valence and radii.
% # <atomic_scattering_factors.html atomic_scattering_factors(Atom_label,lambda,twotheta,DW)> % Retrieve atomic scattering factors vs. 2theta.
% # <ave_atom.html ave_atom(atom)> % Calculate the mean of the atom coordinates.
% # <Box_dim2Cell.html Box_dim2Cell(Box_dim)> % Transform the Box_dim variable to the Cell variable.
% # <Cell2Box_dim.html Cell2Box_dim(Cell)> % Transform the Cell variable into the Box_dim variable.
% # <COM_atom.html COM_atom(atom,MolID)> % Calculate the center of mass (COM) for specified elements.
% # <COM_molid.html COM_molid(atom,MolID)> % Calculate the COM for specific elements.
% # <COM_SOL.html COM_SOL(MolID,XYZ_data,Atom_label,XYZ_labels,Box_dim)> % Compute the COM of SPC water.
% # <composition_atom.html composition_atom(atom)> % Analyze the composition of the atom struct.
% # <density_atom.html density_atom(atom,Box_dim)> % Calculate concentration, electron density profiles, and charge density.
% # <density_atom_fast.html density_atom_fast(atom,Box_dim)> % A faster version of the density_atom function.
% # <draw_box_atom.html draw_box_atom(Box_dim,LineColor,LineThickness)> % Draw a box to visualize dimensions.
% # <element_atom.html element_atom(atom,varargin)> % Convert atom type names to element names.
% # <element_color.html element_color(Atom_label)> % Assign a specific color to each element.
% # <frac2atom.html frac2atom(atom,Box_dim,angleparam,angletype)> % Transform fractional coordinates to Cartesian coordinates.
% # <frame2atom.html frame2atom(atom,traj,frame,Box_dim,varargin)> % Extract a frame to the trajectory matrix.
% # <G2_atom.html G2_atom(atom,Box_dim)> % Calculate the continuous G2 factor.
% # <hist_atom.html hist_atom(atom,s)> % Calculate density profiles in the X, Y, or Z direction.
% # <histz_atom.html histz_atom(atom,s)> % Calculate density profiles in the Z direction.
% # <molecule_atom.html molecule_atom(atom,varargin)> % Set molecule ID, residue name, and element names of the atom struct.
% # <median_atom.html median_atom(atom)> % Calculate the median position of the atom struct.
% # <neutralize_atom.html neutralize_atom(atom)> % Set the charge of all atom types to zero.
% # <orto_atom.html orto_atom(atom,Box_dim)> % Transform a triclinic atom struct to an orthogonal one.
% # <place_atom.html place_atom(atom,position)> % Place the atom struct at a specified position.
% # <plot_density_atom.html plot_density_atom(atom,Box_dim,varargin)> % Draw the atom struct in 3D along with density profiles.
% # <plot_atom.html plot_atom(atom,Box_dim,varargin)> % Draw the atom struct in 3D.
% # <plot_xvg.html plot_xvg(filename)> % Plot .xvg files, typically used in GROMACS simulations.
% # <properties_atom.html properties_atom(atom)> % Analyze properties of the atom struct.
% # <reduced_mass.html reduced_mass(Atom_label1,varargin)> % Calculate the reduced mass.
% # <show_atom.html show_atom(atom,varargin)> % Draw the atom struct in 3D with additional features.
% # <show_arrow.html show_arrow(p1,p2,varargin)> % Plot a 3D arrow as a patch object.
% # <show_axis.html show_axis(varargin)> % Draw the axis in a plot.
% # <show_box.html show_box(Box_dim)> % Draw the simulation box.
% # <show_Hbonds_atom.html show_Hbonds_atom(atom)> % Display or calculate hydrogen bonds in the atom struct.
% # <show_miller.html show_miller(Box_dim)> % Draw the Miller planes of the Box_dim/Cell variables.
% # <show_density_atom.html show_density_atom(atom)> % Display density of atoms.
% # <triclinic_atom.html triclinic_atom(atom,Box_dim,angleparam,angletype)> % Transform an orthogonal atom struct to a triclinic one.
% # <update_atom.html update_atom(atom)> % Update molecule and atom indices in the atom struct.
% # <vmd.html vmd(atom,Box_dim)> % Plot the atom struct using VMD.
% # <closest_atom.html closest_atom(atom,Box_dim,ref_atom)> % Find the closest atom to the reference atom.
% # <insert_atom.html insert_atom(atom,new_atom,position)> % Insert a new atom at the specified position.
% # <remove_duplicate_atom.html remove_duplicate_atom(atom)> % Remove duplicate atoms from the atom struct.
% # <sort_atom.html sort_atom(atom,Box_dim,varargin)> % Sort the atom struct based on specified properties.
% # <split_atom.html split_atom(atom,Box_dim,varargin)> % Split the atom struct into separate parts based on criteria.
% # <translate_atom.html translate_atom(atom,trans_vec)> % Translate the atom struct by a specified vector.
% # <rotate_atom.html rotate_atom(atom,rotation_matrix)> % Rotate the atom struct by a specified matrix.
% # <spiral_atom.html spiral_atom(atom,Box_dim,varargin)> % Create or modify spiral structures in the atom struct.
% # <round_atom.html round_atom(atom)> % Round atom positions or coordinates.
% # <mass_atom.html mass_atom(atom)> % Calculate or define mass for atoms.
% # <radius_crystal.html radius_crystal(Atom_label)> % Fetch or calculate ionic radii for crystal structures.
% # <sigma_vdw.html sigma_vdw(Atom_label)> % Compute sigma values for van der Waals interactions.
% # <Bragg.html Bragg(varargin)> % Calculate Bragg peaks for crystallographic or XRD data.
% # <unreplicate_atom.html unreplicate_atom(atom)> % Remove replicated atoms from the structure.

%% Keep/remove functions
% # <keep_atom.html keep_atom(atom,resname)> % Keep only specified residue names in the atom struct.
% # <keep_resname.html keep_resname(atom,resnames)> % Keep only specified residue names.
% # <remove_molid.html remove_molid(atom,MolID)> % Remove residues with specified molecule IDs.
% # <remove_occypancy_atom.html remove_occypancy_atom(atom)> % Remove particles with identical coordinates to preceding ones.
% # <remove_residues.html remove_residues(atom,resnames,lo,hi,dim)> % Remove residues within specified limits.
% # <remove_resname.html remove_resname(atom,resnames)> % Remove residues by name.
% # <remove_SOL.html remove_SOL(atom,atomname,lo,hi,dim)> % Remove solvent residues within specified limits.
% # <remove_type.html remove_type(atom,typescell)> % Remove specified atom types.

%% Lattice fitting functions
% # <fit2lattice_atom.html fit2lattice_atom(atom,Box_dim)> % Fit atoms to a lattice.
% # <fit2lattice_atom_v2.html fit2lattice_atom_v2(atom,Box

%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
