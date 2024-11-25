%% Complete list of all atom functions
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%

%% All functions
% # <Atom_label_variable.html Atom_label> % 1xn cell array 
% # <Box_dim2Cell.html Box_dim2Cell(Box_dim)> % Transform the Box_dim variable to the Cell variable.
% # <Box_dim_variable.html Box_dim> % A 1x3 or 1x9 array variable holding the Box size parameters
% # <Bragg.html Bragg(varargin)> % Calculate Bragg peaks for crystallographic or XRD data.
% # <COM_SOL.html COM_SOL(MolID,XYZ_data,Atom_label,XYZ_labels,Box_dim)> % Compute the COM of SPC water.
% # <COM_atom.html COM_atom(atom,MolID)> % Calculate the center of mass (COM) for specified elements.
% # <COM_molid.html COM_molid(atom,MolID)> % Calculate the COM for specific elements.
% # <CONECT_atom.html CONECT_atom(atom,Box_dim,filename_out)> % Write CONECT records for a PDB file.
% # <Cell2Box_dim.html Cell2Box_dim(Cell)> % Transform the Cell variable into the Box_dim variable.
% # <Cell_variable.html Cell> % A 1x6 array variable containing the a, b, c cell values and  the alfa, beta, gamma angle values as used in a typical .pdb file
% # <G2_atom.html G2_atom(atom,Box_dim)> % Calculate the continuous G2 factor.
% # <XYZ_data_variable.html XYZ_data> % A nx3 matrix holdnig the XYZ-coordinates
% # <XYZ_labels_variable.html XYZ_labels> % A cell list om atom types
% # <add2atom.html add2atom(XYZ_labels,XYZ_data,varargin)> % Append XYZ atom type labels and XYZ data to an existing atom struct.
% # <adjust_H_atom.html adjust_H_atom(atom,Box_dim)> % Adjust hydrogen atoms in the atom struct.
% # <analyze_atom.html analyze_atom(atom,Box_dim)> % Fetch properties of atoms using methods like bond valence and radii.
% # <atom_variable.html atom> % The main matlab struct variable
% # <atomic_scattering_factors.html atomic_scattering_factors(Atom_label,lambda,twotheta,DW)> % Retrieve atomic scattering factors vs. 2theta.
% # <autofit_2xljcoul.html autofit_2xljcoul(atom,Box_dim,varargin)> % Automated fitting of 2x Lennard-Jones and Coulomb parameters.
% # <autofit_2xljcoul_batch.html autofit_2xljcoul_batch(atom,Box_dim,varargin)> % Batch fitting for 2x Lennard-Jones and Coulomb potentials.
% # <autofit_2xljcoul_func.html autofit_2xljcoul_func(atom,Box_dim,varargin)> % Function for automatic fitting of 2x Lennard-Jones and Coulomb parameters.
% # <autofit_C6C8C10xljcoul.html autofit_C6C8C10xljcoul(atom,Box_dim,varargin)> % Automated fitting for Lennard-Jones Coulombic parameters with C6, C8, C10 terms.
% # <autofit_buckcoul.html autofit_buckcoul(atom,Box_dim,varargin)> % Automated fitting for Buckingham and Coulomb potentials.
% # <autofit_force_2xljcoul.html autofit_force_2xljcoul(atom,Box_dim,varargin)> % Automated fitting of force parameters for 2x Lennard-Jones and Coulomb interactions.
% # <autofit_geometric2LB.html autofit_geometric2LB(atom,Box_dim,varargin)> % Automated geometric fitting for Lennard-Jones potential with Born-Mayer interactions.
% # <autofit_ljcoul.html autofit_ljcoul(atom,Box_dim,varargin)> % Automated fitting for Lennard-Jones and Coulombic interactions.
% # <ave_atom.html ave_atom(atom)> % Calculate the mean of the atom coordinates.
% # <bend_atom.html bend_atom(atom,Box_dim,Radii)> % Bend an atom struct.
% # <bond_angle_atom.html bond_angle_atom(atom,Box_dim,max_short_dist,max_long_dist,varargin)> % Find all bonds and angles of the atom struct. 'More' is an optional argument.
% # <bond_angle_atom.html bond_angle_atom(atom,Box_dim,varargin)> % This function tries to find all bonds and angles of the atom struct 'more' is an optional varargin argument.
% # <bond_angle_dihedral_atom.html bond_angle_dihedral_atom(atom,Box_dim,varargin)> % Find all bonds, angles, and dihedrals of the atom struct. Optional arguments include Rmaxshort and Rmaxlong.
% # <bond_angle_type.html bond_angle_type(atom1,atom2,Box_dim,rmin,rmax,angle_limit,varargin)> % Find all bonds and angles of the atom types.
% # <bond_atom.html bond_atom(atom,Box_dim,max_long_dist)> % This function tries to assign all bonds to a Bond_matrix and a Bond_index variable.
% # <bond_matrix_atom.html bond_matrix_atom(atom,Box_dim)> % This function generates a bond matrix for the atom struct.
% # <bond_valence_atom.html bond_valence_atom(atom,Box_dim,varargin)> % Calculate bond valence values using the bond valence method.
% # <bond_valence_data.html bond_valence_data(ion1,ion2,R,varargin)> % Fetch data to calculate bond valence values for specified atom types.
% # <bonded_parameters.html bonded_parameters(atom,varargin)> % Define bonded parameters for atoms.
% # <buckinghamcoul.html buckinghamcoul(atom,Box_dim,varargin)> % Calculates interactions using the Buckingham potential and Coulombic forces.
% # <buckinghamcoul_objective_func.html buckinghamcoul_objective_func(atom,Box_dim,varargin)> % Objective function for fitting Buckingham and Coulomb potentials.
% # <cat_atom.html cat_atom(atom_1,atom_2)> % Concatenate two atom structs.
% # <cell_list_dist_matrix_atom.html cell_list_dist_matrix_atom(atom,Box_dim,varargin)> % This function calculates the distance matrix from the atom struct, using a cell list algorithm adapted from the Matlab MDtoolbox by Yasuhiro Matsunaga.
% # <center_atom.html center_atom(atom,Box_dim,resname,dim)> % Center the atom with respect to the resname molecule.
% # <change_top.html change_top(atom,Box_dim,varargin)> % Modifies the topology file or its parameters.
% # <charge_atom.html charge_atom(atom,Box_dim,ffname,watermodel,varargin)> % Charge the atom according to Clayff or Interface force fields.
% # <charge_clayff_2004_atom.html charge_clayff_2004_atom(atom,Box_dim,varargin)> % Set the charge for the original Clayff atom types from the Cygan et al., 2004 paper.
% # <charge_clayff_atom.html charge_clayff_atom(atom,Box_dim,varargin)> % Set the charge for Clayff atom types.
% # <charge_interface15_atom.html charge_interface15_atom(atom,Box_dim,varargin)> % Set the charge for Interface 1.5 atom types.
% # <charge_interface_atom.html charge_interface_atom(atom,Box_dim,varargin)> % Set the charge for Interface atom types.
% # <charge_minff_atom.html charge_minff_atom(atom,Box_dim,varargin)> % Set the charge for MINFF atom types.
% # <charge_opls_go_atom.html charge_opls_go_atom(atom,Box_dim,varargin)> % Set the charge for specific OPLS atom types.
% # <check_H2Odens.html check_H2Odens(atom,Box_dim)> % Compute the water density.
% # <check_clayff_2004_charge.html check_clayff_2004_charge(atom)> % Check the charge of the original Clayff atom types.
% # <check_clayff_H2Odens.html check_clayff_H2Odens(atom,Box_dim)> % Check the approximate water density for a Clayff system.
% # <check_clayff_charge.html check_clayff_charge(atom)> % Check the charge of the Clayff atom types.
% # <check_interface15_charge.html check_interface15_charge(atom)> % Check the charge of the INTERFACE 1.5 atom types.
% # <check_interface_charge.html check_interface_charge(atom)> % Check the charge of the INTERFACE atom types.
% # <clayff210_atom.html clayff210_atom(atom,Box_dim,varargin)> % Assign modified Clayff atom types, with edge healing.
% # <clayff211_atom.html clayff211_atom(atom,Box_dim,varargin)> % Assign modified Clayff atom types (faster version), with edge healing.
% # <clayff_2004_atom.html clayff_2004_atom(atom,Box_dim,varargin)> % Assign the original Clayff atom types by Cygan et al., 2004, with edge healing.
% # <clayff_2004_param.html clayff_2004_param(Atom_label,varargin)> % Hold ion and Clayff atom type parameters for the original Clayff force field.
% # <clayff_atom.html clayff_atom(atom,Box_dim,varargin)> % Assign Clayff atom types with edge healing.
% # <clayff_param.html clayff_param(Atom_label,varargin)> % Hold ion and Clayff atom type parameters.
% # <closest_atom.html closest_atom(atom,Box_dim,ref_atom)> % Find the closest atom to the reference atom.
% # <closest_atom.html closest_atom(atom1,atom2,Box_dim)> % This function returns the atom1 struct with the nMolId's in atom1 closest to the atom2 struct.
% # <cn_atom.html cn_atom(atom,Box_dim,rmax)> % This function calculates the coordination number of atoms within a specified radius.
% # <composition_atom.html composition_atom(atom)> % Analyze the composition of the atom struct.
% # <condense_atom.html condense_atom(atom,Box_dim,s)> % Minimize the box size and remove gaps between molids.
% # <copy_atom.html copy_atom(atom,atomtype,new_atomtype,new_resname,trans_vec,varargin)> % Copy and translate atoms in the atom struct.
% # <copy_type.html copy_type(atom,atomtype,new_atomtype,new_resname,trans_vec,varargin)> % Copy and translate types in the atom struct.
% # <create_atom.html create_atom(type,resname,limits,nmax,varargin)> % Create new atoms, useful for adding ions to a system.
% # <create_grid_atom.html create_grid_atom(atom_label,nM,limits,dim,varargin)> % Put ions on a grid plane and add them to an atom struct.
% # <density_atom.html density_atom(atom,Box_dim)> % Calculate concentration, electron density profiles, and charge density.
% # <density_atom_fast.html density_atom_fast(atom,Box_dim)> % A faster version of the density_atom function.
% # <dist_matrix_atom.html dist_matrix_atom(atom,Box_dim)> % Calculate the distance matrix from the atom struct.
% # <dist_matrix_noPBC_atom.html dist_matrix_noPBC_atom(atom,Box_dim)> % This function calculates the distance matrix without applying periodic boundary conditions.
% # <dist_matrix_xyz.html dist_matrix_xyz(XYZ,Box_dim)> % This function calculates the distance matrix from XYZ coordinates.
% # <distance_factor_variable.html distance_factor> % variable related to finding the nearest neighbours or bonds based on different atomtypes vdw radii.
% # <draw_box_atom.html draw_box_atom(Box_dim,LineColor,LineThickness)> % Draw a box to visualize dimensions.
% # <duplicate_atom.html duplicate_atom(atom,molID)> % Duplicate residue with molid MolID.
% # <element_atom.html element_atom(atom,varargin)> % Convert atom type names to element names.
% # <element_color.html element_color(Atom_label)> % Assign a specific color to each element.
% # <export_ndx.html export_ndx(atom,Box_dim,filename_out)> % Export an index file (.ndx) from the atom struct.
% # <find_H2O.html find_H2O(atom)> % This function identifies and returns water molecules (H2O) within the atom struct.
% # <find_bonded_atom.html find_bonded_atom(atom,bond_matrix,label1,label2)> % This function performs a cross-check of the bond matrix.
% # <find_pair_atom.html find_pair_atom(atom,bond_matrix,pair1,pair2)> % This function finds and returns specific atom pairs from the bond matrix.
% # <fit2lattice_atom.html fit2lattice_atom(atom,Box_dim)> % Fit atoms to a lattice.
% # <fit2lattice_atom_v2.html fit2lattice_atom_v2(atom,Box
% # <frac2atom.html frac2atom(atom,Box_dim,angleparam,angletype)> % Transform fractional coordinates to Cartesian coordinates.
% # <frame2atom.html frame2atom(atom,traj,frame,Box_dim,varargin)> % Extract a frame to the trajectory matrix.
% # <fuse_atom.html fuse_atom(atom,Box_dim,varargin)> % Fuse all sites within a certain cutoff distance.
% # <harmonic_bond.html harmonic_bond(atom,Box_dim,varargin)> % Defines or computes harmonic bond forces or parameters.
% # <heal_atom.html heal_atom(atom,Box_dim,ind,varargin)> % Heal sites in the atom struct by adding a certain atom type.
% # <hist_atom.html hist_atom(atom,s)> % Calculate density profiles in the X, Y, or Z direction.
% # <histz_atom.html histz_atom(atom,s)> % Calculate density profiles in the Z direction.
% # <import_atom.html import_atom(filename)> % Import a .xyz, .gro, or .pdb file into a structure variable called atom.
% # <import_atom_car.html import_atom_car(filename,varargin)> % Import .car files from Hendrik Heinz INTERFACE force field distribution, then write out a Gromacs molecular topology file (.itp) and a new .pdb file.
% # <import_atom_gro.html import_atom_gro(filename)> % Import .gro files into the atom struct.
% # <import_atom_gro_fscanf.html import_atom_gro_fscanf(filename)> % Import .gro files using fscanf (alternative method).
% # <import_atom_gro_octave.html import_atom_gro_octave(filename)> % Import .gro files using Octave-compatible code.
% # <import_atom_mol2.html import_atom_mol2(filename)> % Import .mol2 files into the atom struct.
% # <import_atom_pdb.html import_atom_pdb(filename)> % Import .pdb files into the atom struct.
% # <import_atom_poscar.html import_atom_poscar(filename)> % Import a VASP POSCAR file into the atom struct.
% # <import_atom_pqr.html import_atom_pqr(filename)> % Import .pqr files into the atom struct.
% # <import_atom_xyz.html import_atom_xyz(filename)> % Import an .xyz file into the atom struct.
% # <import_ave_gro.html import_ave_gro(filename)> % Import an averaged structure from a .gro trajectory.
% # <import_bar.html import_bar(filename)> % Import a BAR file.
% # <import_cp2k.html import_cp2k(filename)> % Import a CP2K output file.
% # <import_cp2k_resp.html import_cp2k_resp(filename)> % Import RESP charges from CP2K.
% # <import_dat.html import_dat(filename)> % Import a .dat file.
% # <import_ddec_charges.html import_ddec_charges(filename)> % Import DDEC charges from a file.
% # <import_ff_table.html import_ff_table(filename,varargin)> % Imports forcefield parameter tables.
% # <import_gmx_energy.html import_gmx_energy(filename)> % Import a Gromacs energy file.
% # <import_gro_traj.html import_gro_traj(filename,varargin)> % Import a structure and a .gro trajectory file.
% # <import_mc_pdb_traj.html import_mc_pdb_traj(filename,varargin)> % Import a structure and a .pdb trajectory file, handling changing numbers of particles.
% # <import_mclf.html import_mclf(filename)> % Import a multi-configurational London force (MCLF) file.
% # <import_mclf_C10.html import_mclf_C10(filename)> % Import MCLF C10 dispersion parameters.
% # <import_mclf_C6.html import_mclf_C6(filename)> % Import MCLF C6 dispersion parameters.
% # <import_mclf_C8.html import_mclf_C8(filename)> % Import MCLF C8 dispersion parameters.
% # <import_mclf_dispersion.html import_mclf_dispersion(filename)> % Import dispersion parameters for MCLF.
% # <import_pdb_traj.html import_pdb_traj(filename,varargin)> % Import a structure and a .pdb trajectory file.
% # <import_red_charges.html import_red_charges(filename)> % Import reduced charges from a file.
% # <import_traj.html import_traj(filenameconf,filenametraj)> % Import a structure and a .dcd, .trr, .xtc, .xyz, or .gro trajectory file.
% # <import_trr.html import_trr(filenameconf,filenametraj)> % Import a structure and a .trr trajectory file.
% # <import_trrv2.html import_trrv2(filenameconf,filenametraj)> % Import a structure and a .trr trajectory file (version 2).
% # <import_xtc.html import_xtc(filenameconf,filenamextc)> % Import a structure and an .xtc file.
% # <import_xtcv2.html import_xtcv2(filenameconf,filenamextc)> % Import a structure and an .xtc file (version 2).
% # <import_xvg.html import_xvg(filename)> % Import a Gromacs .xvg file.
% # <import_xyz.html import_xyz(filename)> % Import an .xyz file. Atom types should be made of letters, not numbers. Use import_atom_xyz instead.
% # <import_xyz_traj.html import_xyz_traj(filenametraj)> % Import a structure and an .xyz trajectory file.
% # <insert_atom.html insert_atom(atom,new_atom,position)> % Insert a new atom at the specified position.
% # <interface15_atom.html interface15_atom(atom,Box_dim,varargin)> % Assign atoms according to the Interface 1.5 atom types, with modifications for edges.
% # <interface15_param.html interface15_param(Atom_label,water_model)> % Hold extended INTERFACE 1.5 force field parameters.
% # <interface15_silica_atom.html interface15_silica_atom(atom,Box_dim,varargin)> % Assign atom types for the Interface 1.5 force field, specific to silica.
% # <interface_atom.html interface_atom(atom,Box_dim,varargin)> % Assign atoms according to the Interface atom types, with modifications for edges.
% # <interface_param.html interface_param(Atom_label,water_model)> % Hold extended INTERFACE force field parameters.
% # <ionize_atom.html ionize_atom(type,resname,limits,nmax,varargin)> % Add ions within a certain region defined by limits.
% # <keep_atom.html keep_atom(atom,resname)> % Keep only specified residue names in the atom struct.
% # <keep_resname.html keep_resname(atom,resnames)> % Keep only specified residue names.
% # <limits_variable.html limits> % A 1x6 array defining a volumetric region
% # <list_bonds_atom.html list_bonds_atom(atom,bond_matrix)> % This function lists all bonds in the atom struct based on the bond matrix.
% # <ljcoul.html ljcoul(atom,Box_dim,varargin)> % General Lennard-Jones and Coulomb interaction function.
% # <ljcoul_12_6.html ljcoul_12_6(atom,Box_dim,varargin)> % Handles Lennard-Jones 12-6 potential along with Coulombic interactions.
% # <ljcoul_2x.html ljcoul_2x(atom,Box_dim,varargin)> % Variation of Lennard-Jones Coulomb with 2x factors.
% # <ljcoul_2x2x.html ljcoul_2x2x(atom,Box_dim,varargin)> % Lennard-Jones and Coulomb potential with 2x factors.
% # <ljcoul_2x_force.html ljcoul_2x_force(atom,Box_dim,varargin)> % Force calculation involving 2x Lennard-Jones and Coulomb potentials.
% # <ljcoul_2x_objective_func.html ljcoul_2x_objective_func(atom,Box_dim,varargin)> % Objective function for 2x Lennard-Jones and Coulomb potentials.
% # <ljcoul_C12C6.html ljcoul_C12C6(atom,Box_dim,varargin)> % Lennard-Jones and Coulomb potential with C12, C6 terms.
% # <ljcoul_C12C6C4.html ljcoul_C12C6C4(atom,Box_dim,varargin)> % Lennard-Jones and Coulomb potential with C12, C6, C4 terms.
% # <ljcoul_force.html ljcoul_force(atom,Box_dim,varargin)> % Computes forces based on Lennard-Jones and Coulomb potentials.
% # <ljcoul_force_C12C6C4.html ljcoul_force_C12C6C4(atom,Box_dim,varargin)> % Force calculation for Lennard-Jones potential with C12, C6, C4 terms.
% # <ljcoul_force_objective_func.html ljcoul_force_objective_func(atom,Box_dim,varargin)> % Objective function for force calculations with Lennard-Jones and Coulomb potentials.
% # <ljcoul_objective_func.html ljcoul_objective_func(atom,Box_dim,varargin)> % Objective function for Lennard-Jones and Coulomb potentials.
% # <mass_atom.html mass_atom(atom)> % Calculate or define mass for atoms.
% # <mass_atom_clayff.html mass_atom_clayff(atom)> % Defines or calculates masses for Clayff atoms.
% # <median_atom.html median_atom(atom)> % Calculate the median position of the atom struct.
% # <merge_atom.html merge_atom(atom1,Box1,atom2,type,Atom_label,r)> % Merge atom structs based on distance criteria.
% # <minff_atom.html minff_atom(atom,Box_dim,varargin)> % Assign minff atom types with edge healing.
% # <molecule_atom.html molecule_atom(atom,varargin)> % Set molecule ID, residue name, and element names of the atom struct.
% # <molid_rotate.html molid_rotate(atom,Box_dim,MolID,rotate_dim)> % Rotate the atom struct based on MolID.
% # <molid_translate.html molid_translate(atom,trans_vec,MolID)> % Translate a specific molecule ID.
% # <neigh_atom.html neigh_atom(atom,Box_dim,rmax,varargin)> % Check neighbors for each atom and output their information.
% # <neighbor_atom.html neighbor_atom(atom,Box_dim,radius)> % This function checks the neighbors of each atom and returns their indices.
% # <neighbor_func.html neighbor_func(solute_index,XYZ_solute,XYZ_data,Box_dim,radius)> % Scan xyz data and check neighbors within a specified radius.
% # <neutralize_atom.html neutralize_atom(atom)> % Set the charge of all atom types to zero.
% # <nonbonded_ff.html nonbonded_ff(atom,varargin)> % Defines nonbonded forcefield parameters.
% # <nonbonded_parameters.html nonbonded_parameters(atom,varargin)> % Define nonbonded parameters for atoms.
% # <noupdate_atom.html noupdate_atom(atom)> % Prevent updating of certain properties in the atom struct.
% # <occupancy_atom.html occupancy_atom(atom,Box_dim)> % Calculate occupancy of atoms within the box dimensions.
% # <opls_go_atom.html opls_go_atom(atom,Box_dim,rmin,rlarge)> % Smear out the charge around -OH and epoxide groups in graphene oxide.
% # <oplsaa_go_param.html oplsaa_go_param(Atom_label,water_model)> % Hold the extended OPLS-AA force field parameters for graphite oxide.
% # <orto_atom.html orto_atom(atom,Box_dim)> % Transform a triclinic atom struct to an orthogonal one.
% # <overwrite_atom.html overwrite_atom(In_atom,atomtype,resname)> % Overwrite atom struct information with new data.
% # <place_atom.html place_atom(atom,position)> % Place the atom struct at a specified position.
% # <plot_atom.html plot_atom(atom,Box_dim,varargin)> % Draw the atom struct in 3D.
% # <plot_density_atom.html plot_density_atom(atom,Box_dim,varargin)> % Draw the atom struct in 3D along with density profiles.
% # <plot_xvg.html plot_xvg(filename)> % Plot .xvg files, typically used in GROMACS simulations.
% # <position_molid.html position_molid(atom,position_vec,MolID)> % Move a molecule ID to a certain position.
% # <print_top.html print_top(atom,Box_dim,varargin)> % Prints or generates topology-related data.
% # <properties_atom.html properties_atom(atom)> % Analyze properties of the atom struct.
% # <protonate_atom.html protonate_atom(atom,Box_dim,varargin)> % Protonate specified sites in the atom struct.
% # <radius_atom.html radius_atom(atom,ffname,watermodel)> % Fetch the ion radius from Clayff, Interface, or Interface2015 force fields.
% # <radius_crystal.html radius_crystal(Atom_label)> % Fetch or calculate ionic radii for crystal structures.
% # <radius_ion.html radius_ion(Atom_label)> % Fetch the ionic radius.
% # <radius_vdw.html radius_vdw(Atom_label)> % Fetch the van der Waals radius.
% # <rdf_atom.html rdf_atom(atom,Box_dim,varargin)> % Calculate the radial distribution function and coordination number.
% # <recalc_bond_atom.html recalc_bond_atom(atom,bond_matrix,varargin)> % This function recalculates bonds for the atom struct.
% # <reduced_mass.html reduced_mass(Atom_label1,varargin)> % Calculate the reduced mass.
% # <remove_H2O.html remove_H2O(atom)> % This function removes water molecules (H2O) from the atom struct.
% # <remove_SOL.html remove_SOL(atom,atomname,lo,hi,dim)> % Remove solvent residues between specified limits.
% # <remove_duplicate_atom.html remove_duplicate_atom(atom)> % Remove duplicate atoms from the atom struct.
% # <remove_molid.html remove_molid(atom,MolID)> % Remove residue with a specific molecule ID.
% # <remove_occypancy_atom.html remove_occypancy_atom(atom)> % Remove particles with identical coordinates to preceding ones.
% # <remove_residues.html remove_residues(atom,resnames,lo,hi,dim)> % Remove residues between specified limits in the simulation box.
% # <remove_resname.html remove_resname(atom,resnames)> % Remove residues by name.
% # <remove_sametype_bond.html remove_sametype_bond(atom,bond_matrix)> % This function removes bonds between atoms of the same type.
% # <remove_type.html remove_type(atom,typescell)> % Remove atom types specified in typescell.
% # <rename_atom.html rename_atom(atom,old_name,new_name)> % Rename atoms in the atom struct.
% # <rename_type.html rename_type(atom,atomtype,new_atomtype,varargin)> % Rename atom types in the atom struct.
% # <replace_atom.html replace_atom(new_atom,prev_atom,molid_index)> % Replace molecule ID in an atom struct with a new atom struct.
% # <replace_string.html replace_string(filename_in,filename_out,old_string,new_string)> % Replace strings in files.
% # <replicate_atom.html replicate_atom(atom,Box_dim,replicate)> % Replicate the atom struct along orthogonal dimensions.
% # <resname_atom.html resname_atom(atom)> % Guess residue names for all atom types.
% # <rmaxshort_variable.html rmaxshort> % Maximum H-related bond radius
% # <rmaxslong_variable.html rmaxlong> % Maximum non-H-related bond/neighbor radius
% # <rotate_atom.html rotate_atom(atom,rotation_matrix)> % Rotate the atom struct by a specified matrix.
% # <round_atom.html round_atom(atom)> % Round atom positions or coordinates.
% # <show_Hbonds_atom.html show_Hbonds_atom(atom)> % Display or calculate hydrogen bonds in the atom struct.
% # <show_arrow.html show_arrow(p1,p2,varargin)> % Plot a 3D arrow as a patch object.
% # <show_atom.html show_atom(atom,varargin)> % Draw the atom struct in 3D with additional features.
% # <show_axis.html show_axis(varargin)> % Draw the axis in a plot.
% # <show_box.html show_box(Box_dim)> % Draw the simulation box.
% # <show_density_atom.html show_density_atom(atom)> % Display density of atoms.
% # <show_miller.html show_miller(Box_dim)> % Draw the Miller planes of the Box_dim/Cell variables.
% # <sigma_vdw.html sigma_vdw(Atom_label)> % Compute sigma values for van der Waals interactions.
% # <slice_atom.html slice_atom(atom,limits,invert)> % Slice the atom struct within specified limits.
% # <slice_box.html slice_box(atom,Box_dim,limits)> % Slice a simulation box within given limits.
% # <slice_molid.html slice_molid(atom,limits,invert)> % Slice molecules within specified limits.
% # <slice_triclinic_atom.html slice_triclinic_atom(atom,limits,invert)> % Slice a triclinic atom struct within limits.
% # <smear_charge.html smear_charge(atom,Box_dim,varargin)> % Distributes charge across atoms, possibly using charge smearing techniques.
% # <solvate_atom.html solvate_atom(limits,density,r,maxsol,solute_atom,varargin)> % Generate a solvent structure within specified limits.
% # <sort_atom.html sort_atom(atom,Box_dim,varargin)> % Sort the atom struct based on specified properties.
% # <spc2tip4p.html spc2tip4p(atom)> % Convert SPC water molecules to TIP4P model.
% # <spc2tip4p.html spc2tip4p(filename)> % Convert a .gro or .pdb file with spc water to tip4p water.
% # <spc2tip5p.html spc2tip5p(atom)> % Convert SPC water molecules to TIP5P model.
% # <spc2tip5p.html spc2tip5p(filename)> % Convert a .gro or .pdb file with spc water to tip5p water.
% # <spce2tip4p.html spce2tip4p(atom)> % Convert SPC/E water molecules to TIP4P model.
% # <spce2tip4p.html spce2tip4p(filename)> % Convert a .gro or .pdb file with spce water to tip4p water.
% # <sphere_atom.html sphere_atom(atom,Box_dim,center,radius)> % Create a spherical region of atoms.
% # <spiral_atom.html spiral_atom(atom,Box_dim,varargin)> % Create or modify spiral structures in the atom struct.
% # <split_atom.html split_atom(atom,Box_dim,varargin)> % Split the atom struct into separate parts based on criteria.
% # <substitute_NonCentroSymm_atom.html substitute_NonCentroSymm_atom(atom,Box_dim,replace_type,varargin)> % Substitute non-centrosymmetric atoms.
% # <substitute_atom.html substitute_atom(atom,Box_dim,NumOctSubst,O1,O2,minO2O2_dist,varargin)> % Perform isomorphous substitution of atoms.
% # <tile_atom.html tile_atom(atom,scale_vec,Box_dim,Resname)> % Tile the atom struct in a specific direction.
% # <tip3p2tip4p.html tip3p2tip4p(atom)> % Convert TIP3P water molecules to TIP4P model.
% # <tip3p2tip4p.html tip3p2tip4p(filename)> % Convert a .gro file with tip3p water to tip4p water.
% # <translate_atom.html translate_atom(atom,trans_vec)> % Translate the atom struct by a specified vector.
% # <translate_atom.html translate_atom(atom,trans_vec,Resname)> % Translate a residue by a specified vector.
% # <translate_molid.html translate_molid(atom,trans_vec,molid)> % Translate a molecule ID by a specified vector.
% # <triclinic_atom.html triclinic_atom(atom,Box_dim,angleparam,angletype)> % Transform an orthogonal atom struct to a triclinic one.
% # <tube_atom.html tube_atom(atom,scale_vec,Box_dim,Resname)> % Create a nanotube structure from the atom struct.
% # <tweak_charge_atom.html tweak_charge_atom(atom)> % Tweak the charge of the atom struct to correct rounding errors.
% # <unreplicate_atom.html unreplicate_atom(atom)> % Remove replicated atoms from the structure.
% # <unwrap_atom.html unwrap_atom(atom,Box_dim,dim)> % Unwrap the atom struct along the specified dimension.
% # <update_atom.html update_atom(atom)> % Update molecule and atom indices in the atom struct.
% # <vmd.html vmd(atom,Box_dim)> % Plot the atom struct using VMD.
% # <write_atom.html write_atom(atom,Box_dim,filename_out,varargin)> % Write different file types (.gro, .pdb, .xyz, .itp, etc.) based on filename and parameters.
% # <write_atom_all.html write_atom_all(atom,Box_dim,filename_out,varargin)> % Write various file types for the atom struct, best suited for Clayff systems.
% # <write_atom_cif.html write_atom_cif(atom,Box_dim,filename_out)> % Write a basic .cif file from the atom struct.
% # <write_atom_dodecahedron_gro.html write_atom_dodecahedron_gro(atom,Box_dim,filename_out)> % Write a .gro file using a dodecahedron-shaped simulation box.
% # <write_atom_gro.html write_atom_gro(atom,Box_dim,filename_out)> % Write a .gro file from the atom struct, optionally including velocities.
% # <write_atom_itp.html write_atom_itp(atom,Box_dim,filename_out,varargin)> % Create and print a Gromacs .itp file for Clayff or Interface force fields.
% # <write_atom_lmp.html write_atom_lmp(atom,Box_dim,filename_out,varargin)> % Create and print a LAMMPS data file (.lj) for Clayff systems.
% # <write_atom_mol2.html write_atom_mol2(atom,Bond_index,Box_dim,filename_out)> % Write a .mol2 file from the atom struct.
% # <write_atom_multiple_gro.html write_atom_multiple_gro(atom,traj,filename_out)> % Write multiple .gro files for trajectory output.
% # <write_atom_oplsaa_go_itp.html write_atom_oplsaa_go_itp(atom,Box_dim,filename_out,varargin)> % Create and print a Gromacs .itp file for OPLS-AA or GO systems.
% # <write_atom_pdb.html write_atom_pdb(atom,Box_dim,filename_out)> % Write a .pdb file from the atom struct using Gromacs.
% # <write_atom_pqr.html write_atom_pqr(atom,Box_dim,filename_out,varargin)> % Write a .pqr file from the atom struct.
% # <write_atom_psf.html write_atom_psf(atom,Box_dim,filename_out,varargin)> % Write a .psf file from the atom struct.
% # <write_atom_sdf.html write_atom_sdf(atom,Box_dim,filename_out)> % Write an SDF file from the atom struct.
% # <write_atom_top.html write_atom_top(atom,Box_dim,filename_out)> % Write a topology file (.top) from the atom struct.
% # <write_atom_xyz.html write_atom_xyz(atom,Box_dim,filename_out)> % Write an XYZ file from the atom struct.
% # <write_ave_gro.html write_ave_gro(atom,traj,Box_dim,filename_out)> % Write an average structure from a .gro trajectory.
% # <write_ave_pdb.html write_ave_pdb(atom,traj,Box_dim,filename_out)> % Write an average structure from a .pdb trajectory.
% # <write_ff.html write_ff(atom,filename_out)> % Write force field parameters.
% # <write_ffnonbonded.html write_ffnonbonded(atom,filename_out)> % Write non-bonded force field parameters.
% # <write_ffnonbonded_C6C12.html write_ffnonbonded_C6C12(atom,filename_out)> % Write non-bonded parameters (C6, C12) for a force field.
% # <write_gro_traj.html write_gro_traj(atom,traj,Box_dim,filename_out)> % Write a .gro trajectory file.
% # <write_pdb_traj.html write_pdb_traj(atom,traj,Box_dim,filename_out)> % Write a .pdb trajectory file.
% # <write_tabulated_potentials.html write_tabulated_potentials(filename, data)> % Write tabulated potential files.
% # <write_xvg.html write_xvg(filename, data)> % Export data in Gromacs .xvg format.
% # <write_xyz_traj.html write_xyz_traj(atom,traj,Box_dim,filename_out)> % Write a .xyz trajectory file.
% # <xrd_atom.html xrd_atom(varargin)> % Calculate theoretical XRD patterns from a .pdb, .gro file, or atom struct.
%

%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se