# Atom MATLAB Library

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![MATLAB Version](https://img.shields.io/badge/MATLAB-2024-orange.svg)

## Overview

The **Atom MATLAB Library** is a comprehensive collection of MATLAB functions designed for handling and analyzing atomic structures, molecular simulations, and related computational tasks. Whether you're importing data from various file formats, performing structural manipulations, calculating distances and bonds, or exporting data for simulation tools like GROMACS and LAMMPS, this library provides the necessary tools to streamline your workflow.

## Version

**3.00**

## Contact

For any issues, bugs, or feature requests, please contact:

📧 [michael.holmboe@umu.se](mailto:michael.holmboe@umu.se)

## Table of Contents

- [Main Types of Variables](#main-types-of-variables)
- [Import Functions](#import-functions)
  - [Data Files](#import-of-data-files)
  - [Structure Files](#import-of-structure-files)
  - [Trajectory Files](#import-of-trajectory-files)
  - [Miscellaneous Import Functions](#miscellaneous-import-functions)
- [Export Functions](#export-functions)
  - [Structure/Topology Files](#export-of-structuretopology-files)
  - [Trajectory Files](#export-of-trajectory-files)
  - [Other Formats and Data](#export-of-other-formats-and-data)
  - [Miscellaneous Export Functions](#miscellaneous-export-functions)
- [Neighbor/Distance Functions](#neighbordistance-functions)
- [Specific Atom Struct Functions](#specific-atom-struct-functions)
- [Add/Create/Replicate/Overwrite Atoms](#addcreatereplicateoverwrite-atoms)
- [Slice Out a Region of the Box](#slice-out-a-region-of-the-box)
- [Translate or Rotate Functions](#translate-or-rotate-functions)
- [Make Triclinic/Orthogonal Box](#make-triclinicorthogonal-box)
- [Wrap/Unwrap Functions](#wrapunwrap-functions)
- [Keep/Remove Functions](#keeprm-functions)
- [List of Available Solvents](#list-of-available-solvents)
- [Conversion Functions](#conversion-functions)
- [Custom Topology Tools](#custom-topology-tools)
- [Miscellaneous FF Functions](#miscellaneous-ff-functions)
- [Force Fields](#force-fields)
  - [MINFF](#minff-with-atomtypes-by-mholmboe)
  - [CLAYFF](#clayff-with-atomtypes-by-mholmboe)
  - [INTERFACE](#interface-from-heinz-2005-2013-with-atomtypes-by-mholmboe)
  - [Graphene Oxide with OPLS/aa](#graphene-oxide-modeled-with-oplsaa)
- [Writing Topology Files](#writing-topology-files)
- [Bonded and Nonbonded Parameters](#bonded-and-nonbonded-parameters)
- [Lennard-Jones and Coulomb Potentials](#lennard-jones-and-coulomb-potentials)
- [Objective Functions and Force Calculations](#objective-functions-and-force-calculations)
- [Automated Fitting Tools](#automated-fitting-tools)

---

## Main Types of Variables

- **[Atom_label](Atom_label_variable.html)**: `1xn cell array`
- **[atom](atom_variable.html)**: The main MATLAB struct variable
- **[Box_dim](Box_dim_variable.html)**: `1x3` or `1x9` array holding the box size parameters
- **[Cell](Cell_variable.html)**: `1x6` array containing a, b, c cell values and alpha, beta, gamma angles (as in a typical `.pdb` file)
- **[distance_factor](distance_factor_variable.html)**: Variable for finding nearest neighbors or bonds based on different atom types' VDW radii
- **[limits](limits_variable.html)**: `1x6` array defining a volumetric region
- **[rmaxshort](rmaxshort_variable.html)**: Maximum H-related bond radius
- **[rmaxlong](rmaxslong_variable.html)**: Maximum non-H-related bond/neighbor radius
- **[XYZ_data](XYZ_data_variable.html)**: `nx3` matrix holding XYZ coordinates
- **[XYZ_labels](XYZ_labels_variable.html)**: Cell list of atom types

---

## Import Functions

### Import of Data Files

- **[import_xvg](import_xvg.html)**: `import_xvg(filename)`  
  Import a GROMACS `.xvg` file.

- **[import_gmx_energy](import_gmx_energy.html)**: `import_gmx_energy(filename)`  
  Import a GROMACS energy file.

- **[import_bar](import_bar.html)**: `import_bar(filename)`  
  Import a BAR file.

- **[import_dat](import_dat.html)**: `import_dat(filename)`  
  Import a `.dat` file.

- **[import_red_charges](import_red_charges.html)**: `import_red_charges(filename)`  
  Import reduced charges from a file.

- **[import_ddec_charges](import_ddec_charges.html)**: `import_ddec_charges(filename)`  
  Import DDEC charges from a file.

### Import of Structure Files

- **[import_atom_car](import_atom_car.html)**: `import_atom_car(filename, varargin)`  
  Import `.car` files from Hendrik Heinz INTERFACE force field distribution and write out GROMACS `.itp` and `.pdb` files.

- **[import_atom_gro](import_atom_gro.html)**: `import_atom_gro(filename)`  
  Import `.gro` files into the atom struct.

- **[import_atom_gro_fscanf](import_atom_gro_fscanf.html)**: `import_atom_gro_fscanf(filename)`  
  Import `.gro` files using `fscanf` (alternative method).

- **[import_atom_gro_octave](import_atom_gro_octave.html)**: `import_atom_gro_octave(filename)`  
  Import `.gro` files using Octave-compatible code.

- **[import_atom_mol2](import_atom_mol2.html)**: `import_atom_mol2(filename)`  
  Import `.mol2` files into the atom struct.

- **[import_atom_pdb](import_atom_pdb.html)**: `import_atom_pdb(filename)`  
  Import `.pdb` files into the atom struct.

- **[import_atom_pqr](import_atom_pqr.html)**: `import_atom_pqr(filename)`  
  Import `.pqr` files into the atom struct.

- **[import_atom_xyz](import_atom_xyz.html)**: `import_atom_xyz(filename)`  
  Import an `.xyz` file into the atom struct.

- **[import_atom_poscar](import_atom_poscar.html)**: `import_atom_poscar(filename)`  
  Import a VASP POSCAR file into the atom struct.

- **[import_atom](import_atom.html)**: `import_atom(filename)`  
  Import a `.xyz`, `.gro`, or `.pdb` file into a structure variable called `atom`.

- **[import_xyz](import_xyz.html)**: `import_xyz(filename)`  
  Import an `.xyz` file. Atom types should be letters, not numbers. Use `import_atom_xyz` instead.

### Import of Trajectory Files

- **[import_gro_traj](import_gro_traj.html)**: `import_gro_traj(filename, varargin)`  
  Import a structure and a `.gro` trajectory file.

- **[import_xyz_traj](import_xyz_traj.html)**: `import_xyz_traj(filenametraj)`  
  Import a structure and an `.xyz` trajectory file.

- **[import_mc_pdb_traj](import_mc_pdb_traj.html)**: `import_mc_pdb_traj(filename, varargin)`  
  Import a structure and a `.pdb` trajectory file, handling changing numbers of particles.

- **[import_pdb_traj](import_pdb_traj.html)**: `import_pdb_traj(filename, varargin)`  
  Import a structure and a `.pdb` trajectory file.

- **[import_traj](import_traj.html)**: `import_traj(filenameconf, filenametraj)`  
  Import a structure and a `.dcd`, `.trr`, `.xtc`, `.xyz`, or `.gro` trajectory file.

- **[import_trr](import_trr.html)**: `import_trr(filenameconf, filenametraj)`  
  Import a structure and a `.trr` trajectory file.

- **[import_trrv2](import_trrv2.html)**: `import_trrv2(filenameconf, filenametraj)`  
  Import a structure and a `.trr` trajectory file (version 2).

- **[import_xtc](import_xtc.html)**: `import_xtc(filenameconf, filenamextc)`  
  Import a structure and an `.xtc` file.

- **[import_xtcv2](import_xtcv2.html)**: `import_xtcv2(filenameconf, filenamextc)`  
  Import a structure and an `.xtc` file (version 2).

### Miscellaneous Import Functions

- **[import_ave_gro](import_ave_gro.html)**: `import_ave_gro(filename)`  
  Import an averaged structure from a `.gro` trajectory.

- **[import_mclf](import_mclf.html)**: `import_mclf(filename)`  
  Import a multi-configurational London force (MCLF) file.

- **[import_mclf_C6](import_mclf_C6.html)**: `import_mclf_C6(filename)`  
  Import MCLF C6 dispersion parameters.

- **[import_mclf_C8](import_mclf_C8.html)**: `import_mclf_C8(filename)`  
  Import MCLF C8 dispersion parameters.

- **[import_mclf_C10](import_mclf_C10.html)**: `import_mclf_C10(filename)`  
  Import MCLF C10 dispersion parameters.

- **[import_mclf_dispersion](import_mclf_dispersion.html)**: `import_mclf_dispersion(filename)`  
  Import dispersion parameters for MCLF.

- **[import_cp2k](import_cp2k.html)**: `import_cp2k(filename)`  
  Import a CP2K output file.

- **[import_cp2k_resp](import_cp2k_resp.html)**: `import_cp2k_resp(filename)`  
  Import RESP charges from CP2K.

---

## Export Functions

### Export of Structure/Topology Files

- **[write_atom_all](write_atom_all.html)**: `write_atom_all(atom, Box_dim, filename_out, varargin)`  
  Write various file types for the atom struct, best suited for Clayff systems.

- **[write_atom_cif](write_atom_cif.html)**: `write_atom_cif(atom, Box_dim, filename_out)`  
  Write a basic `.cif` file from the atom struct.

- **[write_atom_gro](write_atom_gro.html)**: `write_atom_gro(atom, Box_dim, filename_out)`  
  Write a `.gro` file from the atom struct, optionally including velocities.

- **[write_atom_itp](write_atom_itp.html)**: `write_atom_itp(atom, Box_dim, filename_out, varargin)`  
  Create and print a GROMACS `.itp` file, best for Clayff or interface force fields.

- **[write_atom_lmp](write_atom_lmp.html)**: `write_atom_lmp(atom, Box_dim, filename_out, varargin)`  
  Create and print a LAMMPS data file (`.lj`), suited for Clayff systems.

- **[write_atom_mol2](write_atom_mol2.html)**: `write_atom_mol2(atom, Bond_index, Box_dim, filename_out)`  
  Write a `.mol2` file from the atom struct.

- **[write_atom_oplsaa_go_itp](write_atom_oplsaa_go_itp.html)**: `write_atom_oplsaa_go_itp(atom, Box_dim, filename_out, varargin)`  
  Create and print a GROMACS `.itp` file for OPLS-AA or GO systems.

- **[write_atom_pdb](write_atom_pdb.html)**: `write_atom_pdb(atom, Box_dim, filename_out)`  
  Write a `.pdb` file from the atom struct using GROMACS.

- **[write_atom_pqr](write_atom_pqr.html)**: `write_atom_pqr(atom, Box_dim, filename_out, varargin)`  
  Write a `.pqr` file from the atom struct.

- **[write_atom_psf](write_atom_psf.html)**: `write_atom_psf(atom, Box_dim, filename_out, varargin)`  
  Write a `.psf` file from the atom struct.

- **[write_atom_xyz](write_atom_xyz.html)**: `write_atom_xyz(atom, Box_dim, filename_out)`  
  Write an `.xyz` file from the atom struct.

- **[write_atom](write_atom.html)**: `write_atom(atom, Box_dim, filename_out, varargin)`  
  Write different file types (`.gro`, `.pdb`, `.xyz`, `.itp`, etc.) based on filename and parameters.

- **[write_atom_sdf](write_atom_sdf.html)**: `write_atom_sdf(atom, Box_dim, filename_out)`  
  Write an `.sdf` file from the atom struct.

- **[write_atom_dodecahedron_gro](write_atom_dodecahedron_gro.html)**: `write_atom_dodecahedron_gro(atom, Box_dim, filename_out)`  
  Write a `.gro` file using a dodecahedron-shaped simulation box.

- **[write_atom_top](write_atom_top.html)**: `write_atom_top(atom, Box_dim, filename_out)`  
  Write a topology file (`.top`) from the atom struct.

- **[CONECT_atom](CONECT_atom.html)**: `CONECT_atom(atom, Box_dim, filename_out)`  
  Write CONECT records for a PDB file.

### Export of Trajectory Files

- **[write_gro_traj](write_gro_traj.html)**: `write_gro_traj(atom, traj, Box_dim, filename_out)`  
  Write a `.gro` trajectory file.

- **[write_atom_multiple_gro](write_atom_multiple_gro.html)**: `write_atom_multiple_gro(atom, traj, filename_out)`  
  Write multiple `.gro` files for trajectory output.

- **[write_pdb_traj](write_pdb_traj.html)**: `write_pdb_traj(atom, traj, Box_dim, filename_out)`  
  Write a `.pdb` trajectory file.

- **[write_xyz_traj](write_xyz_traj.html)**: `write_xyz_traj(atom, traj, Box_dim, filename_out)`  
  Write a `.xyz` trajectory file.

- **[write_ave_gro](write_ave_gro.html)**: `write_ave_gro(atom, traj, Box_dim, filename_out)`  
  Write an average structure from a `.gro` trajectory.

- **[write_ave_pdb](write_ave_pdb.html)**: `write_ave_pdb(atom, traj, Box_dim, filename_out)`  
  Write an average structure from a `.pdb` trajectory.

### Export of Other Formats and Data

- **[export_ndx](export_ndx.html)**: `export_ndx(atom, Box_dim, filename_out)`  
  Export an index file (`.ndx`) from the atom struct.

- **[write_xvg](write_xvg.html)**: `write_xvg(filename, data)`  
  Export data in GROMACS `.xvg` format.

- **[write_tabulated_potentials](write_tabulated_potentials.html)**: `write_tabulated_potentials(filename, data)`  
  Write tabulated potential files.

- **[write_ff](write_ff.html)**: `write_ff(atom, filename_out)`  
  Write force field parameters.

- **[write_ffnonbonded](write_ffnonbonded.html)**: `write_ffnonbonded(atom, filename_out)`  
  Write non-bonded force field parameters.

- **[write_ffnonbonded_C6C12](write_ffnonbonded_C6C12.html)**: `write_ffnonbonded_C6C12(atom, filename_out)`  
  Write non-bonded parameters (C6, C12) for a force field.

### Miscellaneous Export Functions

- **[replace_string](replace_string.html)**: `replace_string(filename_in, filename_out, old_string, new_string)`  
  Replace strings in files.

---

## Neighbor/Distance Functions

- **[bond_atom](bond_atom.html)**: `bond_atom(atom, Box_dim, max_long_dist)`  
  Assign all bonds to a `Bond_matrix` and `Bond_index` variable.

- **[bond_angle_atom](bond_angle_atom.html)**: `bond_angle_atom(atom, Box_dim, varargin)`  
  Find all bonds and angles of the atom struct. `'More'` is an optional argument.

- **[bond_angle_dihedral_atom](bond_angle_dihedral_atom.html)**: `bond_angle_dihedral_atom(atom, Box_dim, varargin)`  
  Find all bonds, angles, and dihedrals of the atom struct. `Rmaxshort` and `Rmaxlong` as well as `'more'` are optional arguments.

- **[bond_angle_type](bond_angle_type.html)**: `bond_angle_type(atom1, atom2, Box_dim, rmin, rmax, angle_limit, varargin)`  
  Find all bonds and angles of the atom types.

- **[bond_valence_atom](bond_valence_atom.html)**: `bond_valence_atom(atom, Box_dim, varargin)`  
  Calculate bond valence values using the bond valence method.

- **[bond_valence_data](bond_valence_data.html)**: `bond_valence_data(ion1, ion2, R, varargin)`  
  Fetch data to calculate bond valence values for specified atom types.

- **[cell_list_dist_matrix_atom](cell_list_dist_matrix_atom.html)**: `cell_list_dist_matrix_atom(atom, Box_dim, varargin)`  
  Calculate the distance matrix from the atom struct using a cell list algorithm.

- **[closest_atom](closest_atom.html)**: `closest_atom(atom1, atom2, Box_dim)`  
  Return the `atom1` struct with the `nMolId's` in `atom1` closest to the `atom2` struct.

- **[dist_matrix_atom](dist_matrix_atom.html)**: `dist_matrix_atom(atom, Box_dim)`  
  Calculate the distance matrix from the atom struct.

- **[dist_matrix_noPBC_atom](dist_matrix_noPBC_atom.html)**: `dist_matrix_noPBC_atom(atom, Box_dim)`  
  Calculate the distance matrix without applying periodic boundary conditions.

- **[dist_matrix_xyz](dist_matrix_xyz.html)**: `dist_matrix_xyz(XYZ, Box_dim)`  
  Calculate the distance matrix from XYZ coordinates.

- **[find_bonded_atom](find_bonded_atom.html)**: `find_bonded_atom(atom, bond_matrix, label1, label2)`  
  Perform a cross-check of the bond matrix.

- **[find_pair_atom](find_pair_atom.html)**: `find_pair_atom(atom, bond_matrix, pair1, pair2)`  
  Find and return specific atom pairs from the bond matrix.

- **[list_bonds_atom](list_bonds_atom.html)**: `list_bonds_atom(atom, bond_matrix)`  
  List all bonds in the atom struct based on the bond matrix.

- **[neigh_atom](neigh_atom.html)**: `neigh_atom(atom, Box_dim, rmax, varargin)`  
  Check which neighbors each atom has and output their information.

- **[neighbor_func](neighbor_func.html)**: `neighbor_func(solute_index, XYZ_solute, XYZ_data, Box_dim, radius)`  
  Scan XYZ data and check which atoms are within a certain radius. Outputs the neighbor index.

- **[neighbor_atom](neighbor_atom.html)**: `neighbor_atom(atom, Box_dim, radius)`  
  Check the neighbors of each atom and return their indices.

- **[rdf_atom](rdf_atom.html)**: `rdf_atom(atom, Box_dim, varargin)`  
  Calculate the radial distribution function and the coordination number, with optional Gaussian smoothing.

- **[cn_atom](cn_atom.html)**: `cn_atom(atom, Box_dim, rmax)`  
  Calculate the coordination number of atoms within a specified radius.

- **[recalc_bond_atom](recalc_bond_atom.html)**: `recalc_bond_atom(atom, bond_matrix, varargin)`  
  Recalculate bonds for the atom struct.

- **[remove_H2O](remove_H2O.html)**: `remove_H2O(atom)`  
  Remove water molecules (`H2O`) from the atom struct.

- **[remove_sametype_bond](remove_sametype_bond.html)**: `remove_sametype_bond(atom, bond_matrix)`  
  Remove bonds between atoms of the same type.

- **[find_H2O](find_H2O.html)**: `find_H2O(atom)`  
  Identify and return water molecules (`H2O`) within the atom struct.

- **[bond_matrix_atom](bond_matrix_atom.html)**: `bond_matrix_atom(atom, Box_dim)`  
  Generate a bond matrix for the atom struct.

---

## Specific Atom Struct Functions

- **[add2atom](add2atom.html)**: `add2atom(XYZ_labels, XYZ_data, varargin)`  
  Append XYZ atom type labels and XYZ data to an existing atom struct.

- **[adjust_H_atom](adjust_H_atom.html)**: `adjust_H_atom(atom, Box_dim)`  
  Adjust hydrogen atoms in the atom struct.

- **[cat_atom](cat_atom.html)**: `cat_atom(atom_1, atom_2)`  
  Concatenate two atom structs.

- **[closest_atom](closest_atom.html)**: `closest_atom(atom, Box_dim, ref_atom)`  
  Find the closest atom to the reference atom.

- **[copy_atom](copy_atom.html)**: `copy_atom(atom, atomtype, new_atomtype, new_resname, trans_vec, varargin)`  
  Copy and translate atoms in the atom struct.

- **[create_atom](create_atom.html)**: `create_atom(type, resname, limits, nmax, varargin)`  
  Create new atoms, useful for adding ions to a system.

- **[create_grid_atom](create_grid_atom.html)**: `create_grid_atom(atom_label, nM, limits, dim, varargin)`  
  Put ions on a grid plane and add them to an atom struct.

- **[duplicate_atom](duplicate_atom.html)**: `duplicate_atom(atom, molID)`  
  Duplicate residue with `molid` MolID.

- **[fuse_atom](fuse_atom.html)**: `fuse_atom(atom, Box_dim, varargin)`  
  Fuse all sites within a certain cutoff distance.

- **[heal_atom](heal_atom.html)**: `heal_atom(atom, Box_dim, ind, varargin)`  
  Heal sites in the atom struct by adding a certain atom type.

- **[ionize_atom](ionize_atom.html)**: `ionize_atom(type, resname, limits, nmax, varargin)`  
  Add ions within a certain region defined by limits.

- **[insert_atom](insert_atom.html)**: `insert_atom(atom, new_atom, position)`  
  Insert a new atom at the specified position.

- **[merge_atom](merge_atom.html)**: `merge_atom(atom1, Box1, atom2, type, Atom_label, r)`  
  Merge atom structs based on distance criteria.

- **[molid_rotate](molid_rotate.html)**: `molid_rotate(atom, Box_dim, MolID, rotate_dim)`  
  Rotate the atom struct based on MolID.

- **[molid_translate](molid_translate.html)**: `molid_translate(atom, trans_vec, MolID)`  
  Translate a specific molecule ID.

- **[noupdate_atom](noupdate_atom.html)**: `noupdate_atom(atom)`  
  Prevent updating of certain properties in the atom struct.

- **[occupancy_atom](occupancy_atom.html)**: `occupancy_atom(atom, Box_dim)`  
  Calculate occupancy of atoms within the box dimensions.

- **[overwrite_atom](overwrite_atom.html)**: `overwrite_atom(In_atom, atomtype, resname)`  
  Overwrite atom struct information with new data.

- **[place_atom](place_atom.html)**: `place_atom(atom, position)`  
  Place the atom struct at the specified position.

- **[position_molid](position_molid.html)**: `position_molid(atom, position_vec, MolID)`  
  Move a molecule ID to a certain position.

- **[protonate_atom](protonate_atom.html)**: `protonate_atom(atom, Box_dim, varargin)`  
  Protonate specified sites in the atom struct.

- **[remove_molid](remove_molid.html)**: `remove_molid(atom, MolID)`  
  Remove residue with a specific molecule ID.

- **[remove_occypancy_atom](remove_occypancy_atom.html)**: `remove_occypancy_atom(atom)`  
  Remove particles with identical coordinates to preceding ones.

- **[remove_residues](remove_residues.html)**: `remove_residues(atom, resnames, lo, hi, dim)`  
  Remove residues between specified limits in the simulation box.

- **[remove_resname](remove_resname.html)**: `remove_resname(atom, resnames)`  
  Remove residues with specified names.

- **[remove_SOL](remove_SOL.html)**: `remove_SOL(atom, atomname, lo, hi, dim)`  
  Remove solvent residues between specified limits.

- **[remove_type](remove_type.html)**: `remove_type(atom, typescell)`  
  Remove atom types specified in `typescell`.

- **[rename_atom](rename_atom.html)**: `rename_atom(atom, old_name, new_name)`  
  Rename atoms in the atom struct.

- **[rename_type](rename_type.html)**: `rename_type(atom, atomtype, new_atomtype, varargin)`  
  Rename atom types in the atom struct.

- **[replicate_atom](replicate_atom.html)**: `replicate_atom(atom, Box_dim, replicate)`  
  Replicate the atom struct along orthogonal dimensions.

- **[replace_atom](replace_atom.html)**: `replace_atom(new_atom, prev_atom, molid_index)`  
  Replace molecule ID in an atom struct with a new atom struct.

- **[resname_atom](resname_atom.html)**: `resname_atom(atom)`  
  Guess residue names for all atom types.

- **[slice_atom](slice_atom.html)**: `slice_atom(atom, limits, invert)`  
  Slice the atom struct within specified limits.

- **[slice_box](slice_box.html)**: `slice_box(atom, Box_dim, limits)`  
  Slice a simulation box within given limits.

- **[slice_molid](slice_molid.html)**: `slice_molid(atom, limits, invert)`  
  Slice molecules within specified limits.

- **[slice_triclinic_atom](slice_triclinic_atom.html)**: `slice_triclinic_atom(atom, limits, invert)`  
  Slice a triclinic atom struct within limits.

- **[solvate_atom](solvate_atom.html)**: `solvate_atom(limits, density, r, maxsol, solute_atom, varargin)`  
  Generate a solvent structure within specified limits.

- **[spc2tip4p](spc2tip4p.html)**: `spc2tip4p(atom)`  
  Convert SPC water molecules to TIP4P model.

- **[spc2tip5p](spc2tip5p.html)**: `spc2tip5p(atom)`  
  Convert SPC water molecules to TIP5P model.

- **[spce2tip4p](spce2tip4p.html)**: `spce2tip4p(atom)`  
  Convert SPC/E water molecules to TIP4P model.

- **[sphere_atom](sphere_atom.html)**: `sphere_atom(atom, Box_dim, center, radius)`  
  Create a spherical region of atoms.

- **[substitute_atom](substitute_atom.html)**: `substitute_atom(atom, Box_dim, NumOctSubst, O1, O2, minO2O2_dist, varargin)`  
  Perform isomorphous substitution of atoms.

- **[substitute_NonCentroSymm_atom](substitute_NonCentroSymm_atom.html)**: `substitute_NonCentroSymm_atom(atom, Box_dim, replace_type, varargin)`  
  Substitute non-centrosymmetric atoms.

- **[tile_atom](tile_atom.html)**: `tile_atom(atom, scale_vec, Box_dim, Resname)`  
  Tile the atom struct in a specific direction.

- **[tip3p2tip4p](tip3p2tip4p.html)**: `tip3p2tip4p(atom)`  
  Convert TIP3P water molecules to TIP4P model.

- **[translate_atom](translate_atom.html)**: `translate_atom(atom, trans_vec, Resname)`  
  Translate a residue by a specified vector.

- **[translate_molid](translate_molid.html)**: `translate_molid(atom, trans_vec, molid)`  
  Translate a molecule ID by a specified vector.

- **[tube_atom](tube_atom.html)**: `tube_atom(atom, scale_vec, Box_dim, Resname)`  
  Create a nanotube structure from the atom struct.

- **[update_atom](update_atom.html)**: `update_atom(atom)`  
  Update molecule and atom indices in the atom struct.

---

## Add/Create/Replicate/Overwrite Atoms

*(This section overlaps with "Specific Atom Struct Functions" and "Add/Create/Replicate/Overwrite Atoms" categories. For brevity, refer to the "Specific Atom Struct Functions" section above.)*

---

## Slice Out a Region of the Box

*(This section overlaps with "Slice Out a Region of the Box" under "Specific Atom Struct Functions". Refer to the relevant functions listed above.)*

---

## Translate or Rotate Functions

- **[bend_atom](bend_atom.html)**: `bend_atom(atom, Box_dim, Radii)`  
  Bend an atom struct.

- **[center_atom](center_atom.html)**: `center_atom(atom, Box_dim, resname, dim)`  
  Center the atom with respect to the `resname` molecule.

- **[condense_atom](condense_atom.html)**: `condense_atom(atom, Box_dim, s)`  
  Minimize the box size and remove gaps between `molids`.

- **[molid_rotate](molid_rotate.html)**: `molid_rotate(atom, Box_dim, MolID, rotate_dim)`  
  Rotate the atom struct based on MolID.

- **[molid_translate](molid_translate.html)**: `molid_translate(atom, trans_vec, MolID)`  
  Translate a specific molecule ID.

- **[place_atom](place_atom.html)**: `place_atom(atom, position)`  
  Place the atom struct at the specified position.

- **[position_molid](position_molid.html)**: `position_molid(atom, position_vec, MolID)`  
  Move a molecule ID to a certain position.

- **[rotate_atom](rotate_atom.html)**: `rotate_atom(atom, rotation_matrix)`  
  Rotate the atom struct by a specified matrix.

- **[translate_atom](translate_atom.html)**: `translate_atom(atom, trans_vec)`  
  Translate the atom struct by a specified vector.

- **[translate_molid](translate_molid.html)**: `translate_molid(atom, trans_vec, molid)`  
  Translate the molecule ID by a specified vector.

---

## Make Triclinic/Orthogonal Box

- **[frac2atom](frac2atom.html)**: `frac2atom(atom, Box_dim, angleparam, angletype)`  
  Transform fractional coordinates to Cartesian coordinates.

- **[orto_atom](orto_atom.html)**: `orto_atom(atom, Box_dim)`  
  Transform a triclinic atom struct to an orthogonal one.

- **[triclinic_atom](triclinic_atom.html)**: `triclinic_atom(atom, Box_dim, angleparam, angletype)`  
  Transform an orthogonal atom struct to a triclinic one.

---

## Wrap/Unwrap Functions

- **[unwrap_atom](unwrap_atom.html)**: `unwrap_atom(atom, Box_dim, dim)`  
  Unwrap the atom struct along the specified dimension.

---

## Keep/Remove Functions

*(This section overlaps with "Specific Atom Struct Functions" and "Keep/Remove Functions" categories. Refer to the relevant functions listed above.)*

---

## List of Available Solvents

### Water

Use the following `.pdb` or `.gro` files with the `solvate_atom` function. You can also use custom solvent boxes to solvate a simulation cell.

- `864_spc.gro` | `.pdb` – Equilibrated SPC water box
- `864_spce.gro` | `.pdb` – Equilibrated SPC/E water box
- `864_tip3p.gro` | `.pdb` – Equilibrated TIP3P water box
- `864_tip4p.gro` | `.pdb` – Equilibrated TIP4P water box
- `864_tip5p.gro` | `.pdb` – Equilibrated TIP5P water box
- `96spc_hex_ice_h.gro` | `.pdb` – Equilibrated SPC hex-ice water box
- `96tip4p_hex_ice_h.gro` | `.pdb` – Equilibrated TIP4P hex-ice water box
- `864_swm4_ndp.gro` | `.pdb` – Polarizable water v1
- `864_swm4_ndp_vds.gro` | `.pdb` – Polarizable water v2

### Organics

- `500xEtOH.gro` | `.pdb` – Equilibrated SPC water box

### Mineral

- `1xPyro_Lee_Guggenheim_1981_alfabeta90.pdb`
- `Pyrophyllite.pdb`
- Hexagonal layered particles:
  - Pyrophyllite / Montmorillonite
  - Talc / Laponite

---

## Conversion Functions

- **[spc2tip4p](spc2tip4p.html)**: `spc2tip4p(filename)`  
  Convert a `.gro` or `.pdb` file with SPC water to TIP4P water.

- **[spc2tip5p](spc2tip5p.html)**: `spc2tip5p(filename)`  
  Convert a `.gro` or `.pdb` file with SPC water to TIP5P water.

- **[spce2tip4p](spce2tip4p.html)**: `spce2tip4p(filename)`  
  Convert a `.gro` or `.pdb` file with SPC/E water to TIP4P water.

- **[tip3p2tip4p](tip3p2tip4p.html)**: `tip3p2tip4p(filename)`  
  Convert a `.gro` file with TIP3P water to TIP4P water.

---

## Custom Topology Tools

*(No functions listed under this category. Add relevant functions if available.)*

---

## Miscellaneous FF Functions

- **[print_top](print_top.html)**: `print_top(atom, Box_dim, varargin)`  
  Print or generate topology-related data.

- **[import_ff_table](import_ff_table.html)**: `import_ff_table(filename, varargin)`  
  Import force field parameter tables.

- **[change_top](change_top.html)**: `change_top(atom, Box_dim, varargin)`  
  Modify the topology file or its parameters.

- **[mass_atom_clayff](mass_atom_clayff.html)**: `mass_atom_clayff(atom)`  
  Define or calculate masses for Clayff atoms.

- **[smear_charge](smear_charge.html)**: `smear_charge(atom, Box_dim, varargin)`  
  Distribute charge across atoms, possibly using charge smearing techniques.

---

## Force Fields

### MINFF, with Atom Types by MHolmboe

- **[minff_atom](minff_atom.html)**: `minff_atom(atom, Box_dim, varargin)`  
  Assign MINFF atom types with edge healing.

- **[charge_minff_atom](charge_minff_atom.html)**: `charge_minff_atom(atom, Box_dim, varargin)`  
  Set the charge for MINFF atom types.

### CLAYFF, with Atom Types by MHolmboe

- **[charge_atom](charge_atom.html)**: `charge_atom(atom, Box_dim, ffname, watermodel, varargin)`  
  Charge the atom according to CLAYFF or INTERFACE force fields.

- **[charge_clayff_2004_atom](charge_clayff_2004_atom.html)**: `charge_clayff_2004_atom(atom, Box_dim, varargin)`  
  Set the charge for the original CLAYFF atom types from the Cygan et al., 2004 paper.

- **[charge_clayff_atom](charge_clayff_atom.html)**: `charge_clayff_atom(atom, Box_dim, varargin)`  
  Set the charge for CLAYFF atom types.

- **[charge_opls_go_atom](charge_opls_go_atom.html)**: `charge_opls_go_atom(atom, Box_dim, varargin)`  
  Set the charge for specific OPLS atom types.

- **[check_clayff_2004_charge](check_clayff_2004_charge.html)**: `check_clayff_2004_charge(atom)`  
  Check the charge of the original CLAYFF atom types.

- **[check_clayff_charge](check_clayff_charge.html)**: `check_clayff_charge(atom)`  
  Check the charge of the CLAYFF atom types.

- **[check_clayff_H2Odens](check_clayff_H2Odens.html)**: `check_clayff_H2Odens(atom, Box_dim)`  
  Check the approximate water density for a CLAYFF system.

- **[check_H2Odens](check_H2Odens.html)**: `check_H2Odens(atom, Box_dim)`  
  Compute the water density.

- **[clayff_2004_atom](clayff_2004_atom.html)**: `clayff_2004_atom(atom, Box_dim, varargin)`  
  Assign the original CLAYFF atom types by Cygan et al., 2004, with edge healing.

- **[clayff_2004_param](clayff_2004_param.html)**: `clayff_2004_param(Atom_label, varargin)`  
  Hold ion and CLAYFF atom type parameters for the original CLAYFF force field.

- **[clayff_atom](clayff_atom.html)**: `clayff_atom(atom, Box_dim, varargin)`  
  Assign CLAYFF atom types with edge healing.

- **[clayff_param](clayff_param.html)**: `clayff_param(Atom_label, varargin)`  
  Hold ion and CLAYFF atom type parameters.

- **[clayff210_atom](clayff210_atom.html)**: `clayff210_atom(atom, Box_dim, varargin)`  
  Assign modified CLAYFF atom types, with edge healing.

- **[clayff211_atom](clayff211_atom.html)**: `clayff211_atom(atom, Box_dim, varargin)`  
  Assign modified CLAYFF atom types (faster version), with edge healing.

- **[tweak_charge_atom](tweak_charge_atom.html)**: `tweak_charge_atom(atom)`  
  Tweak the charge of the atom struct to correct rounding errors.

### INTERFACE from Heinz 2005, 2013, with Atom Types by MHolmboe

- **[charge_interface_atom](charge_interface_atom.html)**: `charge_interface_atom(atom, Box_dim, varargin)`  
  Set the charge for INTERFACE atom types.

- **[charge_interface15_atom](charge_interface15_atom.html)**: `charge_interface15_atom(atom, Box_dim, varargin)`  
  Set the charge for INTERFACE 1.5 atom types.

- **[check_interface_charge](check_interface_charge.html)**: `check_interface_charge(atom)`  
  Check the charge of the INTERFACE atom types.

- **[check_interface15_charge](check_interface15_charge.html)**: `check_interface15_charge(atom)`  
  Check the charge of the INTERFACE 1.5 atom types.

- **[interface_atom](interface_atom.html)**: `interface_atom(atom, Box_dim, varargin)`  
  Assign atoms according to the INTERFACE atom types, with modifications for edges.

- **[interface_param](interface_param.html)**: `interface_param(Atom_label, water_model)`  
  Hold extended INTERFACE force field parameters.

- **[interface15_atom](interface15_atom.html)**: `interface15_atom(atom, Box_dim, varargin)`  
  Assign atoms according to the INTERFACE 1.5 atom types, with modifications for edges.

- **[interface15_param](interface15_param.html)**: `interface15_param(Atom_label, water_model)`  
  Hold extended INTERFACE 1.5 force field parameters.

- **[check_interface_charge](check_interface_charge.html)**: *(Duplicate)*

- **[check_interface15_charge](check_interface15_charge.html)**: *(Duplicate)*

- **[interface15_silica_atom](interface15_silica_atom.html)**: `interface15_silica_atom(atom, Box_dim, varargin)`  
  Assign atom types for the INTERFACE 1.5 force field, specific to silica.

### Graphene Oxide Modeled with OPLS/AA

- **[opls_go_atom](opls_go_atom.html)**: `opls_go_atom(atom, Box_dim, rmin, rlarge)`  
  Smear out the charge around -OH and epoxide groups in graphene oxide.

- **[oplsaa_go_param](oplsaa_go_param.html)**: `oplsaa_go_param(Atom_label, water_model)`  
  Hold the extended OPLS-AA force field parameters for graphite oxide.

- **[charge_opls_go_atom](charge_opls_go_atom.html)**: `charge_opls_go_atom(atom, Box_dim, varargin)`  
  Set the charge for specific OPLS atom types.

---

## Writing Topology Files

- **[write_atom_itp](write_atom_itp.html)**: `write_atom_itp(atom, Box_dim, filename_out, varargin)`  
  Create and print a GROMACS `.itp` file for CLAYFF or INTERFACE force fields.

- **[write_atom_lmp](write_atom_lmp.html)**: `write_atom_lmp(atom, Box_dim, filename_out, varargin)`  
  Create and print a LAMMPS data file (`.lj`) for CLAYFF systems.

- **[write_atom_oplsaa_go_itp](write_atom_oplsaa_go_itp.html)**: `write_atom_oplsaa_go_itp(atom, Box_dim, filename_out, varargin)`  
  Create and print a GROMACS `.itp` file for OPLS-AA or GO systems.

---

## Bonded and Nonbonded Parameters

- **[bonded_parameters](bonded_parameters.html)**: `bonded_parameters(atom, varargin)`  
  Define bonded parameters for atoms.

- **[nonbonded_parameters](nonbonded_parameters.html)**: `nonbonded_parameters(atom, varargin)`  
  Define nonbonded parameters for atoms.

- **[nonbonded_ff](nonbonded_ff.html)**: `nonbonded_ff(atom, varargin)`  
  Define nonbonded force field parameters.

---

## Lennard-Jones and Coulomb Potentials

- **[buckinghamcoul](buckinghamcoul.html)**: `buckinghamcoul(atom, Box_dim, varargin)`  
  Calculate interactions using the Buckingham potential and Coulombic forces.

- **[ljcoul_12_6](ljcoul_12_6.html)**: `ljcoul_12_6(atom, Box_dim, varargin)`  
  Handle Lennard-Jones 12-6 potential along with Coulombic interactions.

- **[ljcoul_C12C6C4](ljcoul_C12C6C4.html)**: `ljcoul_C12C6C4(atom, Box_dim, varargin)`  
  Lennard-Jones and Coulomb potential with C12, C6, C4 terms.

- **[ljcoul_C12C6](ljcoul_C12C6.html)**: `ljcoul_C12C6(atom, Box_dim, varargin)`  
  Lennard-Jones and Coulomb potential with C12, C6 terms.

- **[ljcoul](ljcoul.html)**: `ljcoul(atom, Box_dim, varargin)`  
  General Lennard-Jones and Coulomb interaction function.

- **[ljcoul_2x2x](ljcoul_2x2x.html)**: `ljcoul_2x2x(atom, Box_dim, varargin)`  
  Lennard-Jones and Coulomb potential with 2x factors.

- **[ljcoul_2x](ljcoul_2x.html)**: `ljcoul_2x(atom, Box_dim, varargin)`  
  Variation of Lennard-Jones Coulomb with 2x factors.

---

## Objective Functions and Force Calculations

- **[buckinghamcoul_objective_func](buckinghamcoul_objective_func.html)**: `buckinghamcoul_objective_func(atom, Box_dim, varargin)`  
  Objective function for fitting Buckingham and Coulomb potentials.

- **[ljcoul_objective_func](ljcoul_objective_func.html)**: `ljcoul_objective_func(atom, Box_dim, varargin)`  
  Objective function for Lennard-Jones and Coulomb potentials.

- **[ljcoul_2x_objective_func](ljcoul_2x_objective_func.html)**: `ljcoul_2x_objective_func(atom, Box_dim, varargin)`  
  Objective function for 2x Lennard-Jones and Coulomb potentials.

- **[ljcoul_force](ljcoul_force.html)**: `ljcoul_force(atom, Box_dim, varargin)`  
  Compute forces based on Lennard-Jones and Coulomb potentials.

- **[ljcoul_2x_force](ljcoul_2x_force.html)**: `ljcoul_2x_force(atom, Box_dim, varargin)`  
  Force calculation involving 2x Lennard-Jones and Coulomb potentials.

- **[ljcoul_force_C12C6C4](ljcoul_force_C12C6C4.html)**: `ljcoul_force_C12C6C4(atom, Box_dim, varargin)`  
  Force calculation for Lennard-Jones potential with C12, C6, C4 terms.

- **[ljcoul_force_objective_func](ljcoul_force_objective_func.html)**: `ljcoul_force_objective_func(atom, Box_dim, varargin)`  
  Objective function for force calculations with Lennard-Jones and Coulomb potentials.

---

## Automated Fitting Tools

- **[autofit_C6C8C10xljcoul](autofit_C6C8C10xljcoul.html)**: `autofit_C6C8C10xljcoul(atom, Box_dim, varargin)`  
  Automated fitting for Lennard-Jones Coulombic parameters with C6, C8, C10 terms.

- **[autofit_2xljcoul_func](autofit_2xljcoul_func.html)**: `autofit_2xljcoul_func(atom, Box_dim, varargin)`  
  Function for automatic fitting of 2x Lennard-Jones and Coulomb parameters.

- **[autofit_2xljcoul_batch](autofit_2xljcoul_batch.html)**: `autofit_2xljcoul_batch(atom, Box_dim, varargin)`  
  Batch fitting for 2x Lennard-Jones and Coulomb potentials.

- **[autofit_2xljcoul](autofit_2xljcoul.html)**: `autofit_2xljcoul(atom, Box_dim, varargin)`  
  Automated fitting of 2x Lennard-Jones and Coulomb parameters.

- **[autofit_geometric2LB](autofit_geometric2LB.html)**: `autofit_geometric2LB(atom, Box_dim, varargin)`  
  Automated geometric fitting for Lennard-Jones potential with Born-Mayer interactions.

- **[autofit_buckcoul](autofit_buckcoul.html)**: `autofit_buckcoul(atom, Box_dim, varargin)`  
  Automated fitting for Buckingham and Coulomb potentials.

- **[autofit_force_2xljcoul](autofit_force_2xljcoul.html)**: `autofit_force_2xljcoul(atom, Box_dim, varargin)`  
  Automated fitting of force parameters for 2x Lennard-Jones and Coulomb interactions.

- **[autofit_ljcoul](autofit_ljcoul.html)**: `autofit_ljcoul(atom, Box_dim, varargin)`  
  Automated fitting for Lennard-Jones and Coulombic interactions.

---

## Additional Resources

- **Documentation**: Each function is documented in its respective HTML file. Refer to the function links above for detailed usage instructions and examples.

- **Examples**: Check the `examples` directory for sample scripts demonstrating how to use the various functions within the library.

- **Contributing**: Contributions are welcome! Please fork the repository and submit a pull request with your enhancements.

- **License**: This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Acknowledgements

Special thanks to Hendrik Heinz for the INTERFACE force field distribution and to M. Holmboe for developing the MINFF and CLAYFF atom types.

---

*Happy Computing!*