# Atom MATLAB Library

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![MATLAB Version](https://img.shields.io/badge/MATLAB-2024-orange.svg)

## Overview

The **Atom MATLAB Library** is a comprehensive collection of MATLAB functions designed for handling and analyzing atomic structures, molecular simulations, and related computational tasks. Whether you're importing data from various file formats, performing structural manipulations, calculating distances and bonds, or exporting data for simulation tools like GROMACS and LAMMPS, this quirky MATLAB library/toolbox may (not) streamline your workflow.

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

- **[import_xvg](import_functions/import_xvg.m)**: `import_xvg(filename)`  
  Import a GROMACS `.xvg` file.

- **[import_gmx_energy](import_functions/import_gmx_energy.m)**: `import_gmx_energy(filename)`  
  Import a GROMACS energy file.

- **[import_bar](import_functions/import_bar.m)**: `import_bar(filename)`  
  Import a BAR file.

- **[import_dat](import_functions/import_dat.m)**: `import_dat(filename)`  
  Import a `.dat` file.

- **[import_red_charges](import_functions/import_red_charges.m)**: `import_red_charges(filename)`  
  Import reduced charges from a file.

- **[import_ddec_charges](import_functions/import_ddec_charges.m)**: `import_ddec_charges(filename)`  
  Import DDEC charges from a file.

### Import of Structure Files

- **[import_atom_car](import_functions/import_atom_car.m)**: `import_atom_car(filename, varargin)`  
  Import `.car` files from Hendrik Heinz INTERFACE force field distribution and write out GROMACS `.itp` and `.pdb` files.

- **[import_atom_gro](import_functions/import_atom_gro.m)**: `import_atom_gro(filename)`  
  Import `.gro` files into the atom struct.

- **[import_atom_gro_fscanf](import_functions/import_atom_gro_fscanf.m)**: `import_atom_gro_fscanf(filename)`  
  Import `.gro` files using `fscanf` (alternative method).

- **[import_atom_gro_octave](import_functions/import_atom_gro_octave.m)**: `import_atom_gro_octave(filename)`  
  Import `.gro` files using Octave-compatible code.

- **[import_atom_mol2](import_functions/import_atom_mol2.m)**: `import_atom_mol2(filename)`  
  Import `.mol2` files into the atom struct.

- **[import_atom_pdb](import_functions/import_atom_pdb.m)**: `import_atom_pdb(filename)`  
  Import `.pdb` files into the atom struct.

- **[import_atom_pqr](import_functions/import_atom_pqr.m)**: `import_atom_pqr(filename)`  
  Import `.pqr` files into the atom struct.

- **[import_atom_xyz](import_functions/import_atom_xyz.m)**: `import_atom_xyz(filename)`  
  Import an `.xyz` file into the atom struct.

- **[import_atom_poscar](import_functions/import_atom_poscar.m)**: `import_atom_poscar(filename)`  
  Import a VASP POSCAR file into the atom struct.

- **[import_atom](import_functions/import_atom.m)**: `import_atom(filename)`  
  Import a `.xyz`, `.gro`, or `.pdb` file into a structure variable called `atom`.

- **[import_xyz](import_functions/import_xyz.m)**: `import_xyz(filename)`  
  Import an `.xyz` file. Atom types should be letters, not numbers. Use `import_atom_xyz` instead.

### Import of Trajectory Files

- **[import_gro_traj](import_functions/import_gro_traj.m)**: `import_gro_traj(filename, varargin)`  
  Import a structure and a `.gro` trajectory file.

- **[import_xyz_traj](import_functions/import_xyz_traj.m)**: `import_xyz_traj(filenametraj)`  
  Import a structure and an `.xyz` trajectory file.

- **[import_mc_pdb_traj](import_functions/import_mc_pdb_traj.m)**: `import_mc_pdb_traj(filename, varargin)`  
  Import a structure and a `.pdb` trajectory file, handling changing numbers of particles.

- **[import_pdb_traj](import_functions/import_pdb_traj.m)**: `import_pdb_traj(filename, varargin)`  
  Import a structure and a `.pdb` trajectory file.

- **[import_traj](import_functions/import_traj.m)**: `import_traj(filenameconf, filenametraj)`  
  Import a structure and a `.dcd`, `.trr`, `.xtc`, `.xyz`, or `.gro` trajectory file.

- **[import_trr](import_functions/import_trr.m)**: `import_trr(filenameconf, filenametraj)`  
  Import a structure and a `.trr` trajectory file.

- **[import_trrv2](import_functions/import_trrv2.m)**: `import_trrv2(filenameconf, filenametraj)`  
  Import a structure and a `.trr` trajectory file (version 2).

- **[import_xtc](import_functions/import_xtc.m)**: `import_xtc(filenameconf, filenamextc)`  
  Import a structure and an `.xtc` file.

- **[import_xtcv2](import_functions/import_xtcv2.m)**: `import_xtcv2(filenameconf, filenamextc)`  
  Import a structure and an `.xtc` file (version 2).

### Miscellaneous Import Functions

- **[import_ave_gro](import_functions/import_ave_gro.m)**: `import_ave_gro(filename)`  
  Import an averaged structure from a `.gro` trajectory.

- **[import_mclf](import_functions/import_mclf.m)**: `import_mclf(filename)`  
  Import a multi-configurational London force (MCLF) file.

- **[import_mclf_C6](import_functions/import_mclf_C6.m)**: `import_mclf_C6(filename)`  
  Import MCLF C6 dispersion parameters.

- **[import_mclf_C8](import_functions/import_mclf_C8.m)**: `import_mclf_C8(filename)`  
  Import MCLF C8 dispersion parameters.

- **[import_mclf_C10](import_functions/import_mclf_C10.m)**: `import_mclf_C10(filename)`  
  Import MCLF C10 dispersion parameters.

- **[import_mclf_dispersion](import_functions/import_mclf_dispersion.m)**: `import_mclf_dispersion(filename)`  
  Import dispersion parameters for MCLF.

- **[import_cp2k](import_functions/import_cp2k.m)**: `import_cp2k(filename)`  
  Import a CP2K output file.

- **[import_cp2k_resp](import_functions/import_cp2k_resp.m)**: `import_cp2k_resp(filename)`  
  Import RESP charges from CP2K.

---

## Export Functions

### Export of Structure/Topology Files

- **[write_atom_all](export_functions/write_atom_all.m)**: `write_atom_all(atom, Box_dim, filename_out, varargin)`  
  Write various file types for the atom struct, best suited for Clayff systems.

- **[write_atom_cif](export_functions/write_atom_cif.m)**: `write_atom_cif(atom, Box_dim, filename_out)`  
  Write a basic `.cif` file from the atom struct.

- **[write_atom_gro](export_functions/write_atom_gro.m)**: `write_atom_gro(atom, Box_dim, filename_out)`  
  Write a `.gro` file from the atom struct, optionally including velocities.

- **[write_atom_itp](export_functions/write_atom_itp.m)**: `write_atom_itp(atom, Box_dim, filename_out, varargin)`  
  Create and print a GROMACS `.itp` file, best for MINFF, CLAYFF (INTERFACE maybe) force fields.

- **[write_atom_lmp](export_functions/write_atom_lmp.m)**: `write_atom_lmp(atom, Box_dim, filename_out, varargin)`  
  Create and print a LAMMPS data file (`.lj`), suited for Clayff systems.

- **[write_atom_mol2](export_functions/write_atom_mol2.m)**: `write_atom_mol2(atom, Bond_index, Box_dim, filename_out)`  
  Write a `.mol2` file from the atom struct.

- **[write_atom_oplsaa_go_itp](export_functions/write_atom_oplsaa_go_itp.m)**: `write_atom_oplsaa_go_itp(atom, Box_dim, filename_out, varargin)`  
  Create and print a GROMACS `.itp` file for OPLS-AA or GO systems.

- **[write_atom_pdb](export_functions/write_atom_pdb.m)**: `write_atom_pdb(atom, Box_dim, filename_out)`  
  Write a `.pdb` file from the atom struct using GROMACS.

- **[write_atom_pqr](export_functions/write_atom_pqr.m)**: `write_atom_pqr(atom, Box_dim, filename_out, varargin)`  
  Write a `.pqr` file from the atom struct.

- **[write_atom_psf](export_functions/write_atom_psf.m)**: `write_atom_psf(atom, Box_dim, filename_out, varargin)`  
  Write a `.psf` file from the atom struct.

- **[write_atom_xyz](export_functions/write_atom_xyz.m)**: `write_atom_xyz(atom, Box_dim, filename_out)`  
  Write an `.xyz` file from the atom struct.

- **[write_atom](export_functions/write_atom.m)**: `write_atom(atom, Box_dim, filename_out, varargin)`  
  Write different file types (`.gro`, `.pdb`, `.xyz`, `.itp`, etc.) based on filename and parameters.

- **[write_atom_sdf](export_functions/write_atom_sdf.m)**: `write_atom_sdf(atom, Box_dim, filename_out)`  
  Write an `.sdf` file from the atom struct.

- **[write_atom_dodecahedron_gro](export_functions/write_atom_dodecahedron_gro.m)**: `write_atom_dodecahedron_gro(atom, Box_dim, filename_out)`  
  Write a `.gro` file using a dodecahedron-shaped simulation box.

- **[write_atom_top](export_functions/write_atom_top.m)**: `write_atom_top(atom, Box_dim, filename_out)`  
  Write a topology file (`.top`) from the atom struct.

- **[CONECT_atom](export_functions/CONECT_atom.m)**: `CONECT_atom(atom, Box_dim, filename_out)`  
  Write CONECT records for a PDB file.

### Export of Trajectory Files

- **[write_gro_traj](export_functions/write_gro_traj.m)**: `write_gro_traj(atom, traj, Box_dim, filename_out)`  
  Write a `.gro` trajectory file.

- **[write_atom_multiple_gro](export_functions/write_atom_multiple_gro.m)**: `write_atom_multiple_gro(atom, traj, filename_out)`  
  Write multiple `.gro` files for trajectory output.

- **[write_pdb_traj](export_functions/write_pdb_traj.m)**: `write_pdb_traj(atom, traj, Box_dim, filename_out)`  
  Write a `.pdb` trajectory file.

- **[write_xyz_traj](export_functions/write_xyz_traj.m)**: `write_xyz_traj(atom, traj, Box_dim, filename_out)`  
  Write a `.xyz` trajectory file.

- **[write_ave_gro](export_functions/write_ave_gro.m)**: `write_ave_gro(atom, traj, Box_dim, filename_out)`  
  Write an average structure from a `.gro` trajectory.

- **[write_ave_pdb](export_functions/write_ave_pdb.m)**: `write_ave_pdb(atom, traj, Box_dim, filename_out)`  
  Write an average structure from a `.pdb` trajectory.

### Export of Other Formats and Data

- **[export_ndx](export_functions/export_ndx.m)**: `export_ndx(atom, Box_dim, filename_out)`  
  Export an index file (`.ndx`) from the atom struct.

- **[write_xvg](export_functions/write_xvg.m)**: `write_xvg(filename, data)`  
  Export data in GROMACS `.xvg` format.

- **[write_tabulated_potentials](export_functions/write_tabulated_potentials.m)**: `write_tabulated_potentials(filename, data)`  
  Write tabulated potential files.

- **[write_ff](export_functions/write_ff.m)**: `write_ff(atom, filename_out)`  
  Write force field parameters.

- **[write_ffnonbonded](export_functions/write_ffnonbonded.m)**: `write_ffnonbonded(atom, filename_out)`  
  Write non-bonded force field parameters.

- **[write_ffnonbonded_C6C12](export_functions/write_ffnonbonded_C6C12.m)**: `write_ffnonbonded_C6C12(atom, filename_out)`  
  Write non-bonded parameters (C6, C12) for a force field.

### Miscellaneous Export Functions

- **[replace_string](export_functions/replace_string.m)**: `replace_string(filename_in, filename_out, old_string, new_string)`  
  Replace strings in files.

---

## Neighbor/Distance Functions

- **[bond_atom](neigh_functions/bond_atom.m)**: `bond_atom(atom, Box_dim, max_long_dist)`  
  Assign all bonds to a `Bond_matrix` and `Bond_index` variable.

- **[bond_angle_atom](neigh_functions/bond_angle_atom.m)**: `bond_angle_atom(atom, Box_dim, varargin)`  
  Find all bonds and angles of the atom struct. `'More'` is an optional argument.

- **[bond_angle_dihedral_atom](neigh_functions/bond_angle_dihedral_atom.m)**: `bond_angle_dihedral_atom(atom, Box_dim, varargin)`  
  Find all bonds, angles, and dihedrals of the atom struct. `Rmaxshort` and `Rmaxlong` as well as `'more'` are optional arguments.

- **[bond_angle_type](neigh_functions/bond_angle_type.m)**: `bond_angle_type(atom1, atom2, Box_dim, rmin, rmax, angle_limit, varargin)`  
  Find all bonds and angles of the atom types.

- **[bond_valence_atom](general_functions/bond_valence_atom.m)**: `bond_valence_atom(atom, Box_dim, varargin)`  
  Calculate bond valence values using the bond valence method.

- **[bond_valence_data](general_functions/bond_valence_data.m)**: `bond_valence_data(ion1, ion2, R, varargin)`  
  Fetch data to calculate bond valence values for specified atom types.

- **[cell_list_dist_matrix_atom](cell_list_dist_matrix_atom.html)**: `cell_list_dist_matrix_atom(atom, Box_dim, varargin)`  
  Calculate the distance matrix from the atom struct using a cell list algorithm.

- **[closest_atom](build_functions/closest_atom.m)**: `closest_atom(atom1, atom2, Box_dim)`  
  Return the `atom1` struct with the `nMolId's` in `atom1` closest to the `atom2` struct.

- **[dist_matrix_atom](neigh_functions/dist_matrix_atom.m)**: `dist_matrix_atom(atom, Box_dim)`  
  Calculate the distance matrix from the atom struct.

- **[dist_matrix_noPBC_atom](neigh_functions/dist_matrix_noPBC_atom.m)**: `dist_matrix_noPBC_atom(atom, Box_dim)`  
  Calculate the distance matrix without applying periodic boundary conditions.

- **[dist_matrix_xyz](neigh_functions/dist_matrix_xyz.m)**: `dist_matrix_xyz(XYZ, Box_dim)`  
  Calculate the distance matrix from XYZ coordinates.

- **[find_bonded_atom](neigh_functions/find_bonded_atom.m)**: `find_bonded_atom(atom, bond_matrix, label1, label2)`  
  Perform a cross-check of the bond matrix.

- **[find_pair_atom](neigh_functions/find_pair_atom.m)**: `find_pair_atom(atom, bond_matrix, pair1, pair2)`  
  Find and return specific atom pairs from the bond matrix.

- **[list_bonds_atom](neigh_functions/list_bonds_atom.m)**: `list_bonds_atom(atom, bond_matrix)`  
  List all bonds in the atom struct based on the bond matrix.

- **[neigh_atom](neigh_functions/neigh_atom.m)**: `neigh_atom(atom, Box_dim, rmax, varargin)`  
  Check which neighbors each atom has and output their information.

- **[neighbor_func](neigh_functions/neighbor_func.m)**: `neighbor_func(solute_index, XYZ_solute, XYZ_data, Box_dim, radius)`  
  Scan XYZ data and check which atoms are within a certain radius. Outputs the neighbor index.

- **[neighbor_atom](neigh_functions/neighbor_atom.m)**: `neighbor_atom(atom, Box_dim, radius)`  
  Check the neighbors of each atom and return their indices.

- **[rdf_atom](general_functions/rdf_atom.m)**: `rdf_atom(atom, Box_dim, varargin)`  
  Calculate the radial distribution function and the coordination number, with optional Gaussian smoothing.

- **[cn_atom](neigh_functions/cn_atom.m)**: `cn_atom(atom, Box_dim, rmax)`  
  Calculate the coordination number of atoms within a specified radius.

- **[recalc_bond_atom](neigh_functions/recalc_bond_atom.m)**: `recalc_bond_atom(atom, bond_matrix, varargin)`  
  Recalculate bonds for the atom struct.

- **[remove_H2O](neigh_functions/remove_H2O.m)**: `remove_H2O(atom)`  
  Remove water molecules (`H2O`) from the atom struct.

- **[remove_sametype_bond](neigh_functions/remove_sametype_bond.m)**: `remove_sametype_bond(atom, bond_matrix)`  
  Remove bonds between atoms of the same type.

- **[find_H2O](neigh_functions/find_H2O.m)**: `find_H2O(atom)`  
  Identify and return water molecules (`H2O`) within the atom struct.

- **[bond_matrix_atom](neigh_functions/bond_matrix_atom.m)**: `bond_matrix_atom(atom, Box_dim)`  
  Generate a bond matrix for the atom struct.

---

## Specific Atom Struct Functions

- **[add2atom](build_functions/add2atom.m)**: `add2atom(XYZ_labels, XYZ_data, varargin)`  
  Append XYZ atom type labels and XYZ data to an existing atom struct.

- **[adjust_H_atom](build_functions/adjust_H_atom.m)**: `adjust_H_atom(atom, Box_dim)`  
  Adjust hydrogen atoms in the atom struct.

- **[cat_atom](build_functions/cat_atom.m)**: `cat_atom(atom_1, atom_2)`  
  Concatenate two atom structs.

- **[closest_atom](build_functions/closest_atom.m)**: `closest_atom(atom, Box_dim, ref_atom)`  
  Find the closest atom to the reference atom.

- **[copy_atom](build_functions/copy_atom.m)**: `copy_atom(atom, atomtype, new_atomtype, new_resname, trans_vec, varargin)`  
  Copy and translate atoms in the atom struct.

- **[create_atom](build_functions/create_atom.m)**: `create_atom(type, resname, limits, nmax, varargin)`  
  Create new atoms, useful for adding ions to a system.

- **[create_grid_atom](build_functions/create_grid_atom.m)**: `create_grid_atom(atom_label, nM, limits, dim, varargin)`  
  Put ions on a grid plane and add them to an atom struct.

- **[duplicate_atom](build_functions/duplicate_atom.m)**: `duplicate_atom(atom, molID)`  
  Duplicate residue with `molid` MolID.

- **[fuse_atom](build_functions/fuse_atom.m)**: `fuse_atom(atom, Box_dim, varargin)`  
  Fuse all sites within a certain cutoff distance.

- **[heal_atom](build_functions/heal_atom.m)**: `heal_atom(atom, Box_dim, ind, varargin)`  
  Heal sites in the atom struct by adding a certain atom type.

- **[ionize_atom](build_functions/ionize_atom.m)**: `ionize_atom(type, resname, limits, nmax, varargin)`  
  Add ions within a certain region defined by limits.

- **[insert_atom](build_functions/insert_atom.m)**: `insert_atom(atom, new_atom, position)`  
  Insert a new atom at the specified position.

- **[merge_atom](build_functions/merge_atom.m)**: `merge_atom(atom1, Box1, atom2, type, Atom_label, r)`  
  Merge atom structs based on distance criteria.

- **[molid_rotate](build_functions/molid_rotate.m)**: `molid_rotate(atom, Box_dim, MolID, rotate_dim)`  
  Rotate the atom struct based on MolID.

- **[molid_translate](build_functions/molid_translate.m)**: `molid_translate(atom, trans_vec, MolID)`  
  Translate a specific molecule ID.

- **[noupdate_atom](build_functions/noupdate_atom.m)**: `noupdate_atom(atom)`  
  Prevent updating of certain properties in the atom struct.

- **[occupancy_atom](build_functions/occupancy_atom.m)**: `occupancy_atom(atom, Box_dim)`  
  Calculate occupancy of atoms within the box dimensions.

- **[overwrite_atom](build_functions/overwrite_atom.m)**: `overwrite_atom(In_atom, atomtype, resname)`  
  Overwrite atom struct information with new data.

- **[place_atom](build_functions/place_atom.m)**: `place_atom(atom, position)`  
  Place the atom struct at the specified position.

- **[position_molid](build_functions/position_molid.m)**: `position_molid(atom, position_vec, MolID)`  
  Move a molecule ID to a certain position.

- **[protonate_atom](build_functions/protonate_atom.m)**: `protonate_atom(atom, Box_dim, varargin)`  
  Protonate specified sites in the atom struct.

- **[remove_molid](build_functions/remove_molid.m)**: `remove_molid(atom, MolID)`  
  Remove residue with a specific molecule ID.

- **[remove_occypancy_atom](build_functions/remove_occypancy_atom.m)**: `remove_occypancy_atom(atom)`  
  Remove particles with identical coordinates to preceding ones.

- **[remove_residues](build_functions/remove_residues.m)**: `remove_residues(atom, resnames, lo, hi, dim)`  
  Remove residues between specified limits in the simulation box.

- **[remove_resname](build_functions/remove_resname.m)**: `remove_resname(atom, resnames)`  
  Remove residues with specified names.

- **[remove_SOL](build_functions/remove_SOL.m)**: `remove_SOL(atom, atomname, lo, hi, dim)`  
  Remove solvent residues between specified limits.

- **[remove_type](build_functions/remove_type.m)**: `remove_type(atom, typescell)`  
  Remove atom types specified in `typescell`.

- **[rename_atom](build_functions/rename_atom.m)**: `rename_atom(atom, old_name, new_name)`  
  Rename atoms in the atom struct.

- **[rename_type](build_functions/rename_type.m)**: `rename_type(atom, atomtype, new_atomtype, varargin)`  
  Rename atom types in the atom struct.

- **[replicate_atom](build_functions/replicate_atom.m)**: `replicate_atom(atom, Box_dim, replicate)`  
  Replicate the atom struct along orthogonal dimensions.

- **[replace_atom](build_functions/replace_atom.m)**: `replace_atom(new_atom, prev_atom, molid_index)`  
  Replace molecule ID in an atom struct with a new atom struct.

- **[resname_atom](build_functions/resname_atom.m)**: `resname_atom(atom)`  
  Guess residue names for all atom types.

- **[slice_atom](build_functions/slice_atom.m)**: `slice_atom(atom, limits, invert)`  
  Slice the atom struct within specified limits.

- **[slice_box](build_functions/slice_box.m)**: `slice_box(atom, Box_dim, limits)`  
  Slice a simulation box within given limits.

- **[slice_molid](build_functions/slice_molid.m)**: `slice_molid(atom, limits, invert)`  
  Slice molecules within specified limits.

- **[slice_triclinic_atom](build_functions/slice_triclinic_atom.m)**: `slice_triclinic_atom(atom, limits, invert)`  
  Slice a triclinic atom struct within limits.

- **[solvate_atom](build_functions/solvate_atom.m)**: `solvate_atom(limits, density, r, maxsol, solute_atom, varargin)`  
  Generate a solvent structure within specified limits.

- **[spc2tip4p](build_functions/spc2tip4p.m)**: `spc2tip4p(atom)`  
  Convert SPC water molecules to TIP4P model.

- **[spc2tip5p](build_functions/spc2tip5p.m)**: `spc2tip5p(atom)`  
  Convert SPC water molecules to TIP5P model.

- **[spce2tip4p](build_functions/spce2tip4p.m)**: `spce2tip4p(atom)`  
  Convert SPC/E water molecules to TIP4P model.

- **[sphere_atom](build_functions/sphere_atom.m)**: `sphere_atom(atom, Box_dim, center, radius)`  
  Create a spherical region of atoms.

- **[substitute_atom](build_functions/substitute_atom.m)**: `substitute_atom(atom, Box_dim, NumOctSubst, O1, O2, minO2O2_dist, varargin)`  
  Perform isomorphous substitution of atoms.

- **[substitute_NonCentroSymm_atom](build_functions/substitute_NonCentroSymm_atom.m)**: `substitute_NonCentroSymm_atom(atom, Box_dim, replace_type, varargin)`  
  Substitute non-centrosymmetric atoms.

- **[tile_atom](build_functions/tile_atom.m)**: `tile_atom(atom, scale_vec, Box_dim, Resname)`  
  Tile the atom struct in a specific direction.

- **[tip3p2tip4p](build_functions/tip3p2tip4p.m)**: `tip3p2tip4p(atom)`  
  Convert TIP3P water molecules to TIP4P model.

- **[translate_atom](build_functions/translate_atom.m)**: `translate_atom(atom, trans_vec, Resname)`  
  Translate a residue by a specified vector.

- **[translate_molid](build_functions/translate_molid.m)**: `translate_molid(atom, trans_vec, molid)`  
  Translate a molecule ID by a specified vector.

- **[tube_atom](build_functions/tube_atom.m)**: `tube_atom(atom, scale_vec, Box_dim, Resname)`  
  Create a nanotube structure from the atom struct.

- **[update_atom](build_functions/update_atom.m)**: `update_atom(atom)`  
  Update molecule and atom indices in the atom struct.

---

## Translate or Rotate Functions

- **[bend_atom](build_functions/bend_atom.m)**: `bend_atom(atom, Box_dim, Radii)`  
  Bend an atom struct.

- **[center_atom](build_functions/center_atom.m)**: `center_atom(atom, Box_dim, resname, dim)`  
  Center the atom with respect to the `resname` molecule.

- **[condense_atom](build_functions/condense_atom.m)**: `condense_atom(atom, Box_dim, s)`  
  Minimize the box size and remove gaps between `molids`.

- **[molid_rotate](build_functions/molid_rotate.m)**: `molid_rotate(atom, Box_dim, MolID, rotate_dim)`  
  Rotate the atom struct based on MolID.

- **[molid_translate](build_functions/molid_translate.m)**: `molid_translate(atom, trans_vec, MolID)`  
  Translate a specific molecule ID.

- **[place_atom](build_functions/place_atom.m)**: `place_atom(atom, position)`  
  Place the atom struct at the specified position.

- **[position_molid](build_functions/position_molid.m)**: `position_molid(atom, position_vec, MolID)`  
  Move a molecule ID to a certain position.

- **[rotate_atom](general_functions/rotate_atom.m)**: `rotate_atom(atom, rotation_matrix)`  
  Rotate the atom struct by a specified matrix.

- **[translate_atom](build_functions/translate_atom.m)**: `translate_atom(atom, trans_vec)`  
  Translate the atom struct by a specified vector.

- **[translate_molid](build_functions/translate_molid.m)**: `translate_molid(atom, trans_vec, molid)`  
  Translate the molecule ID by a specified vector.

---

## Make Triclinic/Orthogonal Box

- **[frac2atom](general_functions/frac2atom.m)**: `frac2atom(atom, Box_dim, angleparam, angletype)`  
  Transform fractional coordinates to Cartesian coordinates.

- **[orto_atom](general_functions/orto_atom.m)**: `orto_atom(atom, Box_dim)`  
  Transform a triclinic atom struct to an orthogonal one.

- **[triclinic_atom](general_functions/triclinic_atom.m)**: `triclinic_atom(atom, Box_dim, angleparam, angletype)`  
  Transform an orthogonal atom struct to a triclinic one.

---

## Wrap/Unwrap Functions

- **[unwrap_atom](general_functions/unwrap_atom.m)**: `unwrap_atom(atom, Box_dim, dim)`  
  Unwrap the atom struct along the specified dimension.

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

- **[spc2tip4p](build_functions/spc2tip4p.m)**: `spc2tip4p(filename)`  
  Convert a `.gro` or `.pdb` file with SPC water to TIP4P water.

- **[spc2tip5p](build_functions/spc2tip5p.m)**: `spc2tip5p(filename)`  
  Convert a `.gro` or `.pdb` file with SPC water to TIP5P water.

- **[spce2tip4p](build_functions/spce2tip4p.m)**: `spce2tip4p(filename)`  
  Convert a `.gro` or `.pdb` file with SPC/E water to TIP4P water.

- **[tip3p2tip4p](build_functions/tip3p2tip4p.m)**: `tip3p2tip4p(filename)`  
  Convert a `.gro` file with TIP3P water to TIP4P water.

---

## Custom Topology Tools

*(No functions listed under this category. Add relevant functions if available.)*

---

## Miscellaneous FF Functions

- **[print_top](forcefield_functions/print_top.m)**: `print_top(atom, Box_dim, varargin)`  
  Print or generate topology-related data.

- **[import_ff_table](forcefield_functions/import_ff_table.m)**: `import_ff_table(filename, varargin)`  
  Import force field parameter tables.

- **[change_top](forcefield_functions/change_top.m)**: `change_top(atom, Box_dim, varargin)`  
  Modify the topology file or its parameters.

- **[mass_atom_clayff](forcefield_functions/mass_atom_clayff.m)**: `mass_atom_clayff(atom)`  
  Define or calculate masses for Clayff atoms.

- **[smear_charge](forcefield_functions/smear_charge.m)**: `smear_charge(atom, Box_dim, varargin)`  
  Distribute charge across atoms, possibly using charge smearing techniques.

---

## Force Fields

### MINFF, with Atom Types by MHolmboe

- **[minff_atom](forcefield_functions/minff_atom.m)**: `minff_atom(atom, Box_dim, varargin)`  
  Assign MINFF atom types with edge healing.

- **[charge_minff_atom](forcefield_functions/charge_minff_atom.m)**: `charge_minff_atom(atom, Box_dim, varargin)`  
  Set the charge for MINFF atom types.

### CLAYFF, with Atom Types by MHolmboe

- **[charge_atom](forcefield_functions/charge_atom.m)**: `charge_atom(atom, Box_dim, ffname, watermodel, varargin)`  
  Charge the atom according to CLAYFF or INTERFACE force fields.

- **[charge_clayff_2004_atom](forcefield_functions/charge_clayff_2004_atom.m)**: `charge_clayff_2004_atom(atom, Box_dim, varargin)`  
  Set the charge for the original CLAYFF atom types from the Cygan et al., 2004 paper.

- **[charge_clayff_atom](forcefield_functions/charge_clayff_atom.m)**: `charge_clayff_atom(atom, Box_dim, varargin)`  
  Set the charge for CLAYFF atom types.

- **[charge_opls_go_atom](forcefield_functions/charge_opls_go_atom.m)**: `charge_opls_go_atom(atom, Box_dim, varargin)`  
  Set the charge for specific OPLS atom types.

- **[check_clayff_2004_charge](forcefield_functions/check_clayff_2004_charge.m)**: `check_clayff_2004_charge(atom)`  
  Check the charge of the original CLAYFF atom types.

- **[check_clayff_charge](forcefield_functions/check_clayff_charge.m)**: `check_clayff_charge(atom)`  
  Check the charge of the CLAYFF atom types.

- **[check_clayff_H2Odens](forcefield_functions/check_clayff_H2Odens.m)**: `check_clayff_H2Odens(atom, Box_dim)`  
  Check the approximate water density for a CLAYFF system.

- **[check_H2Odens](forcefield_functions/check_H2Odens.m)**: `check_H2Odens(atom, Box_dim)`  
  Compute the water density.

- **[clayff_2004_atom](forcefield_functions/clayff_2004_atom.m)**: `clayff_2004_atom(atom, Box_dim, varargin)`  
  Assign the original CLAYFF atom types by Cygan et al., 2004, with edge healing.

- **[clayff_2004_param](forcefield_functions/clayff_2004_param.m)**: `clayff_2004_param(Atom_label, varargin)`  
  Hold ion and CLAYFF atom type parameters for the original CLAYFF force field.

- **[clayff_atom](forcefield_functions/clayff_atom.m)**: `clayff_atom(atom, Box_dim, varargin)`  
  Assign CLAYFF atom types with edge healing.

- **[clayff_param](forcefield_functions/clayff_param.m)**: `clayff_param(Atom_label, varargin)`  
  Hold ion and CLAYFF atom type parameters.

- **[clayff210_atom](forcefield_functions/clayff210_atom.m)**: `clayff210_atom(atom, Box_dim, varargin)`  
  Assign modified CLAYFF atom types, with edge healing.

- **[clayff211_atom](forcefield_functions/clayff211_atom.m)**: `clayff211_atom(atom, Box_dim, varargin)`  
  Assign modified CLAYFF atom types (faster version), with edge healing.

- **[tweak_charge_atom](forcefield_functions/tweak_charge_atom.m)**: `tweak_charge_atom(atom)`  
  Tweak the charge of the atom struct to correct rounding errors.

### INTERFACE from Heinz 2005, 2013, with Atom Types by MHolmboe

Note that the support for the INTERFACE forcefield by the atom MATLAB functions is limited, and may not work for your systems...

- **[charge_interface_atom](forcefield_functions/charge_interface_atom.m)**: `charge_interface_atom(atom, Box_dim, varargin)`  
  Tries to set the charge for INTERFACE atom types.

- **[charge_interface15_atom](forcefield_functions/charge_interface15_atom.m)**: `charge_interface15_atom(atom, Box_dim, varargin)`  
  Tries to set the charge for INTERFACE 1.5 atom types.

- **[check_interface_charge](forcefield_functions/check_interface_charge.m)**: `check_interface_charge(atom)`  
  Tries to check the charge of the INTERFACE atom types.

- **[check_interface15_charge](forcefield_functions/check_interface15_charge.m)**: `check_interface15_charge(atom)`  
  Check the charge of the INTERFACE 1.5 atom types.

- **[interface_atom](forcefield_functions/interface_atom.m)**: `interface_atom(atom, Box_dim, varargin)`  
  Assign atoms according to the INTERFACE atom types, with modifications for edges.

- **[interface_param](forcefield_functions/interface_param.m)**: `interface_param(Atom_label, water_model)`  
  Hold extended INTERFACE force field parameters.

- **[interface15_atom](forcefield_functions/interface15_atom.m)**: `interface15_atom(atom, Box_dim, varargin)`  
  Assign atoms according to the INTERFACE 1.5 atom types, with modifications for edges.

- **[interface15_param](forcefield_functions/interface15_param.m)**: `interface15_param(Atom_label, water_model)`  
  Hold extended INTERFACE 1.5 force field parameters.

- **[check_interface_charge](forcefield_functions/check_interface_charge.m)**: *(Duplicate)*

- **[check_interface15_charge](forcefield_functions/check_interface15_charge.m)**: *(Duplicate)*

- **[interface15_silica_atom](forcefield_functions/interface15_silica_atom.m)**: `interface15_silica_atom(atom, Box_dim, varargin)`  
  Assign atom types for the INTERFACE 1.5 force field, specific to silica.

### Graphene Oxide Modeled with OPLS/AA

- **[opls_go_atom](forcefield_functions/opls_go_atom.m)**: `opls_go_atom(atom, Box_dim, rmin, rlarge)`  
  Smear out the charge around -OH and epoxide groups in graphene oxide.

- **[oplsaa_go_param](forcefield_functions/oplsaa_go_param.m)**: `oplsaa_go_param(Atom_label, water_model)`  
  Hold the extended OPLS-AA force field parameters for graphite oxide.

- **[charge_opls_go_atom](forcefield_functions/charge_opls_go_atom.m)**: `charge_opls_go_atom(atom, Box_dim, varargin)`  
  Set the charge for specific OPLS atom types.

---

## Writing Topology Files

- **[write_atom_itp](export_functions/write_atom_itp.m)**: `write_atom_itp(atom, Box_dim, filename_out, varargin)`  
  Create and print a GROMACS `.itp` file for MINFF, CLAYFF (or INTERFACE, maybe) force fields.

- **[write_atom_lmp](export_functions/write_atom_lmp.m)**: `write_atom_lmp(atom, Box_dim, filename_out, varargin)`  
  Create and print a LAMMPS data file (`.lj`) for CLAYFF systems.

- **[write_atom_oplsaa_go_itp](export_functions/write_atom_oplsaa_go_itp.m)**: `write_atom_oplsaa_go_itp(atom, Box_dim, filename_out, varargin)`  
  Create and print a GROMACS `.itp` file for OPLS-AA or GO systems.

---

## Bonded and Nonbonded Parameters

- **[bonded_parameters](bonded_parameters.html)**: `bonded_parameters(atom, varargin)`  
  Define bonded parameters for atoms.

- **[nonbonded_parameters](nonbonded_parameters.html)**: `nonbonded_parameters(atom, varargin)`  
  Define nonbonded parameters for atoms.

- **[nonbonded_ff](forcefield_functions/nonbonded_ff.m)**: `nonbonded_ff(atom, varargin)`  
  Define nonbonded force field parameters.

---

## Lennard-Jones and Coulomb Potentials

- **[buckinghamcoul](forcefield_functions/buckinghamcoul.m)**: `buckinghamcoul(atom, Box_dim, varargin)`  
  Calculate interactions using the Buckingham potential and Coulombic forces.

- **[ljcoul_12_6](forcefield_functions/ljcoul_12_6.m)**: `ljcoul_12_6(atom, Box_dim, varargin)`  
  Handle Lennard-Jones 12-6 potential along with Coulombic interactions.

- **[ljcoul_C12C6C4](forcefield_functions/ljcoul_C12C6C4.m)**: `ljcoul_C12C6C4(atom, Box_dim, varargin)`  
  Lennard-Jones and Coulomb potential with C12, C6, C4 terms.

- **[ljcoul_C12C6](forcefield_functions/ljcoul_C12C6.m)**: `ljcoul_C12C6(atom, Box_dim, varargin)`  
  Lennard-Jones and Coulomb potential with C12, C6 terms.

- **[ljcoul](forcefield_functions/ljcoul.m)**: `ljcoul(atom, Box_dim, varargin)`  
  General Lennard-Jones and Coulomb interaction function.

- **[ljcoul_2x2x](forcefield_functions/ljcoul_2x2x.m)**: `ljcoul_2x2x(atom, Box_dim, varargin)`  
  Lennard-Jones and Coulomb potential with 2x factors.

- **[ljcoul_2x](forcefield_functions/ljcoul_2x.m)**: `ljcoul_2x(atom, Box_dim, varargin)`  
  Variation of Lennard-Jones Coulomb with 2x factors.

---

## Additional Resources

- **Documentation**: Each function is documented in its respective HTML file. Refer to the function links above for detailed usage instructions and examples.

- **Examples**: Check the `examples` directory for sample scripts demonstrating how to use the various functions within the library.

- **Contributing**: Contributions are welcome! Please fork the repository and submit a pull request with your enhancements.

- **License**: This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---
