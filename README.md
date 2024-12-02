
# Atomistic Topology Operations in MATLAB

MATLAB scripts for building and manipulating molecular simulation cells to be used with, for instance, molecular dynamics or Monte Carlo simulations.

## Overview
The purpose of the atom scripts and functions is to automate and enable the construction and analysis of complex and anisotropic, multicomponent molecular systems, and generate topological information with bonds, angles, and (optionally) dihedrals.

Download the whole function library from the:
- [Matlab File Exchange](https://se.mathworks.com/matlabcentral/fileexchange/59622-atom)
- [GitHub](https://github.com/mholmboe/atom)

For a comprehensive list of the different functions, look here
- [List all functions by topic](List_all_functions.md)

See also the link to the:
- [HTML documentation](http://moleculargeo.chem.umu.se/wp-content/uploads/file-manager/atom/Documentation/)

---

## What Can These Scripts/Functions Do?

You can easily build new custom-made molecular systems by issuing a few MATLAB commands and calling any of the >100 functions (see the List_all_functions_topics.m|.html)

```matlab

open import_atom; % Tip: see the function file/s for their example input arguments

[oneunitcell,Box_dim] = import_atom('someunitcell.pdb'); % Take some unit cell for instance 

mineral = replicate_atom(oneunitcell, Box_dim, [6 4 1]); % Replicate the unit cell 6x4x1 times 

mineral_ff = minff_atom(mineral, Box_dim); % Assign, for instance, the atom names according to the MINFF or CLAYFF (INTERFACE may work) force fields 

ions = ionize_atom(arguments); % Add ions homogeneously or with preference to a surface 

water = solvate_atom(arguments); % Create solvating water molecules 

FinalSystem = update_atom({mineral_ff, ions, water}); % Add all components together 

% Print the final system to file (pdb|gro|xyz and some other formats supported)
write_atom_pdb(FinalSystem,Box_dim,'filename_out.pdb'); % Write a new output file
write_atom_gro(FinalSystem,Box_dim,'filename_out.gro'); % Write a new output file
write_atom_xyz(FinalSystem,Box_dim,'filename_out.xyz'); % Write a new output file

```

Furthermore, bond/neighbour distances can be compared to the atom-specific distances by comparison to the Revised Shannon radii and the valence of atoms calculated using the Bond Valence Method, using the command:

```matlab

properties = analyze_atom(atom, Box_dim, 2.45); % Where 2.45 [Å] is the maximum neighbor cutoff distance

```

The common theme for the scripts is the use of the MATLAB variable type called struct, which acts as a variable container for the atoms' and molecules' different attributes/properties and uses the so-called '.' (dot) notation. By default, an input structure is put into the variable atom (as in the acronym for the MATLAB library) along with a Box_dim variable, in which atom(1) holds the static attributes/properties of the first atom/particle/site in the molecule, for example, .molid, .type, .resname, .index, .neigh, .bond, .angle, .x, .y, .z, .vx, .vy, .vz (these are the default).

The scripts can be used to manipulate/transform the structures in various ways, making use of the MATLAB struct variable and indexing. Thanks to the atom struct dot notation, one can also filter and manipulate the atom struct with respect to the different attributes using MATLAB's own index-related functions like find, strcmp, ismember, and so on. This even allows us to use dynamic indices of groups of atom.molid, resname, type, index, x, etc., if analyzing a trajectory.

Since the scripts can be used to translate, rotate, and slice molecules or entire mineral slabs, these atom scripts are ideally used to build or manipulate multicomponent systems by adding molecules, ions, and SPC/TIP3P/TIP4P water molecules (i.e., solvating an existing molecule/slab) into a simulation box (and removing any molecular overlaps). For plotting, one can issue the function vmd(atom, Box_dim) if the VMD software is installed and the PATH2VMD() function is properly set. Alternatively, the show_atom(atom, Box_dim) function can plut a full simulation box with thousands of atoms in MATLAB. Most functions, like dist_matrix_atom() and bond_angle_atom(), take PBC into account, which allows for the generation of topologies of molecules with bonds across PBC. There is also some support for triclinic systems using the 'tilt vectors' xy, xz, yz as defined in GROMACS and LAMMPS. Note that several of the functions use two types of neighbor cutoffs to find interacting atoms, one short cutoff for bonded H's (default 1.25Å) and a larger cutoff for all other interactions (default 2.25Å).

Apart from scripts for importing/exporting .xyz and .gro trajectory files, the atom scripts can read and write basic .pdb or .xyz and GROMACS .gro structure files as well as output basic .itp and .psf topology files with bonds and angles across the PBC, for the MINFF forcefield, CLAYFF (Cygan et al., 2004). There is also some support for the  INTERFACE (Heinz, 2005, 2015) force fields.

---

### Working with MATLAB Struct Variables

Import the molecular info into the MATLAB workspace. Here the atom variable becomes an indexed MATLAB struct variable containing all your molecule's attributes like atom names, coordinates, and more. Note that you also get the box dimension variable Box_dim holding the box size if one exists.

```matlab
atom = import_atom_pdb('filename.pdb'); % Example of the pdb import command
```
The scripts also support importing (and some exporting) of certain trajectory formats like .xtc, .trr, .dcd, .pdb, and .gro (the latter two are text formats). See documentation in the import_traj() function, or try, for instance:

```matlab
atom = import_xyz_traj('filename.xyz'); 
atom = import_gro_traj('filename.gro'); 
write_xyz_traj(atom, traj, Box_dim, 'filename_out.xyz'); 
write_gro_traj(atom, traj, Box_dim, 'filename_out.xyz');
```

To write an atom struct to a structure file, run, for instance:

```matlab
write_atom_pdb(atom, Box_dim, 'filename_out.pdb'); 
write_atom_gro(atom, Box_dim, 'filename_out.gro'); 
write_atom_xyz(atom, Box_dim, 'filename_out.xyz');
```
---

### Index Filtering and Dynamic Selection:

```matlab
index = ismember([atom.type], [{'Al', 'Alt', 'Mgo'}]); % Gives a binary (1/0) logical array index = strcmp([atom.type], 'Al'); 

% Try also strncmp or strncmpi 
index = find(strncmpi([atom.type], 'al', 2)); % Will find the indices of 'Al', 'Alt' 

new_atom = atom(index); % Creates a new_atom struct with the filtered/selected atom types 
positive_z_atom = atom([atom.z] > 0); % Finds all atoms with a positive z-coordinate 
first100_atom = atom([atom.index] < 101); % Finds the first 100 atoms in the atom struct 
first100_v2_atom = atom(1:100); % Also finds the first 100 atoms in the atom struct
```
---

### Assign atomtypes according to different forcefields

Label the atomtypes according to different forcefields. Currently supports MINFF and CLAYFF, and the INTERface forcefield to some extent.

```matlab
atom_minff = minff_atom(atom, Box_dim); 
atom_clayff = clayff_atom(atom, Box_dim);
atom_interface = interface_atom(atom, Box_dim); 
```
---

### Write topology files
```matlab
write_atom_itp(atom, Box_dim, filename, 1.2, 1.2, 'minff', 'spce'); % GROMACS topology file, note only bonds and angles
write_atom_psf(atom, Box_dim, filename, 1.2, 1.2, 'clayff'); % Note: only bonds and angles 
```
---

### Adding Water:
Adding water to a box - This function solvates a certain region defined by limits with a water structure with a specified density. r and r-0.5 are the closest distances of Ow and Hw to the solute atoms.

```matlab
SOL_atom = solvate_atom(limits, density, r, maxsol, solute_atom, watermodel); % watermodel (optional) is SPC, TIP3P, TIP4P, TIP5P
```
---
### Merge several Atom Sets or molecular components:

Merging two different atom structs - This function returns the second atom set with non-overlapping atoms that are at least a distance r away from the atoms in the first atom set.

```matlab
new_atom = merge_atom(atom1, Box1, atom2, Box2, type, Atom_label, r);
```
```matlab
% Merge different atom structs, will not check for overlapping atoms
new_atom = update_atom({atom1 atom2 atom3});
```
---

### Manipulation and Transformation:
Supports many different types of operations/manipulatins like translate, rotate, slice etc., using a syntax like:

```matlab
atom = translate_atom(atom, trans_vec, Resname);
atom = wrap_atom(atom, Box_dim);
atom = triclinic_atom(atom, Box_dim, angle_param, angle_type); 
```

---
### Calculating a distance matrix/matrices - This function calculates the distance matrix from the atom struct:

```matlab
dist_matrix = dist_matrix_atom(atom, Box_dim); dist_matrix = dist_matrix_atom(atom1, atom2, Box_dim); % The same but for the distance matrix between particles in atom1 and atom2
```

---


Plus many other functions... see the HTML documentation.

### Original reference to cite:
> Holmboe, M., "Atom: A MATLAB Package for Manipulation of Molecular Systems," Clays and Clay Minerals, Accepted November 2019. DOI:10.1007/s42860-019-00043-y

---

### Contact:
Any issues/concerns/bugs or suggestions for improvements/contributions, etc., email me.
**Michael Holmboe**  
michael.holmboe@umu.se
