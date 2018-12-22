Atomistic Topology Operations in Matlab, i.e. matlab scripts for building and manipulating molecular simulation cells to be used with for instance dynamics or monte carlo simulations.

The purpose of the atom scripts and functions is to automate and enable construction and analysis of complex and anisotropic multicomponent molecular systems, and generate topological information with bonds and angles.

See also the html documentation...

What can these scripts/functions do for me?

You can easily build new custom made molecular systems by issuing a few Matlab commands, typically by issuing something like:

oneunitcell=import_atom('someunitcell.pdb') % Take some unit cell for instance
mineral=replicate_atom(oneunitcell,Box_dim,[6 4 1]) % replicate the unit cell 6x4x1 times
mineral_ff=clayff_atom(mineral,...) % Assign for instance the atom names according to the CLAYFF or INTERFACE forcefields
ions=ionize_atom(...) % add ions homogeneously or with preference to a surface
water=solvate_atom(...) % Create solvating water molecules Final
FinalSystem=update_atom({mineral_ff ions water}) % Add all components together
write_atom_pdb|gro(FinalSystem,...) % Write a new output file

Furthermore, bond/neighbour distances can be compared to the atom-specific distances by comparison to the Revised Shannon radii and the valence of atoms calculated using the Bond valence method, using the command:

properties = properties_atom(atom,Box_dim,2.25) % Where 2.25 [Å] is the maximum neighbor cutoff distance

The common theme for the scripts is the use of the Matlab variable type called struct, which acts as a variable container for the atoms and molecules different attributes/properties, which uses the so-called '.' (dot) notation. By default an input structure is put into the variable atom (as in the acronym for the matlab library) along a Box_dim variable, inwhere atom(1) holds the static attributes/properties of the first atom/particle/site in the molecule, for example  .molid|.type|.resname|.index|.neigh|.bond|.angle|.x|.y|.z|.vx|.vy|.vz (these are  default).

The scripts can be used to manipulate/transform the structures in various ways making use of the Matlab struct variable and indexing. Thanks to the atom struct dot notation, one can also filter and manipulate the atom struct with respect to the different attributes using matlab own index-related functions like find, strcmp, ismember and so on. This even allows us to use dynamic indexes of groups of  atom.molid|resname|type|index|x etc. if analyzing a trajectory.

Since the scripts can be used to translate, rotate and slice molecules or entire mineral slabs, these atom scripts are ideally used to build or manipulate multicomponent systems by adding molecules, ions and SPC/TIP3P/TIP4P water molecules (ie solvating an existing molecule/slab) into a simulation box (and remove any molecular overlaps). For plotting, one can issue the function vmd(atom,Box_dim) if the VMD software is installed and the PATH2VMD() function is properly set. Else the quick-and-dirty plot_atom(atom,Box_dim) function rapidly plots a full simulation box with thousands of atoms. Most functions, like dist_matrix_atom() and bond_angle_atom() takes PBC into account, which allows for generation of topologies of molecules with bonds across PBC. There is also some support for triclinic support using the ?tilt vectors? xy, xz, yz as defined in Gromacs and Lammps. Note that several of the functions uses two types of neighbor cutoff's to find interacting atoms, one short cutoff for bonded H's (default 1.25Å) and a larger cutoff for all other interactions (default 2.25Å).

Apart from scripts for import/exporting .xyz and .gro trajectory files, the atom scripts can read and write basic .pdb or .xyz and Gromacs .gro structure files as well as output basic .itp and .psf topology files with bonds and angles across the PBC, for the Clayff (Cygan et al 2004) and the INTERFACE (Heinz, 2005) forcefields.

Any issues/concerns/bugs or suggestions for improvements/contributions etc, email me.

Michael Holmboe
michael.holmboe@umu.se

So... how to use?

To read a structure file into Matlab (check Matlab's variable explorer)

atom=import_atom(filename) % where filename is a .pdb, .xyz or .gro file

or...

atom=import_atom_pdb(filenamepdb)
atom=import_atom_xyz(filenamexyz)
atom=import_atom_gro(filenamegro)

Where the atom variable becomes an indexed Matlab struct variable - containing all your molecules attributes like atom names, coordinates and more. Note that you also get the box dimension variable Box_dim holding the box size if such exists.

The  scripts also support importing (and some exporting) of certain trajectory formats like .xtc, .trr, .dcd , .pdb and .gro (the latter two are text formats). See documentation in the import_traj() function, or try for instance:

atom = import_xyz_traj(filenamexyz)
atom = import_gro_traj(filenamegro)
write_xyz_traj(atom,traj,Box_dim,filenameout)
write_gro_traj(atom,traj,Box_dim,filenameout)

To write a atom struct to a structure file, run for instance:

write_atom_pdb(atom,Box_dim,filename)
write_atom_gro(atom,Box_dim,filename)
write_atom_xyz(atom,Box_dim,filename)

Some basic examples on using 'index' filtering...:

index=ismember([atom.type],[{'Al' 'Alt' 'Mgo'}]) % gives a binary (1/0) logical array
index=strcmp([atom.type],'Al') % try also strncmp or strncmpi?
index=find(strncmpi([atom.type],'al',2) % Will find the indexes of 'Al' 'Alt'
new_atom=atom(index) % This creates a new_atom struct with the filtered/selected atomtypes
positive_z_atom=atom([atom.z]>0) % finds all atoms with a positve z-coordinate
first100_atom=atom([atom.index]<101) % finds the first 100 atoms in the atom struct
first100_v2_atom=atom(1:100) % also finds the first 100 atoms in the atom struct

Adding water to a box - This function SOLvates a certain region defined by limits with a water structure with density. r and r-0.5 is the closest distance of Ow and Hw to the solute atoms
SOL_atom = solvate_atom(limits,density,r,maxsol,solute_atom,watermodel) % watermodel (optional) is spc, tip3p, tip4p, tip5p

Merging two different atom structs  - This function returns the second atom set with non-overlapping atoms with a distance r away from the atoms in the first atom set
new_atom = merge_atom(atom1,Box1,atom2,Box2,type,Atom_label,r)

Calculating a distance matrix/es  - This function calculates the distance matrix from the atom struct

dist_matrix = dist_matrix_atom(atom,Box_dim) 
dist_matrix = dist_matrix_atom(atom1,atom2,Box_dim) % The same but for the distance matrix between particles in atom1 and atom2

In case you want to use the Clayff forcefield for instance, you could assign the clayff atom types with clayff_atom and write topology files with:
write_atom_psf(atom,Box_dim,filename,1.2,1.2,'clayff') % note only bonds and angles write_atom_itp(atom,Box_dim,filename,1.2,1.2,'clayff','spce') % Gromacs topology file, note only bonds and angles

Other example functions (see complete list in function_descriptions.m)

atom = translate_atom(atom,trans_vec,Resname)
atom = wrap_atom(atom,Box_dim)
atom = triclinic_atom(atom,Box_dim,angleparam,angletype)
.
.
.

plus many other functions... see the html documentation.
