Atomistic Topology Operations in Matlab, i.e. matlab scripts for building and manipulating molecular simulation cells to be used with for instance dynamics or monte carlo simulations.

The purpose of the atom scripts and functions is to automate and enable construction and analysis of complex and anisotropic multicomponent molecular systems, and generate topological information with bonds and angles.

See also the html documentation...

What can these scripts/functions do for me?

You can easily build new custom made molecular systems by issuing a few Matlab commands. Furthermore, bond/neighbour distances can be compared to the atom-specific distances by comparison to the Revised Shannon radii and the valence of atoms calculated using the Bond valence method.

The common theme for the scripts is the use of the Matlab variable type called struct, which acts as a variable container for the atoms and molecules different attributes/properties, which uses the so-called '.' (dot) notation. By default an input structure is put into the variable atom (as in the acronym for the matlab library) along a Box_dim variable, inwhere atom(1) holds the static attributes/properties of the first atom/particle/site in the molecule, for example  .molid|.type|.resname|.index|.neigh|.bond|.angle|.x|.y|.z|.vx|.vy|.vz (these are  default).

The scripts can be used to manipulate/transform the structures in various ways making use of the Matlab struct variable and indexing. Thanks to the atom struct dot notation, one can also filter and manipulate the atom struct with respect to the different attributes using matlab own index-related functions like find, strcmp, ismember and so on. This even allows us to use dynamic indexes of groups of  atom.molid|resname|type|index|x etc. if analyzing a trajectory.

Since the scripts can be used to translate, rotate and slice molecules or entire mineral slabs, these atom scripts are ideally used to build or manipulate multicomponent systems by adding molecules, ions and SPC/TIP3P/TIP4P water molecules (ie solvating an existing molecule/slab) into a simulation box (and remove any molecular overlaps). For plotting, one can issue the function vmd(atom,Box_dim) if the VMD software is installed and the PATH2VMD() function is properly set. Else the quick-and-dirty plot_atom(atom,Box_dim) function rapidly plots a full simulation box with thousands of atoms. Most functions, like dist_matrix_atom() and bond_angle_atom() takes PBC into account, which allows for generation of topologies of molecules with bonds across PBC. There is also some support for triclinic support using the ?tilt vectors? xy, xz, yz as defined in Gromacs and Lammps. Note that several of the functions uses two types of neighbor cutoff's to find interacting atoms, one short cutoff for bonded H's (default 1.25Å) and a larger cutoff for all other interactions (default 2.25Å).

Apart from scripts for import/exporting .xyz and .gro trajectory files, the atom scripts can read and write basic .pdb or .xyz and Gromacs .gro structure files as well as output basic .itp and .psf topology files with bonds and angles across the PBC, for the Clayff (Cygan et al 2004) and the INTERFACE (Heinz, 2005) forcefields.

Any issues/concerns/bugs or suggestions for improvements/contributions etc, email me.

Michael Holmboe
michael.holmboe@umu.se
http://moleculargeo.chem.umu.se/codes/atom-scripts/


