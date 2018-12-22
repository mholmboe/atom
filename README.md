# atom
Matlab scripts for reading/writing/creating/manipulating molecular systems

% Atomistic Topology Operations in Matlab, scripts for manipulation of molecular dynamics or monte carlo simulation systems.
 
% The purpose of these atom scripts and functions is to automate and enable efficient construction/manipulation and analysis of complex and multicomponent molecular systems, and generate topological information with bonds and angles etc. In the near future, I will also add scripts for trajectory analysis.

% For lists of all available functions by category, see inside these files:
List_all_functions.m
List_build_functions.m
List_export_functions.m
List_general_functions.m
List_import_functions.m
List_forcefield_functions.m
 
% The atom scripts can read and write basic .pdb .xyz and gromacs .gro structure files as well as write basic .itp and .psf topology files with bonds and angles across the PBC. The can also manipulate/transform the structures in various ways making use of the Matlab struct variable and indexing. The atom scripts can be used to build and plot multicomponent systems, by adding molecules, ions and SPC/TIP3P/TIP4P water molecules or other solvents (ie solvating an existing molecule/slab) into a simulation box, and remove molecular overlap. For plotting one can call vmd(atom,Box_dim) if the VMD software is also installed and the PATH2VMD() function is properly set. Else the very quick-and-dirty plot_atom(atom,Box_dim) or plot_density_atom(atom,Box_dim) functions which rapidly plots a full simulation box with thousands of atoms. Most functions, like dist_matrix_atom() and bond_angle_atom() takes PBC into account, which allows for generation of topologies of molecules with bonds across PBC. There is also some support for triclinic support using the tilt vectors xy, xz, yz.

% Any issues/concerns/bugs, email me. 
 
% Michael Holmboe 
% michael.holmboe@umu.se 
% Chemistry department
% Umeå University 
% Sweden
 
% How to use? 
% 
% To read a structure file into matlab (check the variable explorer) 
atom=import_atom(filename) % where filename is a .pdb | .xyz  | .pdb file 
 
% or... 
atom=import_atom_pdb(filenamepdb)
atom=import_atom_gro(filenamegro) 
atom=import_atom_xyz(filenamexyz)
% Note that you get a lot more info then just the atom struct variable, like the box dimension variable Box_dim 
 
% To write a atom struct to a new topology or structure file
write_atom_psf(atom,Box_dim,filename,1.2,1.2,'clayff') % note only bonds and angles
write_atom_itp(atom,Box_dim,filename,1.2,1.2,'clayff','spce') % Gromacs topology file, note only bonds and angles
write_atom_pdb(atom,Box_dim,filename)
write_atom_cif(atom,Box_dim,filename) % Prints only fractional coordinates
write_atom_gro(atom,Box_dim,filename) 
write_atom_xyz(atom,Box_dim,filename)
% Adding water to a box 
% - This function SOLvates a certain region defined by limits with a water 
% structure with density. r (and r-0.5 for H) is the closest distance of solvent atoms
% to the (optional) solute atoms
SOL_atom = solvate_atom(limits,density,r,maxsol) % limits can be [10] | [10 20 30] | [10 20 30 40 50 60]
SOL_atom = solvate_atom(limits,density,r,maxsol,solute_atom) % solute_atom struct can be empty, like []
SOL_atom = solvate_atom(limits,density,r,maxsol,solute_atom,'tip4p') % spc | tip3p | tip4p | tip5p
 
% One can filter the atom struct with respect to molid, resname, atomtype, index, coordinates and so on. This allows us to manipulate an atom struct on % the atomic, molecule and molecular type level. This also allows us to use  'dynamic indexs' of groups of atom.{molid/resname/type/index/} when analyzing a trajectory for instance. Some basic examples: 
index=ismember([atom.type],[{'Al' 'Alt' 'Mgo'}]) % gives a binary (1/0) logical array 
index=strcmp([atom.type],'Al') % try also strncmp or strncmpi? 
index=find(strncmpi([atom.type],'al',2) % Will find the indexes of 'Al' 'Alt? 
new_atom=atom(index) % This creates a new_atom struct with the filtered/selected atomtypes 
 
positive_z_atom=atom([atom.z]>0) % finds all atoms with a positve z-coordinate
first100_atom=atom([atom.index]<101) % finds the first 100 atoms in the atom struct 
first100_v2_atom=atom(1:100) % also finds the first 100 atoms in the atom struct
 
% Merging two different atom structs 
% - This function returns the second atom set with non-overlapping atoms
% with a distance r away from the atoms in the first atom set
new_atom = merge_atom(atom1,Box1,atom2,Box2,type,Atom_label,r)
 
% Calculating a distance matrix/es 
% - This function calculates the distance matrix from the atom struct . Other versions using cell lists also exist.
dist_matrix = dist_matrix_atom(atom,Box_dim)
dist_matrix = dist_matrix_atom(atom1,atom2,Box_dim)
 
% Other example functions (see complete list in function_descriptions.m)
[twotheta,intensity] = xrd_atom( filename |  atom,Box_dim )
atom = translate_atom(atom,trans_vec,Resname)
atom = wrap_atom(atom,Box_dim)
atom = orto_atom(atom,Box_dim)
atom = triclinic_atom(atom,Box_dim,angleparam,angletype)

