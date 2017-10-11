%% Atomistic Trajectory Operations in Matlab, i.e. matlab scripts for manipulation of molecular dynamics or monte carlo simulation systems and trajectories. 

%% The purpose of these 50+ 'atom' scripts and functions is to automate and enable construction and analysis of complex and anisotropic multicomponent molecular systems, and generate topological information with bonds and angles. Note that 'atom' in this context represents a Matlab struct variable, holding the information of a molecule/structure, where each particle can be accessed through the '.' syntax, i.e. [atom(1).x] gives the x-coordinate of the first particle/atom. 

%% The atom scripts can import/exporting basic .xyz, .pdb .xyz and gromacs .gro structure files, as well as output simple .itp and .psf topology files with bonds and angles across the PBC. They can also be used to manipulate/transform the molecules/structures in various ways making use of the Matlab struct variable type with logical indexing. The atom scripts can be used to build and plot multicomponent systems, by adding molecules, ions and SPC/TIP3P/TIP4P/TIP5P water molecules (ie solvating an existing molecule/slab) into a simulation box, and remove molecular overlap. For plotting one can call vmd(atom,Box_dim) if the VMD software is installed and the PATH2VMD() function is properly set. Else the quick-and-dirty plot_atom(atom,Box_dim) function rapidly plots a full simulation box with thousands of atoms. Most functions, like dist_matrix_atom() and bond_angle_atom() takes PBC into account. There is also some support for triclinic support using the ’tilt vectors’ xy, xz, yz as defined in Gromacs and Lammps. 

%For a more complete collection of atom scripts, or to report any issues/concerns/bugs/ suggestions for improvements/contributions etc, email me... 

% Michael Holmboe 
% michael.holmboe@umu.se 
% Chemistry department 
% Umeå University 
% Sweden 

%% How to use? 
% To read a structure file into matlab (check the variable explorer) 
atom=import_atom(filename) % where filename is a .pdb, .xyz or .gro file 

% or... 
atom=import_atom_pdb(filenamepdb) 
atom=import_atom_xyz(filenamexyz) 
atom=import_atom_gro(filenamegro) 
% Note that you also get the box dimension variable Box_dim 
% or... 
import_xyz_traj(filenamexyz) | write_xyz_traj(atom,traj,Box_dim,filenameout) 
import_gro_traj(filenamegro) | write_gro_traj(atom,traj,Box_dim,filenameout) 

%% To write a atom struct to a new topology or structure file 
write_atom_psf(atom,Box_dim,filename,1.2,1.2,'clayff') % note only bonds and angles 
write_atom_itp(atom,Box_dim,filename,1.2,1.2,'clayff','spce') % Gromacs topology file, note only bonds and angles 
write_atom_pdb(atom,Box_dim,filename) 
write_atom_gro(atom,Box_dim,filename) 
write_atom_xyz(atom,Box_dim,filename) 

% One can filter the atom struct with respect to molid, resname, atomtype, 
% index, coordinates and so on. This allows us to manipulate an atom struct on % the atomic, molecule and molecular type level. This also allows us to use 
% 'dynamic indexs' of groups of atom.{molid/resname/type/index/} when 
% analyzing a trajectory for instance. Some basic examples: 
index=ismember([atom.type],[{'Al' 'Alt' 'Mgo'}]) % gives a binary (1/0) logical array 
index=strcmp([atom.type],'Al') % try also strncmp or strncmpi? 
index=find(strncmpi([atom.type],'al',2) % Will find the indexes of 'Al' 'Alt? 
new_atom=atom(index) % This creates a new_atom struct with the filtered/selected atomtypes 
positive_z_atom=atom([atom.z]>0) % finds all atoms with a positve z-coordinate 
first100_atom=atom([atom.index]<101) % finds the first 100 atoms in the atom struct 
first100_v2_atom=atom(1:100) % also finds the first 100 atoms in the atom struct 

%% Adding water to a box 
% - This function SOLvates a certain region defined by limits with a water 
% structure with density. r and r-0.5 is the closest distance of Ow and Hw 
% to the solute atoms 
SOL_atom = solvate_atom(limits,density,r,maxsol,solute_atom,watermodel) % watermodel (optional) is spc, tip3p, tip4p, tip5p 

%% Merging two different atom structs 
% - This function returns the second atom set with non-overlapping atoms 
% with a distance r away from the atoms in the first atom set 
new_atom = merge_atom(atom1,Box1,atom2,Box2,type,Atom_label,r) 

%% Calculating a distance matrix/es 
% - This function calculates the distance matrix from the atom struct 
dist_matrix = dist_matrix_atom(atom,Box_dim) 
% - The same but for the distance matrix between atoms in atom1 and atom2 
dist_matrix = dist_matrix_atom(atom1,atom2,Box_dim) 

%% Other example functions (see complete list in function_descriptions.m) 
atom = place_atom(atom,position) 
atom = translate_atom(atom,trans_vec,Resname) 
atom = wrap_atom(atom,Box_dim) 
atom = triclinic_atom(atom,Box_dim,angleparam,angletype) 
% plus 50+ other functions...
