%%  function descriptions 
% 864_spc.gro % equilibrated SPC water box
% 864_spce.gro % equilibrated SPC/E water box
% 864_tip3p.gro % equilibrated TIP3P water box
% 864_tip4p.gro % equilibrated TIP4P water box
% 864_tip5p.gro % equilibrated TIP5P water box
% 96spc_hex_ice_h.gro % equilibrated SPC hex-ice water box
% 96tip4p_hex_ice_h.gro % equilibrated TIP4P hex-ice water box
add2atom(XYZ_labels,XYZ_data,varargin) % This function appends so-called XYZ atomtype labels and XYZ data to an existing atom struct
analyze_atom(atom,Box_dim,max_H_dist,max_dist) % This function analyzes variofus things of the MMT atom struct
atom2make_ndx(filename,groupname,atomtypes,molid) % This little script can be used to write custom made Gromacs index files
ave_atom(atom) % This function calculates the mean of the atom coordinates
bond_angle_atom(atom,Box_dim,max_short_dist,max_long_dist,varargin) % This function tries to find all bonds and angles of the atom struct 'more' is an optional varargin argument
bond_angle_type(atom1,atom2,Box_dim,rmin,rmax,angle_limit,varargin) % This tries to find all bonds and angles of the atom types
cat_atom.m % This is a special script (and not a function) that imports and appends atom structs into a .gro trajectory file, useful to make a trajectory with varying number of particles
center_atom(atom,Box_dim,resname,dim) % This function centers the atom with respect to the resname molecule
charge_atom(atom,Box_dim,ffname,watermodel,varargin) % This function tries to charge the atom accorind to clayff or interface ff
charge_clayff_atom(atom,Box_dim,varargin) % Sets the charge for Clayff atomtypes 
charge_interface_atom(atom,Box_dim,varargin) % Sets the charge for Interface atomtypes 
charge_interface15_atom(atom,Box_dim,varargin) % Sets the charge for Interface 1.5 atomtypes
charge_opls_go_atom(atom,Box_dim,varargin) % Sets the charge for some specific OPLS atomtypes
check_clayff_charge(atom) % check_clayff_charge.m - This checks the charge of the Clayff atomtypes by Mholmboe
check_clayff_H2Odens(atom,Box_dim) % Specific function that calculates the approx. Wtaer density in heterogenic system
check_H2Odens(atom,Box_dim) % Computes the water density
check_interface_charge(atom) %  This checks the charge of the INTERFACE atomtypes by Mholmboe
check_interface15_charge(atom) %  This checks the charge of the INTERFACE 1.5 atomtypes by Mholmboe
clayff_atom(atom,Box_dim,varargin) % Assigns the Clayff atom types by MHolmboe. Can also 'heal' edges 
clayff_param(Atom_label,varargin) % Holds the ion and Clayff atomtype parameters
clayff_param(Atom_label,water_model) % This function holds the extended Clayff ff parameters
clayffmod_atom(atom,Box_dim,varargin) % Assigns the modififed Clayff atom types. Can also 'heal' edges 
COM_atom(atom,MolID) % This function calculates the COM for certain elements
COM_func(MolID,XYZ_data,Atom_label,XYZ_labels,Box_dim) % This calculates the center of mass for water. Slow due to pbc...
COM_molid(atom,MolID) % This function calculates the COM for certain elements
COM_SOL(MolID,XYZ_data,Atom_label,XYZ_labels,Box_dim) % Computes the COM of SPC water?
composition_atom(atom) % This function looks at the composition of the atom struct
concatenate_atom(atom_1,atom_2) % This function concatenats atom sections.
condense_atom(atom,Box_dim,s) % This function tries to minimize the box size and remove gaps between molids?
CONECT_atom(atom,Box_dim,short_r,long_r) % This prints conect records for pdb files
copy_atom(atom,atomtype,new_atomtype,new_resname,trans_vec,varargin) % This function copies and translates atoms in the atom struct
copy_type(atom,atomtype,new_atomtype,new_resname,trans_vec,varargin) % This function copies and translates types in the atom struct
create_atom(type,resname,limits,nmax,varargin) % Creates new atoms, good for adding ions to a system. Creates atoms within a certain region defined by <limits>
create_grid_atom(atom_label,nM,limits,dim,varargin) % This old function puts ions on a grid plane and adds it to an atom struct
dipoles_atom(Elements,Box_dim) % This function calculates the dipole vector of water. Similar to the COM_func
dist_matrix_atom(atom,Box_dim) % This calculates the distance matrix from the atom struct
draw_box_atom(Box_dim,LineColor,LineThickness) % Draws a box
duplicate_atom(atom,molID) % This function duplicates residue with molid MolID
element_atom(atom,varargin)  % Converts atomtypes to element types. This function replaces the atomtypes names with the element names
element_color(Atom_label) % This function assigns a certain color to each element. Estethic improvements are welcome...
find_bonded_atom(atom,bond_matrix,label1,label2) % This function does a cross check of the bond matrix
frac2atom(atom,Box_dim,angleparam,angletype) % This function transforms fractional coordinates to cartesian
frame2atom(atom,traj,frame,Box_dim,varargin) % This function extracts a frame to the trajectory matrix
gmx_make_ndx(groupname,ind) % This function helps you print custom gromacs .ndx files
grid2atom(atom_label,nM,limits,dim,varargin) %  grid2atom.m - This puts particles such as ions on a 2D grid (i.e. a plane) % and adds it to an atom struct
H2Odens = check_clayff_H2Odens(atom,Box_dim) % Check the approx. water density for a clayff system
hist_atom(atom,s) % This function is used to calculate density profiles in the Z-direction
import_atom_car(filename,varargin) % This function imports .car files from Hendrik Heinz INTERFACE ff distribution, and then tries to write out a Gromacs molecular topology file (.itp) and a new .pdb file
import_atom_gro(filename) % This function imports .gro files into the atom struct
import_atom_multiple_gro(filename,nFrames) % This function import multiple .gro files
import_atom_pdb(filename) % This function imports .pdb files into the atom struct
import_atom_xyz_matlab(filename) % This function imports an .xyz file into the atom struct
import_atom_xyz(filename) % This imports an .xyz file into the atom struct
import_atom(filename) % import_atom.m - This imports a .xyz/.gro/.pdb file and puts the data in a structure variable called atom
import_gro_traj(filename,varargin) % This function imports an strcture and an .gro trajectory file
import_pdb_traj(filename,varargin) % This function imports an strcture and an .pdb trajectory file.
import_traj(filenameconf,filenametraj) % This function imports an strcture and an dcd, trr, xtc, xyz or gro  trajectory file.
import_trr(filenameconf,filenametraj) % This function imports an strcture and an trr  trajectory file
import_trrv2(filenameconf,filenametraj) % This function imports an strcture and an trr  trajectory file
import_xtc(filenameconf,filenamextc) % import_atom_xtc.m - This imports a structure file and a xtc file
import_xtcv2(filenameconf,filenamextc) % import_atom_xtc.m - This imports a structure file and a xtc file
import_xvg(filename) % This function imports a Gromacs .xvg file
import_xyz_traj(filenametraj) % import_xyz_traj.m - This imports an strcture and an .xyz trajectory file.
import_xyz(filename) % This function imports an .xyz file. Atom types should be made of letters, not numbers... Try the import_atom_xyz function instead...
insert_atom(atom_in,limits,rotate,r,maxsol,solute_atom,varargin) % - This inserts a molecule from a structure file into a region defined by <limits> with a atom (molecule) % structure
interface_atom(atom,Box_dim,varargin) % This function tries to assign all atoms according to the interface atom types (with modified atom names by MHolmboe), with some modifications for edges...
interface_param(Atom_label,water_model) %  This holds the extended INTERFACE ff parameters
interface15_atom(atom,Box_dim,varargin) % This function tries to assign all atoms according to the interface1.5 atom types (with modified atom names by MHolmboe), with some modifications for edges...
interface15_param(Atom_label,water_model) %  This holds the extended INTERFACE1.5 ff parameters
keep_atom(atom,resname) % keep_atom.m - This removes all but resname
keep_resname(atom,resnames) % keep_resname.m - This removes all but the resnames
lmp_atom_style_full_func(fid,Atom_label,Charge,XYZ_labels,XYZ_data) % This creates and prints the Atoms section properties in the LAMMPS data file .lj file according to atom style full, without image flags
mass_atom_clayff(atom,varargin) % This function fetches the atom weight from the clayff and interface ff's
mass_atom(atom) % This function fetches the mass for each atomtype and put it into atom.mass
median_atom(atom) % This function calculates the median position of the atom struct
merge_atom(atom1,Box1,atom2,type,Atom_label,r) % This function returns the atom2 struct with atoms in the atom2 struct with a distance r [1x1 or 1x2] away from the atoms in the atom1 struct. There is also a possibility to use a twin-range cutoff approach (suitable for OH2), by setting r(2) to a smaller value than r(1)
molid_rotate(atom,Box_dim,MolID,rotate_dim) % This function rotate the atom randomly
molid_translate(atom,trans_vec,MolID) % This translates a certain molid
neigh_atom(atom,Box_dim,rmax,varargin) % This function checks which neighbors each atom has and ouputs their info
new_neigh_atom(atom,Box_dim,rmax,varargin) % Not finished yet...
opls_go_atom(atom,Box_dim,rmin,rlarge) % This function tries to smear out the charge at around -OH and epoxides in GO
oplsaa_go_param(Atom_label,water_model) % This custom function holds the extended oplsaa_aa ff for graphite oxide
orto_atom(atom,Box_dim) % This transforms a triclinic atom struct to an orthogonal atom struct. Box_dim must look like [lx ly lz 0 0 xy 0 xz yz]
overwrite_atom(In_atom,atomtype,resname) % This function overwrites the atom struct information with new information 
PATH2VMD() % The VMD path on your computer
place_atom(atom,position) % This function places the atom struct according to the position vector called position, trying to use the COM of the molecule
plot_atom(atom,Box_dim,varargin) % This function draws the atom struct in 3D. Its very simplistic with no cool features
position_molid(atom,position_vec,MolID) % This function movies a molid (COM) % to a certain position
radius_atom(atom,ffname,watermodel) % This function fetches the ion radius from clayff or interface or interface2015 ff's and
radius_ion(Atom_label) % This function fetches the ionic radius, originally taken from the link below
remove_molid(atom,MolID) %  remove_molid.m - This removes residue with molid MolID = [1 2 3 .....]
remove_residues(atom,resnames,lo,hi,dim) % This function section is used to remove residues in the simulation box between limits lo and hi
remove_resname(atom,resnames) % This function removes residue with molid MolID, resnames = {'SOL' 'Protein'}
remove_SOL(atom,atomname,lo,hi,dim) %  This section is used to remove residues in the simulation box between limits lo and hi
remove_type(atom,typescell) % This function removes atomtypes with types as in typescell = {'OW' 'HW1' 'HW2'}
rename_type(atom,atomtype,new_atomtype,varargin) % This function renames atoms in the atom
reorder_atom_gro(atom,atomlist,Box_dim,filename_out) % This function reorders the atoms in a .gro file
replicate_atom(atom,Box_dim,replicate) %  replicate_atom.m This replicates the atom struct and the orthogonal box dimensions
resname_atom(atom) % This function tries to guess the resname of all atom types
rotate_atom(atom,Box_dim,alfa,beta,gamma) % This function rotate the atom randomly
scale_atom(atom,scale_vec,Box_dim,Resname) % This function scales the coordinates in the atom struct
slice_atom(atom,limits,invert) % This function checks if the coordinates for each time record in XYZ_data is within the specified limits, and if not sets the x,y,z to nan,nan,nan.
slice_molid(atom,limits,invert) % This function checks if the coordinates is within the specified limits, and if not sets the x,y,z to nan,nan,nan.
solvate_atom(limits,density,r,maxsol,solute_atom,varargin) % This function generates a certain region defined by <limits> with a solvent structure of density <density>
sort_atom(atom) % sort_atom.m - This section orders to atoms with respect to z
sort_molid(Molid) % This function sorts the molecular indexes in an ascending order
spc2tip4p(filename) % This function converts a .gro or .pdb file with spc water to some tip4p water
spc2tip5p(filename) % This function converts a .gro or .pdb file with spc water to some tip5p water
spce2tip4p(filename) % This function converts a .gro or .pdb file with spc water to some tip4p water
substitute_atom(atom,Box_dim,NumOctSubst,O1,O2,minO2O2_dist,varargin) % This scripts performs isomorphous substitution, by replacing some O1->O2 atomtypes and optionally T1->T2 atomtypes
tip3p2tip4p(filename) % This function converts a .gro file with spc water to some tip4p water
translate_atom(atom,trans_vec,Resname) % translate_atom.m - This translates the resname by a vector
translate_molid(atom,trans_vec,molid) %  translate_molid.m - This translates the molid by a vector
triclinic_atom(atom,Box_dim,angleparam,angletype) %  triclinic_atom.m - This transforms an orthogonal atom struct to a triclinic with the angles alfa, beta, gamma or tilt factors xy, xz, yz
tweak_charge_atom(atom) % This function tries to tweak the charge of the atom struct in case of rounding errors
unwrap_atom(atom,Box_dim,dim) % This function unwraps the atom struct along the dimension dim
update_atom(atom) % This function updates the molid index and the atoms index in the atom struct
vmd(atom,Box_dim) % This function plots the atom struct
wrap_atom(atom,Box_dim) % This wraps the atoms into the orthogonal box
Wrap_Coord_func(XYZ_data,Box_dim) % This is an old function that wraps atoms 'sticking out' back into the box
wrap_molid(atom,Box_dim) % This function wraps the atom struct into the box
write_atom_gro(atom,Box_dim,filename_out) %  write_atom_gro.m - This writes a gro file. Does it also write velocities?
write_atom_itp(atom,Box_dim,filename,varargin) % This script creates and prints a gromacs .itp file. Works best for clayff or interface ff with spc, spce or tip3p
write_atom_lmp(atom,Box_dim,filename,varargin) % This script creates and prints a lammps data file (.lj). Works best for Clayff systems
write_atom_mol2(atom,Bond_index,Box_dim,filename_out) % This function writes a .mol2 file from the atom struct
write_atom_multiple_gro(atom,traj,filename_out) % This function writes a .gro trajectory
write_atom_oplsaa_go_itp(atom,Box_dim,filename,varargin) % This custom made script creates and prints a gromacs .itp file for 
write_atom_pdb(atom,Box_dim,filename_out) % This function writes an .pdb file from the atom struct using Gromacs
write_atom_pqr(atom,Box_dim,filename_out,varargin) % This function writes an .pqr file from the atom struct
write_atom_psf(atom,Box_dim,filename_out,varargin) % This function writes an .psf file from the atom struct
write_atom_xyz(atom,Box_dim,filename_out) % This function writes an XYZ file from the atom struct
write_atom(atom,Box_dim,filename_out,varargin) % This function tries to write various files for you. Works best for systems designed for Clayff...
write_gro_traj(atom,traj,Box_dim,filename_out) % This function writes a .gro trajectory 
write_pdb_traj(atom,traj,Box_dim,filename_out) % This function writes a .pdb trajectory 
write_xyz_traj(atom,traj,Box_dim,filename_out) % This function writes a .xyz trajectory 
xyz2atom(XYZ_labels,XYZ_data,Box_dim,resname,in_atom) % This function can be used to add XYZ data (like from a .xyz structure file)to the atom struct format
