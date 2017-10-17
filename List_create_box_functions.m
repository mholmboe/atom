%% List of Box manipulation and build system functions

add2atom(XYZ_labels,XYZ_data,varargin) % This function appends so-called XYZ atomtype labels and XYZ data to an existing atom struct
center_atom(atom,Box_dim,resname,dim) % This function centers the atom with respect to the resname molecule
composition_atom(atom) % This function looks at the composition of the atom struct
concatenate_atom(atom_1,atom_2) % This function concatenats atom sections.
condense_atom(atom,Box_dim,s) % This function tries to minimize the box size and remove gaps between molids?
copy_atom(atom,atomtype,new_atomtype,new_resname,trans_vec,varargin) % This function copies and translates atoms in the atom struct
copy_type(atom,atomtype,new_atomtype,new_resname,trans_vec,varargin) % This function copies and translates types in the atom struct
create_atom(type,resname,limits,nmax,varargin) % Creates new atoms, good for adding ions to a system. Creates atoms within a certain region defined by <limits>
create_grid_atom(atom_label,nM,limits,dim,varargin) % This old function puts ions on a grid plane and adds it to an atom struct
draw_box_atom(Box_dim,LineColor,LineThickness) % Draws a box
duplicate_atom(atom,molID) % This function duplicates residue with molid MolID
element_atom(atom,varargin) % Converts atomtypes to element types. This function replaces the atomtypes names with the element names
frac2atom(atom,Box_dim,angleparam,angletype) % This function transforms fractional coordinates to cartesian
grid2atom(atom_label,nM,limits,dim,varargin) %  grid2atom.m - This puts particles such as ions on a 2D grid (i.e. a plane) % and adds it to an atom struct
insert_atom(atom_in,limits,rotate,r,maxsol,solute_atom,varargin) % - This inserts a molecule from a structure file into a region defined by <limits> with a atom (molecule) % structure
mass_atom(atom) % This function fetches the mass for each atomtype and put it into atom.mass
merge_atom(atom1,Box1,atom2,type,Atom_label,r) % This function returns the atom2 struct with atoms in the atom2 struct with a distance r [1x1 or 1x2] away from the atoms in the atom1 struct. There is also a possibility to use a twin-range cutoff approach (suitable for OH2), by setting r(2) to a smaller value than r(1)
molid_rotate(atom,Box_dim,MolID,rotate_dim) % This function rotate the atom randomly
molid_translate(atom,trans_vec,MolID) % This translates a certain molid
orto_atom(atom,Box_dim) % This transforms a triclinic atom struct to an orthogonal atom struct. Box_dim must look like [lx ly lz 0 0 xy 0 xz yz]
overwrite_atom(In_atom,atomtype,resname) % This function overwrites the atom struct information with new information 
place_atom(atom,position) % This function places the atom struct according to the position vector called position, trying to use the COM of the molecule
position_molid(atom,position_vec,MolID) % This function movies a molid (COM) % to a certain position
reorder_atom_gro(atom,atomlist,Box_dim,filename_out) % This function reorders the atoms in a .gro file
replicate_atom(atom,Box_dim,replicate) % This replicates the atom struct and the orthogonal box dimensions
resname_atom(atom) % This function tries to guess the resname of all atom types
rotate_atom(atom,Box_dim,alfa,beta,gamma) % This function rotate the atom randomly
scale_atom(atom,scale_vec,Box_dim,Resname) % This function scales the coordinates in the atom struct
slice_atom(atom,limits,invert) % This function checks if the coordinates for each time record in XYZ_data is within the specified limits, and if not sets the x,y,z to nan,nan,nan.
slice_molid(atom,limits,invert) % This function checks if the coordinates is within the specified limits, and if not sets the x,y,z to nan,nan,nan.
solvate_atom(limits,density,r,maxsol,solute_atom,varargin) % This function generates a certain region defined by <limits> with a solvent structure of density <density>
substitute_atom(atom,Box_dim,NumOctSubst,O1,O2,minO2O2_dist,varargin) % This scripts performs isomorphous substitution, by replacing some O1->O2 atomtypes and optionally T1->T2 atomtypes
translate_atom(atom,trans_vec,Resname) % This translates the resname by a vector
translate_molid(atom,trans_vec,molid) % This translates the molid by a vector
triclinic_atom(atom,Box_dim,angleparam,angletype) %  triclinic_atom.m - This transforms an orthogonal atom struct to a triclinic with the angles alfa, beta, gamma or tilt factors xy, xz, yz
unwrap_atom(atom,Box_dim,dim) % This function unwraps the atom struct along the dimension dim
update_atom(atom) % This function updates the molid index and the atoms index in the atom struct
wrap_atom(atom,Box_dim) % This wraps the atoms into the orthogonal box
wrap_atom(atom,Box_dim) % This wraps the atoms into the orthogonal box
Wrap_Coord_func(XYZ_data,Box_dim) % This is an old function that wraps atoms 'sticking out' back into the box
wrap_molid(atom,Box_dim) % This function wraps the atom struct into the box
xyz2atom(XYZ_labels,XYZ_data,Box_dim,resname,in_atom) % This function can be used to add XYZ data (like from a .xyz structure file)to the atom struct format

