%% List of building and simulation cell manipulation functions
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%

%% Specific atom struct functions
% # <add2atom.html add2atom(XYZ_labels,XYZ_data,varargin)> % This function appends so-called XYZ atomtype labels and XYZ data to an existing atom struct
% # <xyz2atom.html xyz2atom(XYZ_labels,XYZ_data,Box_dim,resname,in_atom)> % This function can be used to add XYZ data (like from a .xyz structure file)to the atom struct format
% # <analyze_atom.html analyze_atom(atom,Box_dim,max_H_dist,max_dist)> % This function analyzes variofus things of the MMT atom struct
% # <concatenate_atom.html concatenate_atom(atom_1,atom_2)> % This function concatenats atom sections.
% # <composition_atom.html composition_atom(atom)> % This function looks at the composition of the atom struct
% # <number_type.html numer_type(atom,varargin)> % This function numbers the atom types, like H1, H2, H3...
% # <rename_type.html rename_type(atom,atomtype,new_atomtype,varargin)> % This function renames atoms in the atom
% # <reorder_atom_gro.html reorder_atom_gro(atom,atomlist,Box_dim,filename_out)> % This function reorders the atoms in a .gro file
% # <reorder_atom.html reorder_atom(atom,atomlist)> % This function reorders the atoms in the atom struct
% # <resname_atom.html resname_atom(atom)> % This function tries to guess the resname of all atom types
% # <round_atom.html round_atom(atom,Box_dim,varargin)> % This function rounds the coordinates in the atom struct
% # <sort_atom.html sort_atom(atom)> % sort_atom.m - This section orders to atoms with respect to z
% # <sort_molid.html sort_molid(Molid)> % This function sorts the molecular indexes in an ascending order
% # <scale_atom.html scale_atom(atom,scale_vec,Box_dim,Resname)> % This function scales the coordinates in the atom struct
% # <update_atom.html update_atom(atom)> % This function updates the molid index and the atoms index in the atom struct


%% Add/create/replicate/overwrite atoms
% # <copy_atom.html copy_atom(atom,atomtype,new_atomtype,new_resname,trans_vec,varargin)> % This function copies and translates atoms in the atom struct
% # <copy_type.html copy_type(atom,atomtype,new_atomtype,new_resname,trans_vec,varargin)> % This function copies and translates types in the atom struct
% # <create_atom.html create_atom(type,resname,limits,nmax,varargin)> % Creates new atoms, good for adding ions to a system. Creates atoms within a certain region defined by <limits>
% # <create_grid_atom.html create_grid_atom(atom_label,nM,limits,dim,varargin)> % This old function puts ions on a grid plane and adds it to an atom struct
% # <duplicate_atom.html duplicate_atom(atom,molID)> % This function duplicates residue with molid MolID
% # <fuse_atom.html fuse_atom(atom,Box_dim,varargin)> % This function tries to fuse all sites within a certain cutoff rmax, typically 0.85 Å
% # <grid2atom.html grid2atom(atom_label,nM,limits,dim,varargin)> %  grid2atom.m - This puts particles such as ions on a 2D grid (i.e. a plane)> % and adds it to an atom struct
% # <heal_atom.html heal_atom(atom,Box_dim,ind,varargin)> % This function heals sites in the atom struct given by the index vector ind, by adding a certain atomtype to a new atom struct called healed_atom. It does so by placing the new atom type opposite to the mean position of all neighbours within rcut [Å] of the healed site.
% # <ionize_atom.html ionize_atom(type,resname,limits,nmax,varargin)> % This function adds ions within a certain region defined by <limits> or close to a surface in an exponential fashion
% # <merge_atom.html merge_atom(atom1,Box1,atom2,type,Atom_label,r)> % This function returns the atom2 struct with atoms in the atom2 struct with a distance r [1x1 or 1x2] away from the atoms in the atom1 struct. There is also a possibility to use a twin-range cutoff approach (suitable for OH2), by setting r(2) to a smaller value than r(1)
% # <overwrite_atom.html overwrite_atom(In_atom,atomtype,resname)> % This function overwrites the atom struct information with new information 
% # <protonate_atom.html protonate_atom(atom,Box_dim,varargin)> % This function protonates the sites in the atom struct given by the index vector ind by adding a H's to a new H atom struct.
% # <replicate_atom.html replicate_atom(atom,Box_dim,replicate)> %  replicate_atom.m This replicates the atom struct and the orthogonal box dimensions
% # <replace_atom.html replace_atom(new_atom,prev_atom,molid_index)> % This function replaces molid's in an atom struct with a new (single molid) atom struct by placing the latters COM in the formers place
% # <solvate_atom.html solvate_atom(limits,density,r,maxsol,solute_atom,varargin)> % This function generates a certain region defined by <limits> with a solvent structure of density <density>
% # <substitute_atom.html substitute_atom(atom,Box_dim,NumOctSubst,O1,O2,minO2O2_dist,varargin)> % This scripts performs isomorphous substitution, by replacing some O1->O2 atomtypes and optionally T1->T2 atomtypes

%% Slice out a region of the box 
% # <slice_atom.html slice_atom(atom,limits,invert)> % This function checks if the coordinates for each time record in XYZ_data is within the specified limits, and if not sets the x,y,z to nan,nan,nan.
% # <slice_molid.html slice_molid(atom,limits,invert)> % This function checks if the coordinates is within the specified limits, and if not sets the x,y,z to nan,nan,nan.
% # <slice_triclinic_atom.html slice_triclinic_atom(atom,limits,invert)>

%% Translate or rotate functions
% # <bend_atom.html bend_atom(atom,Box_dim,Radii)> % This simple function tries to bend an atom struct
% # <center_atom.html center_atom(atom,Box_dim,resname,dim)> % This function centers the atom with respect to the resname molecule
% # <condense_atom.html condense_atom(atom,Box_dim,s)> % This function tries to minimize the box size and remove gaps between molids?
% # <molid_rotate.html molid_rotate(atom,Box_dim,MolID,rotate_dim)> % This function rotate the atom randomly
% # <molid_translate.html molid_translate(atom,trans_vec,MolID)> % This translates a certain molid
% # <place_atom.html place_atom(atom,position)> % This function places the atom struct according to the position vector called position, trying to use the COM of the molecule
% # <position_molid.html position_molid(atom,position_vec,MolID)> % This function movies a molid (COM)> % to a certain position
% # <rotate_atom.html rotate_atom(atom,Box_dim,alfa,beta,gamma)> % This function rotate the atom randomly
% # <translate_atom.html translate_atom(atom,trans_vec,Resname)> % translate_atom.m - This translates the resname by a vector
% # <translate_molid.html translate_molid(atom,trans_vec,molid)> %  translate_molid.m - This translates the molid by a vector

%% Make triclinic/orthogonal box
% # <frac2atom.html frac2atom(atom,Box_dim,angleparam,angletype)> % This function transforms fractional coordinates to cartesian
% # <orto_atom.html orto_atom(atom,Box_dim)> % This transforms a triclinic atom struct to an orthogonal atom struct. Box_dim must look like [lx ly lz 0 0 xy 0 xz yz]
% # <triclinic_atom.html triclinic_atom(atom,Box_dim,angleparam,angletype)> %  triclinic_atom.m - This transforms an orthogonal atom struct to a triclinic with the angles alfa, beta, gamma or tilt factors xy, xz, yz

%% Wrap/unwrap functions
% # <wrap_atom.html wrap_atom(atom,Box_dim)> % This wraps the atoms into the orthogonal box
% # <Wrap_Coord_func.html Wrap_Coord_func(XYZ_data,Box_dim)> % This is an old function that wraps atoms 'sticking out' back into the box
% # <wrap_molid.html wrap_molid(atom,Box_dim)> % This function wraps the atom struct into the box
% # <unwrap_atom.html unwrap_atom(atom,Box_dim,dim)> % This function unwraps the atom struct along the dimension dim

%% Keep/remove functions
% # <keep_atom.html keep_atom(atom,resname)> % keep_atom.m - This removes all but resname
% # <keep_resname.html keep_resname(atom,resnames)> % keep_resname.m - This removes all but the resnames
% # <remove_molid.html remove_molid(atom,MolID)> %  remove_molid.m - This removes residue with molid MolID = [1 2 3 .....]
% # <remove_occypancy_atom.html remove_occypancy_atom(atom)> % This function removes all succeding particles in the atom struct that has identical coordinates to a preceding particle
% # <remove_residues.html remove_residues(atom,resnames,lo,hi,dim)> % This function section is used to remove residues in the simulation box between limits lo and hi
% # <remove_resname.html remove_resname(atom,resnames)> % This function removes residue with molid MolID, resnames = {'SOL' 'Protein'}
% # <remove_SOL.html remove_SOL(atom,atomname,lo,hi,dim)> %  This section is used to remove residues in the simulation box between limits lo and hi
% # <remove_type.html remove_type(atom,typescell)> % This function removes atomtypes with types as in typescell = {'OW' 'HW1' 'HW2'}

%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se