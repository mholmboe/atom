
[atom,traj] = import_atom_xtc(filenameconf,filenamextc) %  import_atom_xtc.m - This imports a structure file and a xtc file

atom = add2atom(XYZ_labels,XYZ_data,resname,in_atom) %  add2atom.m - This function appends so-called XYZ atomtype labels and XYZ data to the atom struct

atom = analyze_atom(atom,Box_dim,max_short_dist,max_long_dist) %  analyze_atom.m - This analyzes variofus things of the MMT atom struct

atom = ave_atom(atom) %  ave_atom.m - This calculates the mean of the atom coordinates

atom = bond_angle_atom(atom,Box_dim,max_short_dist,max_long_dist,varargin) %  bond_angle_atom.m - This tries to find all bonds and angles of the atom struct 'more' is an optional varargin argument

atom = center_atom(atom,Box_dim,resname,dim) %  center_atom.m - This centers the atom with respect to the resname molecule

atom = charge_atom(atom,Box_dim,ffname,watermodel,varargin) %  check_clayff_charge.m - This tries to charge the atom accorind to clayff or interface ff

atom = check_clayff_charge(atom) %  check_clayff_charge.m - This checks the charge of the Clayff atomtypes

atom = check_INTERFACE_charge(atom) %  check_INTERFACE_charge.m - This checks the charge of the INTERFACE atomtypes

atom = COM_atom(atom,MolID) %  COM_atom.m - This calculates the COM for certain elements

atom = COM_molid(atom,MolID) %  COM_molid.m - This calculates the COM for certain elements

atom = concatenate_atom(atom_1,atom_2) %  concatenate_atom.m - This concatenats atom sections.

atom = condense_atom(atom,Box_dim,s) %  condense_atom.m - This tries to minimize the box size and remove gaps between molids?

atom = duplicate_atom(atom,molID) %  duplicate_atom.m - This duplicates residue with molid MolID

atom = grid2atom(atom_label,nM,limits,dim,varargin) %  grid2atom.m - This puts particles such as ions on a 2D grid (i.e. a plane) % and adds it to an atom struct

atom = import_atom_gro(filenamegro) %  import_atom_gro.m - This imports .gro files into the atom struct

atom = import_atom_xyz(filename) %  import_atom_xyz.m - This imports an .xyz file into the atom struct

atom = import_atom(filename) %  import_atom.m - This imports a .xyz/.gro/.pdb file and puts the data in a structure variable called atom

atom = import_traj(filenameconf,filenametraj) %  import_traj.m - This imports an strcture and an dcd, xtc, xyz or gro  trajectory file.

atom = import_gro_traj(filenametraj) %  import_gro_traj.m - This imports an strcture and an .gro trajectory file.

atom = import_xyz_traj(filenametraj) %  import_xyz_traj.m - This imports an strcture and an .xyz trajectory file.

atom = insert_atom(atom_in,limits,rotate,r,maxsol,solute_atom,varargin) % - This inserts a molecule from a structure file into a region defined by <limits> with a atom (molecule) % structure

atom = keep_atom(atom,resname) %  keep_atom.m - This removes all but resname

atom = keep_resname(atom,resnames) %  keep_resname.m - This removes all but the resnames

atom = median_atom(atom) %  median_atom.m - This calculates the median position of the atom struct

atom = molid_rotate(atom,Box_dim,MolID,rotate_dim) %  molid_rotate.m - This rotate the atom randomly

atom = molid_translate(atom,trans_vec,MolID) %  molid_translate.m - This translates a certain molid

atom = orto_atom(atom,Box_dim) %  orto_atom.m - This transforms a triclinic atom struct to an orthogonal atom struct. Box_dim must look like [lx ly lz 0 0 xy 0 xz yz]

atom = position_molid(atom,position_vec,MolID) %  position_molid.m - This movies a molid (COM) % to a certain position

atom = remove_molid(atom,MolID) %  remove_molid.m - This removes residue with molid MolID = [1 2 3 .....]

atom = remove_residues(atom,resnames,lo,hi,dim) %  remove_residues.m - This section is used to remove residues in the simulation box between limits lo and hi

atom = remove_resname(atom,resnames) %  remove_resname.m - This removes residue with molid MolID, resnames = {'SOL' 'Protein'}

atom = remove_SOL(atom,atomname,lo,hi,dim) %  This section is used to remove residues in the simulation box between limits lo and hi

atom = remove_type(atom,typescell) %  remove_type.m - This removes atomtypes with types as in typescell = {'OW' 'HW1' 'HW2'}

atom = remove_type(atom,typescell) %  remove_type.m - This removes atomtypes with types as in typescell = {'OW' 'HW1' 'HW2'}

atom = replicate_atom(atom,Box_dim,replicate) %  replicate_atom.m This replicates the atom struct and the orthogonal box dimensions

atom = rotate_atom(atom,Box_dim,alfa,beta,gamma) %  atom_rotate.m - This rotate the atom randomly

atom = scale_atom(atom,scale_vec,Box_dim,Resname) %  scale_atom.m - This scales the coordinates in the atom struct

atom = slice_atom(atom,limits,invert) %  slice_atom.m - This checks if the coordinates for each time record in XYZ_data is within the specified limits, and if not sets the x,y,z to nan,nan,nan.

atom = slice_molid(atom,limits,invert) %  slice_molid.m - This checks if the coordinates is within the specified limits, and if not sets the x,y,z to nan,nan,nan.

SOL = solvate_atom(limits,density,r,maxsol,solute_atom) %  solvate_atom.m - This SOLvates a certain region defined by limits with a water structure with density density. r and r-0.5 is the closest distance of Ow and Hw to the solute atoms

atom = sort_atom(atom) %  sort_atom.m - This section orders to atoms with respect to z

atom = translate_atom(atom,trans_vec,Resname) % translate_atom.m - This translates the resname by a vector

atom = translate_molid(atom,trans_vec,molid) %  translate_molid.m - This translates the molid by a vector

atom = triclinic_atom(atom,Box_dim,angleparam,angletype) %  triclinic_atom.m - This transforms an orthogonal atom struct to a triclinic with the angles alfa, beta, gamma or tilt factors xy, xz, yz

atom = triclinic_atom(atom,Box_dim,angleparam,angletype) %  triclinic_atom.m - This transforms an orthogonal atom struct to a triclinic with the angles alfa, beta, gamma or tilt factors xy, xz, yz

atom = unwrap_atom(atom,Box_dim,dim) %  unwrap_atom.m - This unwraps the atom struct along the dimension dim

atom = update_atom(atom) %  update_atom.m - This updates the molid index and the atoms index in the atom struct

atom = wrap_atom(atom,Box_dim) % wrap_atom.m - This wraps the atoms into the orthogonal box

Atom_prop = lmp_atom_style_full_func(fid,Atom_label,Charge,XYZ_labels,XYZ_data) %  lmp_atom_style_full_func.m - This creates and prints the Atoms section properties in the LAMMPS data file .lj file according to atom style full, without image flags

atom = rename_type(atom,atomtype,new_atomtype,varargin) %  rename_type.m - This renames atoms in the atom

AtomCoords_COM = COM_func(MolID,XYZ_data,Atom_label,XYZ_labels,Box_dim) %   COM_func.m This calculates the center of mass for water. Slow due to pbc...

clayff_param(Atom_label,water_model) %  clayff_param.m - This holds the extended Clayff ff parameters

composition_atom(atom) %  composition_atom.m - This looks at the composition of the atom struct

CONECT_atom(atom,Box_dim,short_r,long_r) %  CONECT_atom.m - This prints conect records for pdb files

dipole_vec = dipoles_atom(Elements,Box_dim) %  dipole_atom.m - This calculates the dipole vector of water. Similar to the COM_func

dist_matrix = dist_matrix_atom(atom,Box_dim) %   dist_matrix.m - This calculates the distance matrix from the atom struct

frame = import_atom_multiple_gro(filename,nFrames) %  import_atom_multiple_gro.m - This import multiple .gro files

gmx_make_ndx(groupname,ind) %  gmx_make_ndx.m - This helps you print custom gromacs .ndx files

H2Odens = check_clayff_H2Odens(atom,Box_dim) %  check_clayff_H2Odens.m - Check the approx. water density for a clayff system

Hist = hist_atom(atom,s) %  hist_atom.m - This is used to calculate density profiles in the Z-direction

ind_12 = find_bonded_atom(atom,bond_matrix,label1,label2) %  find_bonded_atom.m - This does a cross check of the bond matrix

interface_param(Atom_label,water_model) %  interface_param.m - This holds the extended INTERFACE ff parameters

new_atom = ions2atom(atom_label,nM,limits,dim) %  ions2atom.m - This puts ions on a grid plane and adds it to an atom struct

new_atom = copy_type(atom,atomtype,new_atomtype,new_resname,trans_vec,varargin) %  copy_atom.m - This copies and translates types in the atom

new_atom = frame2atom(atom,Elements,frame,Box_dim) %  frame2atom.m - This extracts a frame from a trajectory matrix to the atom struct

new_atom = merge_atom(atom1,Box1,atom2,Box2,type,Atom_label,r) %  merge_atom.m - This returns the second atom set with non-overlapping atoms with a distance r away from the atoms in the first atom set

new_atom = copy_atom(atom,atomtype,new_atomtype,new_resname,trans_vec) %  copy_atom.m - This copies and translates atoms in the atom

new_atom = overwrite_atom(In_atom,atomtype,resname) %  overwrite_atom.m - This overwrites the atom struct information with new information 

reorder_atom_gro(atom,atomlist,Box_dim,filename_out) %  reorder_atom_gro.m - This reorders the atoms in a .gro file

SOL_atom = sol2atom(atom,limits,nH2O,dist,CLAYFF_param,ExcludeEdge) %  sol2atom.m - This solvates a certain region with nH2O molecules

sol2xyz(atom,limits,nH2O,dist,CLAYFF_param,ExcludeEdge) %  sol2xyz.m - This solvates a certain region with nH2O molecules

sorted_molid=sort_molid(MolID) %  sort_molid.m - This sorts the molecular indexes in an ascending order

tip4p_atom = spc2tip4p(filename) %  spc2tip4p.m - This converts a .gro file with spc water to some tip4p water

tip4p_atom = spce2tip4p(filename) %  spce2tip4p.m - This converts a .gro file with spc water to some tip4p water

tip4p_atom = tip3p2tip4p(filename) %  tip3p2tip4p.m - This converts a .gro file with spc water to some tip4p water

vmd(atom,Box_dim) %  vmd.m - This plots the atom struct

write_gro_traj(atom,traj,Box_dim,filename_out) %  write_gro_traj.m - This writes a .gro trajectory 

write_xyz_traj(atom,traj,Box_dim,filename_out) %  write_xyz_traj.m - This writes a .xyz trajectory 

write_atom_gro(atom,Box_dim,filename_out) %  write_atom_gro.m - This writes a gro file. Does it also write velocities?

write_atom_mol2(atom,Bond_index,Box_dim,filename_out) %  write_atom_mol2.m - This writes a .mol2 file from the atom struct

write_atom_multiple_gro(atom,traj,filename_out) %  write_atom_multiple_gro.m - This writes a .gro trajectory

write_atom_pdb(atom,Box_dim,filename_out) %  write_atom_pdb.m - This writes an .pdb file from the atom struct using Gromacs

write_atom_xyz(atom,Box_dim,filename_out) %  write_atom_xyz.m - This writes an XYZ file from the atom struct



