%% List of building and simulation cell manipulation functions
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%

%% Specific atom struct functions
% # <add_H_atom.html add_H_atom(atom,Box_dim,ind)> % This function protonates one or two sites in the atom struct 
% # <add2atom.html add2atom(XYZ_labels,XYZ_data,varargin)> % Append XYZ atomtype labels and XYZ data to an existing atom struct.
% # <adjust_H_atom.html adjust_H_atom(atom,Box_dim)> % Adjust hydrogen atoms in the atom struct.
% # <cat_atom.html cat_atom(atom_1,atom_2)> % Concatenate two atom structs.
% # <closest_atom.html closest_atom(atom,Box_dim,ref_atom)> % Find the closest atom to the reference atom.
% # <copy_atom.html copy_atom(atom,atomtype,new_atomtype,new_resname,trans_vec,varargin)> % Copy and translate atoms in the atom struct.
% # <create_atom.html create_atom(type,resname,limits,nmax,varargin)> % Create new atoms, useful for adding ions to a system.
% # <create_grid_atom.html create_grid_atom(atom_label,nM,limits,dim,varargin)> % Put ions on a grid plane and add them to an atom struct.
% # <duplicate_atom.html duplicate_atom(atom,molID)> % Duplicate residue with molid MolID.
% # <fuse_atom.html fuse_atom(atom,Box_dim,varargin)> % Fuse all sites within a certain cutoff distance.
% # <heal_atom.html heal_atom(atom,Box_dim,ind,varargin)> % Heal sites in the atom struct by adding a certain atom type.
% # <ionize_atom.html ionize_atom(type,resname,limits,nmax,varargin)> % Add ions within a certain region defined by limits.
% # <insert_atom.html insert_atom(atom,new_atom,position)> % Insert a new atom at the specified position.
% # <merge_atom.html merge_atom(atom1,Box1,atom2,type,Atom_label,r)> % Merge atom structs based on distance criteria.
% # <molid_rotate.html molid_rotate(atom,Box_dim,MolID,rotate_dim)> % Rotate the atom struct based on MolID.
% # <molid_translate.html molid_translate(atom,trans_vec,MolID)> % Translate a specific molecule ID.
% # <noupdate_atom.html noupdate_atom(atom)> % Prevent updating of certain properties in the atom struct.
% # <occupancy_atom.html occupancy_atom(atom,Box_dim)> % Calculate occupancy of atoms within the box dimensions.
% # <overwrite_atom.html overwrite_atom(In_atom,atomtype,resname)> % Overwrite atom struct information with new data.
% # <place_atom.html place_atom(atom,position)> % Place the atom struct at the specified position.
% # <position_molid.html position_molid(atom,position_vec,MolID)> % Move a molecule ID to a certain position.
% # <protonate_atom.html protonate_atom(atom,Box_dim,varargin)> % Protonate specified sites in the atom struct.
% # <remove_molid.html remove_molid(atom,MolID)> % Remove residue with a specific molecule ID.
% # <remove_occypancy_atom.html remove_occypancy_atom(atom)> % Remove particles with identical coordinates to preceding ones.
% # <remove_residues.html remove_residues(atom,resnames,lo,hi,dim)> % Remove residues between specified limits in the simulation box.
% # <remove_resname.html remove_resname(atom,resnames)> % Remove residues with specified names.
% # <remove_SOL.html remove_SOL(atom,atomname,lo,hi,dim)> % Remove solvent residues between specified limits.
% # <remove_type.html remove_type(atom,typescell)> % Remove atom types specified in typescell.
% # <rename_atom.html rename_atom(atom,old_name,new_name)> % Rename atoms in the atom struct.
% # <rename_type.html rename_type(atom,atomtype,new_atomtype,varargin)> % Rename atom types in the atom struct.
% # <replicate_atom.html replicate_atom(atom,Box_dim,replicate)> % Replicate the atom struct along orthogonal dimensions.
% # <replace_atom.html replace_atom(new_atom,prev_atom,molid_index)> % Replace molecule ID in an atom struct with a new atom struct.
% # <resname_atom.html resname_atom(atom)> % Guess residue names for all atom types.
% # <slice_atom.html slice_atom(atom,limits,invert)> % Slice the atom struct within specified limits.
% # <slice_box.html slice_box(atom,Box_dim,limits)> % Slice a simulation box within given limits.
% # <slice_molid.html slice_molid(atom,limits,invert)> % Slice molecules within specified limits.
% # <slice_triclinic_atom.html slice_triclinic_atom(atom,limits,invert)> % Slice a triclinic atom struct within limits.
% # <solvate_atom.html solvate_atom(limits,density,r,maxsol,solute_atom,varargin)> % Generate a solvent structure within specified limits.
% # <spc2tip4p.html spc2tip4p(atom)> % Convert SPC water molecules to TIP4P model.
% # <spc2tip5p.html spc2tip5p(atom)> % Convert SPC water molecules to TIP5P model.
% # <spce2tip4p.html spce2tip4p(atom)> % Convert SPC/E water molecules to TIP4P model.
% # <sphere_atom.html sphere_atom(atom,Box_dim,center,radius)> % Create a spherical region of atoms.
% # <substitute_atom.html substitute_atom(atom,Box_dim,NumOctSubst,O1,O2,minO2O2_dist,varargin)> % Perform isomorphous substitution of atoms.
% # <substitute_NonCentroSymm_atom.html substitute_NonCentroSymm_atom(atom,Box_dim,replace_type,varargin)> % Substitute non-centrosymmetric atoms.
% # <tile_atom.html tile_atom(atom,scale_vec,Box_dim,Resname)> % Tile the atom struct in a specific direction.
% # <tip3p2tip4p.html tip3p2tip4p(atom)> % Convert TIP3P water molecules to TIP4P model.
% # <translate_atom.html translate_atom(atom,trans_vec,Resname)> % Translate a residue by a specified vector.
% # <translate_molid.html translate_molid(atom,trans_vec,molid)> % Translate a molecule ID by a specified vector.
% # <tube_atom.html tube_atom(atom,scale_vec,Box_dim,Resname)> % Create a nanotube structure from the atom struct.
% # <update_atom.html update_atom(atom)> % Update molecule and atom indices in the atom struct.

%% Add/create/replicate/overwrite atoms
% # <add_H_atom.html add_H_atom(atom,Box_dim,ind)> % This function protonates one or two sites in the atom struct 
% # <copy_atom.html copy_atom(atom,atomtype,new_atomtype,new_resname,trans_vec,varargin)> % Copy and translate atoms in the atom struct.
% # <copy_type.html copy_type(atom,atomtype,new_atomtype,new_resname,trans_vec,varargin)> % Copy and translate types in the atom struct.
% # <create_atom.html create_atom(type,resname,limits,nmax,varargin)> % Create new atoms, useful for adding ions to a system.
% # <create_grid_atom.html create_grid_atom(atom_label,nM,limits,dim,varargin)> % Put ions on a grid plane and add them to an atom struct.
% # <duplicate_atom.html duplicate_atom(atom,molID)> % Duplicate residue with molid MolID.
% # <fuse_atom.html fuse_atom(atom,Box_dim,varargin)> % Fuse all sites within a certain cutoff distance.
% # <heal_atom.html heal_atom(atom,Box_dim,ind,varargin)> % Heal sites in the atom struct by adding a certain atom type.
% # <ionize_atom.html ionize_atom(type,resname,limits,nmax,varargin)> % Add ions within a certain region defined by limits.
% # <merge_atom.html merge_atom(atom1,Box1,atom2,type,Atom_label,r)> % Merge atom structs based on distance criteria.
% # <overwrite_atom.html overwrite_atom(In_atom,atomtype,resname)> % Overwrite atom struct information with new data.
% # <protonate_atom.html protonate_atom(atom,Box_dim,varargin)> % Protonate specified sites in the atom struct.
% # <replicate_atom.html replicate_atom(atom,Box_dim,replicate)> % Replicate the atom struct along orthogonal dimensions.
% # <replace_atom.html replace_atom(new_atom,prev_atom,molid_index)> % Replace molecule ID in an atom struct with a new atom struct.
% # <solvate_atom.html solvate_atom(limits,density,r,maxsol,solute_atom,varargin)> % Generate a solvent structure within specified limits.
% # <tile_atom.html tile_atom(atom,scale_vec,Box_dim,Resname)> % Tile the atom struct in a specific direction.
% # <translate_atom.html translate_atom(atom,trans_vec,Resname)> % Translate a residue by a specified vector.
% # <translate_molid.html translate_molid(atom,trans_vec,molid)> % Translate a molecule ID by a specified vector.
% # <tube_atom.html tube_atom(atom,scale_vec,Box_dim,Resname)> % Create a nanotube structure from the atom struct.

%% Slice out a region of the box
% # <slice_atom.html slice_atom(atom,limits,invert)> % Slice the atom struct within specified limits.
% # <slice_box.html slice_box(atom,Box_dim,limits)> % Slice a simulation box within given limits.
% # <slice_molid.html slice_molid(atom,limits,invert)> % Slice molecules within specified limits.
% # <slice_triclinic_atom.html slice_triclinic_atom(atom,limits,invert)> % Slice a triclinic atom struct within limits.

%% Translate or rotate functions
% # <bend_atom.html bend_atom(atom,Box_dim,Radii)> % Bend an atom struct.
% # <center_atom.html center_atom(atom,Box_dim,resname,dim)> % Center the atom with respect to the resname molecule.
% # <condense_atom.html condense_atom(atom,Box_dim,s)> % Minimize the box size and remove gaps between molids.
% # <molid_rotate.html molid_rotate(atom,Box_dim,MolID,rotate_dim)> % Rotate the atom struct based on MolID.
% # <molid_translate.html molid_translate(atom,trans_vec,MolID)> % Translate a specific molecule ID.
% # <place_atom.html place_atom(atom,position)> % Place the atom struct at the specified position.
% # <position_molid.html position_molid(atom,position_vec,MolID)> % Move a molecule ID to a certain position.
% # <translate_atom.html translate_atom(atom,trans_vec,Resname)> % Translate a residue by a specified vector.
% # <translate_molid.html translate_molid(atom,trans_vec,molid)> % Translate a molecule ID by a specified vector.

%% Make triclinic/orthogonal box
% # <triclinic_atom.html triclinic_atom(atom,Box_dim,angleparam,angletype)> % Transform an orthogonal atom struct to a triclinic one.

%% Wrap/unwrap functions
% # <unwrap_atom.html unwrap_atom(atom,Box_dim,dim)> % Unwrap the atom struct along the specified dimension.

%% Keep/remove functions
% # <remove_molid.html remove_molid(atom,MolID)> % Remove residue with a specific molecule ID.
% # <remove_occypancy_atom.html remove_occypancy_atom(atom)> % Remove particles with identical coordinates to preceding ones.
% # <remove_residues.html remove_residues(atom,resnames,lo,hi,dim)> % Remove residues between specified limits in the simulation box.
% # <remove_resname.html remove_resname(atom,resnames)> % Remove residues with specified names.
% # <remove_SOL.html remove_SOL(atom,atomname,lo,hi,dim)> % Remove solvent residues between specified limits.
% # <remove_type.html remove_type(atom,typescell)> % Remove atom types specified in typescell.
% # <rename_atom.html rename_atom(atom,old_name,new_name)> % Rename atoms in the atom struct.
% # <rename_type.html rename_type(atom,atomtype,new_atomtype,varargin)> % Rename atom types in the atom struct.

%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se