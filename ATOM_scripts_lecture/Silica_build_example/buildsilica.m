%% Copy your original .pdb to a file called silica.pdb
%% Import a the silica structure
atom=import_atom('silica.pdb'); % Imports a montmorillonite unit cell structure file into matlabs variable space
atom=element_atom(atom);

%% Write a clayff (Cygan, 2004) .gro and .itp file
atom=clayff_atom(atom,Box_dim,'clayff_2004') % Assign the clayff atom types to the atomstruct
write_atom_itp(atom,Box_dim,'SIL_clayff_2004.itp',1.25,1.25,'clayff','spc');
write_atom_pdb(atom,Box_dim,'SIL_clayff_2004.pdb'); % Print the clay sheet to a .gro file

%% Write a modified CLAYFF .gro and .itp file (the CLAYFF vesion with other/custom atomnames)
atom=clayff_atom(atom,Box_dim,'clayff','spc'); % Assign the clayff atom types to the atomstruct
% atom=tweak_charge_atom(atom); % This function tries to tweak the charge if it is not exactly an integer. Passing an additional argument such as {'Si'} only tweaks the charge of the atomtype Si 
write_atom_itp(atom,Box_dim,'SIL_clayff.itp',1.25,1.25,'clayff','spc');
write_atom_pdb(atom,Box_dim,'SIL_clayff.pdb'); % Print the clay sheet to a .gro file

%% Write a interface (Heinz, 2005) .gro and .itp file
atom=interface_atom(atom,Box_dim,'interface','tip3p');
% atom=tweak_charge_atom(atom); % This function tries to tweak the charge if it is not exactly an integer. Passing an additional argument such as {'Si'} only tweaks the charge of the atomtype Si 
write_atom_itp(atom,Box_dim,'SIL_interface.itp',1.25,2.25,'interface','spc');
write_atom_pdb(atom,Box_dim,'SIL_interface.pdb'); % Print the clay sheet to a .gro file

%% Write a interface (Heinz, 2013) .gro and .itp file
atom=interface15_atom(atom,Box_dim,'interace15','tip3p',[1],'SILICA');
% atom=tweak_charge_atom(atom); % This function tries to tweak the charge if it is not exactly an integer. Passing an additional argument such as {'Si'} only tweaks the charge of the atomtype Si 
write_atom_itp(atom,Box_dim,'SIL_interface15.itp',1.25,2.25,'interface15','spc','SILICA');
write_atom_pdb(atom,Box_dim,'SIL_interface15.pdb'); % Print the clay sheet to a .gro file

%% Is the charge what you would expect?