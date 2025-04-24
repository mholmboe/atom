%% Import a unit cell structure
atom=import_atom('1xPyro_LeeGuggenheim_1981_centered.pdb'); % Imports a pyrophillite (montmorillonite) unit cell structure file into matlabs variable space
atom=replicate_atom(atom,Box_dim,[6 4 1]); % Replicate the structure by 6x4 into a clay layer
atom=substitute_atom(atom,Box_dim,6*4*2/3,'Al','Mgo',5.5) % Perform octahedral substitutions, 5.5 is the min dist between octahedral Mg's.
[atom.resname]=deal({'MMT'}); % Rename the residue name from PYR to MMT
pause(3)

%% Write a CLAYFF (Cygan, 2004) .gro and .itp file
atom=clayff_2004_atom(atom,Box_dim,'clayff_2004','spc'); % Assign the clayff atom types to the atomstruct
% atom=tweak_charge_atom(atom); % This function tries to tweak the charge if it is not exactly an integer. Passing an additional argument such as {'Si'} only tweaks the charge of the atomtype Si 
write_atom_itp(atom,Box_dim,'MMT_clayff_2004_2.itp',1.25,1.25,'clayff','spc');
write_atom_gro(atom,Box_dim,'MMT_clayff_2004_2'); % Print the clay sheet to a .gro file

%% Write a modified CLAYFF .gro and .itp file (the CLAYFF vesion with other/custom atomnames)
atom=clayff_atom(atom,Box_dim,'clayff','spc'); % Assign the clayff atom types to the atomstruct
% atom=tweak_charge_atom(atom); % This function tries to tweak the charge if it is not exactly an integer. Passing an additional argument such as {'Si'} only tweaks the charge of the atomtype Si 
write_atom_itp(atom,Box_dim,'MMT_clayff_2.itp',1.25,1.25,'clayff','spc');
write_atom_gro(atom,Box_dim,'MMT_clayff_2'); % Print the clay sheet to a .gro file

%% Write a interface (Heinz, 2005) .gro and .itp file
atom=interface_atom(atom,Box_dim,'interface') % Assign the interface atom types to the atomstruct
atom=tweak_charge_atom(atom); % This function tries to tweak the charge if it is not exactly an integer. Passing an additional argument such as {'Si'} only tweaks the charge of the atomtype Si 
write_atom_itp(atom,Box_dim,'MMT_interface_2.itp',1.25,2.25,'interface','spc');
write_atom_gro(atom,Box_dim,'MMT_interface_2'); % Print the clay sheet to a .gro file

% Is the charges correct?
% Before starting the simulation, make sure the topol.top file is correct,
% and note that you must change the moleculename in each .itp file, i.e.
% MMT should be changed to MMT_1 or MMT_2
