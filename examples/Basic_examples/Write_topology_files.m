%% Examples demonstrating how to write topology files
% (For a full list of functions that deal with forcefield dependent stuff,
% go to 
% <List_forcefield_functions.html List_forcefield_functions>

%%
% In this example we will try to write so-called molecular topology files 
% in the .itp format (Gromacs) and in the .psf format (NAMD2). These 
% topology files carries the bonding and angle information needed. Note 
% that there is currently no support for dihedral/torsion interactions.
% Each format has its own function, like
% <import_atom_itp.html import_atom_itp> and 
% <import_atom_psf.html import_atom_psf>. Note also that both the original 
% CLAYFF forcefield is supported, as well as a modified CLAYFF, with other 
% atomnames, allowing for new atomtypes to be used (see the 
% <clayff_atom.html clayff_atom> function). The same goes for the Interface
% FF. 
%
% Important note about the Interface FF implementation - all bonds and 
% angles except the H-interactions are taken as the experimental bond 
% distances (x1.05)and angles as in Hendrik Heinz 2005 paper. In other
% words they are not set to fixed values as in the Heinz et al., 2013.
%
% These topology functions can be invoked by issuing something like:
% write_atom_itp(atom,Box_dim,'filename.itp',rmin,rmax,forcefieldname,watermodel);
% where rmin is the max cutoff for bonded H's and rmax the max cutoff for
% all other M-O bonds. Note the watermodel string (example 'spc') is not
% really needed and may be removed in later versions.
 
%% First set some convenient matlab settings
format compact; set(gcf,'Visible','on');

%% Import a unit cell structure
% First let's import a clay unit cell structure file into matlabs variable 
% space, with the option of replicating it to a larger and proper clay layer. 
atom=import_atom('Pyrophyllite.pdb'); 
atom=replicate_atom(atom,Box_dim,[6 4 1]); 

%% Assign the forcefield atomtypes
% If we have not already done so, we should assign the forcefield specific 
% atomtypes to all atomnames. Here we will do it for original CLAYFF as
% well as the modified (with atomnames made up by me..), and the
% corresponding Interface FF v1.5.
atom_clayff_2004 = clayff_2004_atom(atom,Box_dim,'clayff');
atom_clayff = clayff_atom(atom,Box_dim); % Modifed CLAYFF
atom_interface15 = interface15_atom(atom,Box_dim,'Interface');
atom_interface = interface_atom(atom,Box_dim,'Interface'); % Modified Interface FF

%% Write a CLAYFF (Cygan, 2004).itp/.psf file
write_atom_itp(atom_clayff_2004,Box_dim,'pyro_clayff_2004.itp',1.25,1.25,'clayff_2004','spc');
write_atom_psf(atom_clayff_2004,Box_dim,'pyro_clayff_2004.psf',1.25,1.25,'clayff_2004','spc');

%% Write a modified CLAYFF .itp/.psf file with modified atomnames
write_atom_itp(atom_clayff,Box_dim,'pyro_clayff.itp',1.25,1.25,'clayff','spc');
write_atom_psf(atom_clayff,Box_dim,'pyro_clayff.psf',1.25,1.25,'clayff','spc');

%% Write a Interface FF (Heinz, 2005) .itp/.psf file
write_atom_itp(atom_interface15,Box_dim,'pyro_interface15.itp',1.25,2.25,'interface15','spc');
write_atom_psf(atom_interface15,Box_dim,'pyro_interface15.psf',1.25,2.25,'interface15','spc');

%% Write a Interface FF .itp/.psf file with modified atomnames
write_atom_itp(atom_interface,Box_dim,'pyro_interface.itp',1.25,2.25,'interface','spc');
write_atom_psf(atom_interface,Box_dim,'pyro_interface.psf',1.25,2.25,'interface','spc');


