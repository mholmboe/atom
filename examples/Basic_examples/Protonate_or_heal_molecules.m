%% Examples demonstrating how to protonate, heal, fuse or merge molecules
% (For a full list of functions that can add atoms/ions/molecules, go to 
% <List_build_functions.html List_build_functions> or look into 
% <List_general_functions.html List_general_functions>

%% First set some convenient matlab settings
format compact; set(gcf,'Visible','on');

%% Protonate a molecule with <protonate_atom.html protonate_atom>
% * This function protonates the sites in the atom struct given by the 
% index vector ind by adding H's to a new H atom struct. It does so by 
% placing a H opposite to the mean position of all neighbours within 2.5
% Ångström of the site to be protonated.
%
atom=import_atom('Pyrophyllite.pdb');
atom(strcmp([atom.type],'H'))=[]; % We need to remove the original H's so 
% we have some sites to protonate
%%
% *Examples*
Hatom = protonate_atom(atom,Box_dim) % Will only protonate single bonded O atoms
Hatom = protonate_atom(atom,Box_dim,ind) % Will protonate all atoms with indexes in ind
Hatom = protonate_atom(atom,Box_dim,ind,rmaxlong) % rcut can be used to change the default cutoff 2.5 Ångström
Hatom = protonate_atom(atom,Box_dim,ind,rmaxlong,{'He'}) % {'He'} can be used to change the default atomtype H to He
Hatom = protonate_atom(atom,Box_dim,ind,rmaxlong,{'He'},'minus') % 'minus' or default 'plus' denotes the tilt direction of the added H in the Z-direction

%% Heal a molecule with <heal_atom.html heal_atom>
% * This function heals sites in the atom struct given by the index
% vector ind, by adding a certain atomtype to a new atom struct called
% healed_atom. It does so by placing the new atom type opposite to the
% mean position of all neighbours within rcut [Å] of the healed site.
% * Note that you first need to find which sites that should be healed,
% using for instance the bond_valence_atom function, and then decide with
% what atomtypes the unsaturated sites should be healed with.
%
%%
% *Examples*
healed_atom = heal_atom(atom,Box_dim,[6 16 26 36]) % Will protonate all sites given by the ind array
healed_atom = heal_atom(atom,Box_dim,[6:10:960],3) % Will protonate all sites given by the ind array. 3 is the rcut
healed_atom = heal_atom(atom,Box_dim,[3 202 493],3,'He') % Will add 'He' atoms to all sites given by the ind array. 3 is the rcut
healed_atom = heal_atom(atom,Box_dim,[3 202 493],2.5,'N',1.87) % Will add 'N' atoms 1.87 Å away from all sites given by the ind array
healed_atom = heal_atom(atom,Box_dim,[3 202 493],2.5,'N','Shannon') % Same as above but will use bond distance from the revised Shannon radii set

%% Fuse two molecules with <fuse_atom.html fuse_atom>
% * Sometimes added atoms/sites overlap. This function tries to fuse all 
% sites within a certain rmax. 
%
% *Examples*
atom = fuse_atom(atom,Box_dim)
atom = fuse_atom(atom,Box_dim,0.85)

%% Merge two molecules with <merge_atom.html merge_atom>
% * This function returns the struct atom2 containing atoms
% with a distance r (1x1 or 1x2 array) away from the atoms in the atom1 
% struct. There is also a possibility to use a twin-range cutoff approach 
% (suitable for OH2), by setting r(2) to a smaller value than r(1).
% * You must decide if atom1w should be wrapped or not before running this
% function
% * Do not wrap atomw1 here, do it before calling this function since 
% Box1 does not always equal the full Box_dim, but rather a region in the 
% full Box_dim
%
%%
% *Examples*
atom2w = merge_atom(SOLUTE,Box_dim,SOLVENT,'index','C',1.4)
atom2w = merge_atom(SOLUTE,Box_dim,SOLVENT,'molid','Hw',[1.6 1.0])
