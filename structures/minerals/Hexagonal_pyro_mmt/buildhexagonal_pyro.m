clear all;
format compact;
tic
%% Import_the hexagonal Pyrophyllite structure, cut in Avogadro for instance
filename='hex_n1570_5nm_pyro';
import_atom(strcat(filename,'.pdb'));

%% Assign the first atomtypes
% atom1=clayff_atom(atom,Box_dim,'clayff','spc');
atom1=clayff211_atom(atom,Box_dim);
%% Find B-type edge O-atoms
fb_atom = find_bonded_atom(atom1,Box_dim,'Oalh','H');
ind1=type1_ind; ind2=type2_ind;
%% Heal the total structure
% atom2=clayff_atom(atom1,Box_dim,'clayff','spc',[6:8]);
atom2=clayff211_atom(atom1,Box_dim,'clayff','spc',[6:8]);
%% Find the H's bonded to the B-type edge O-atoms, and remove some the get pH neutral charge
fb_atom = find_bonded_atom(atom2,Box_dim,ind1,'H');
rm_ind=[type2_ind(length(type2_ind)/2+1:2:end)];
atom2(rm_ind)=[];
%% Assign the all the final atomtypes
% atom3=clayff_atom(atom2,Box_dim,'clayff','spc');
atom3=clayff211_atom(atom2,Box_dim);
%% Print the final hexagonal structure
write_atom_gro(atom3,Box_dim,strcat(filename,'_clayff.gro'))
% vmd(atom3,Box_dim)

toc