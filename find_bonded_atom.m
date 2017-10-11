%% find_bonded_atom.m
% * This function does a cross check of the bond matrix
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se


%% Examples
% * ind_12 = find_bonded_atom(atom,bond_matrix,Atom_label1,Atom_label2)

function ind_12 = find_bonded_atom(atom,bond_matrix,Atom_label1,Atom_label2)
%% 

ind_1=find(ismember([atom.type],Atom_label1));
ind_2=find(ismember([atom.type],Atom_label2));
[row,col]=find(bond_matrix(ind_2,:));
ind_12=intersect(col,ind_1);