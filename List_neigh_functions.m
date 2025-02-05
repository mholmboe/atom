%% List of neighbour analysis and distance matrix functions
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%

%% Neighbor/distance functions
% # <bond_atom.html bond_atom(atom,Box_dim,max_long_dist)> % This function tries to assign all bonds to a Bond_matrix and a Bond_index variable.
% # <bond_angle_atom.html bond_angle_atom(atom,Box_dim,varargin)> % This function tries to find all bonds and angles of the atom struct 'more' is an optional varargin argument.
% # <bond_angle_dihedral_atom.html bond_angle_dihedral_atom(atom,Box_dim,varargin)> % This function tries to find all bonds, angles, and dihedrals of the atom struct. Rmaxshort and Rmaxlong as well as 'more' is an optional varargin argument.
% # <bond_angle_type.html bond_angle_type(atom1,atom2,Box_dim,rmin,rmax,angle_limit,varargin)> % This function tries to find all bonds and angles of the atom types.
% # <bond_valence_atom.html bond_valence_atom(atom,Box_dim,varargin)> % This function tries to calculate the bond valence values according to the bond valence method.
% # <bond_valence_data.html bond_valence_data(ion1,ion2,R,varargin)> % This function fetches the data and matches it to the passed atom types used to calculate the bond valence values according to http://www.iucr.org/resources/data/datasets/bond-valence-parameters.
% # <cell_list_dist_matrix_atom.html cell_list_dist_matrix_atom(atom,Box_dim,varargin)> % This function calculates the distance matrix from the atom struct, using a cell list algorithm.
% # <cell_list_dist_matrix_atom1atom2.html cell_list_dist_matrix_atom1atom2(atom,Box_dim,varargin)> % This function calculates the distance matrix from two structs called atom1 and atom2, using a cell list algorithm.
% # <closest_atom.html closest_atom(atom1,atom2,Box_dim)> % This function returns the atom1 struct with the nMolId's in atom1 closest to the atom2 struct.
% # <dist_matrix_atom.html dist_matrix_atom(atom,Box_dim)> % This function calculates the distance matrix from the atom struct.
% # <dist_matrix_noPBC_atom.html dist_matrix_noPBC_atom(atom,Box_dim)> % This function calculates the distance matrix without applying periodic boundary conditions.
% # <dist_matrix_xyz.html dist_matrix_xyz(XYZ,Box_dim)> % This function calculates the distance matrix from XYZ coordinates.
% # <find_bonded_atom.html find_bonded_atom(atom,bond_matrix,label1,label2)> % This function performs a cross-check of the bond matrix.
% # <find_pair_atom.html find_pair_atom(atom,bond_matrix,pair1,pair2)> % This function finds and returns specific atom pairs from the bond matrix.
% # <list_bonds_atom.html list_bonds_atom(atom,bond_matrix)> % This function lists all bonds in the atom struct based on the bond matrix.
% # <neigh_atom.html neigh_atom(atom,Box_dim,rmax,varargin)> % This function checks which neighbors each atom has and outputs their information.
% # <neighbor_func.html neighbor_func(solute_index,XYZ_solute,XYZ_data,Box_dim,radius)> % This function scans XYZ data and checks which atoms are within a certain radius. It outputs the neighbour index.
% # <neighbor_atom.html neighbor_atom(atom,Box_dim,radius)> % This function checks the neighbors of each atom and returns their indices.
% # <rdf_atom.html rdf_atom(atom,Box_dim,varargin)> % This function calculates the radial distribution function and the coordination number, with optional Gaussian smoothing.
% # <cn_atom.html cn_atom(atom,Box_dim,rmax)> % This function calculates the coordination number of atoms within a specified radius.
% # <recalc_bond_atom.html recalc_bond_atom(atom,bond_matrix,varargin)> % This function recalculates bonds for the atom struct.
% # <remove_H2O.html remove_H2O(atom)> % This function removes water molecules (H2O) from the atom struct.
% # <remove_sametype_bond.html remove_sametype_bond(atom,bond_matrix)> % This function removes bonds between atoms of the same type.
% # <find_H2O.html find_H2O(atom)> % This function identifies and returns water molecules (H2O) within the atom struct.
% # <bond_matrix_atom.html bond_matrix_atom(atom,Box_dim)> % This function generates a bond matrix for the atom struct.

%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
