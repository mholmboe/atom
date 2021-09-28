%% List of neighbour analysis and distance matrix functions
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%

%% Neighbor/distance functions
% # <bond_atom.html bond_atom(atom,Box_dim,max_long_dist)> % This function tries to assign all bonds to a Bond_matrix and a Bond_index variable
% # <bond_angle_atom.html bond_angle_atom(atom,Box_dim,varargin)> % This function tries to find all bonds and angles of the atom struct 'more' is an optional varargin argument
% # <bond_angle_dihedral_atom.html bond_angle_dihedral_atom(atom,Box_dim,varargin)> % This function tries to find all bonds, angles and dihedrals of the atom struct. Rmaxshort and Rmaxlong as well as 'more' is an optional varargin argument
% # <bond_angle_type.html bond_angle_type(atom1,atom2,Box_dim,rmin,rmax,angle_limit,varargin)> % This tries to find all bonds and angles of the atom types
% # <bond_valence_atom.html bond_valence_atom(atom,Box_dim,varargin)> % This function tries to calculate the bond valence values according to the bond valence method
% # <bond_valence_data.html bond_valence_data(ion1,ion2,R,varargin)> % This function fetches the data and matches it to the passed atom types used to calculate the bond valence values according to http://www.iucr.org/resources/data/datasets/bond-valence-parameters
% # <cell_list_dist_matrix_atom.html cell_list_dist_matrix_atom(atom,Box_dim,varargin)> % This function calculates the distance matrix from the atom struct, using a cell list algorithm adapted from the Matlab MDtoolbox by Yasuhiro Matsunaga
% # <dist_matrix_atom.html dist_matrix_atom(atom,Box_dim)> % This calculates the distance matrix from the atom struct
% # <find_bonded_atom.html find_bonded_atom(atom,bond_matrix,label1,label2)> % This function does a cross check of the bond matrix
% # <neigh_atom.html neigh_atom(atom,Box_dim,rmax,varargin)> % This function checks which neighbors each atom has and ouputs their info
% # <neighbor_func.html neighbor_func(solute_index,XYZ_solute,XYZ_data,Box_dim,radius)> % This function scans xyz data and checks who is within a certain radius. It outputs neighbour index,
% # <rdf_atom.html rdf_atom(atom,Box_dim,varargin)> % This function calculates the radial distributtion function and the coordination number. Can also do Gaussion smoothing.
