%% protonate_matrix.m
% * This function protonates the sites in the atom struct given by the index vector ind by adding a H's to a new H atom struct. It does so by placing the H opposite to the mean position of all neughbours within 2.5 Ångrtröm of the site to be protonated
% * Not finished yet
% * Please report bugs to michael.holmboe@umu.se


%% Examples
% * Hatom = protonate_atom(atom,Box_dim,ind)
% * Hatom = protonate_atom(atom,Box_dim,ind,{'He'}) % {'He'} can be used to change the default atomtype H to He
% * Hatom = protonate_atom(atom,Box_dim,ind,{'He'},rcut) % rcut can be used to change the default cutoff 2.5 Ångström
% * Hatom = protonate_atom(atom,Box_dim,ind,{'He'},rcut,'minus') % 'minus' or default 'plus' denotes the tilt direction of the added H in the Z-direction

function H_atom = protonate_atom(atom,Box_dim,ind,varargin)