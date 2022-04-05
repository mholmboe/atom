%% sigma_vdw.m
% * This function fetches the sigma radius parameter, originally taken from the link below
% * from 'A cartography of the van der Waals territories'
% * Santiago Alvarez doi:10.1039/c3dt50599e
% * Z 61 and 84-88 are made up....
%
%% Version
% 2.11
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # sigma = sigma_vdw({'O'})
% # sigma = sigma_vdw('O')

function sigma = sigma_vdw(Atom_label)

rvdw = radius_vdw(Atom_label);

sigma = rvdw/(2^(1/6));


