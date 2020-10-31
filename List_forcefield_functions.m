%% List of topology functions
%
%% Version
% 2.08
%
%% Custom toplogy tools
%% Clayff, with atomtypes by MHolmboe
% # <charge_atom.html charge_atom(atom,Box_dim,ffname,watermodel,varargin)> % This function tries to charge the atom accorind to clayff or interface ff
% # <charge_clayff_2004_atom.html charge_clayff_2004_atom(atom,Box_dim,varargin)> % Sets the charge for the original Clayff atomtypes for the original Clayff ff from the Cygan et al., 2004 paper
% # <charge_clayff_atom.html charge_clayff_atom(atom,Box_dim,varargin)> % Sets the charge for Clayff atomtypes 
% # <charge_opls_go_atom.html charge_opls_go_atom(atom,Box_dim,varargin)> % Sets the charge for some specific OPLS atomtypes
% # <check_clayff_2004_charge.html check_clayff_2004_charge(atom)> % check_clayff_2004_charge.m - This checks the charge of the original Clayff atomtypes by Mholmboe
% # <check_clayff_charge.html check_clayff_charge(atom)> % check_clayff_charge.m - This checks the charge of the Clayff atomtypes by Mholmboe
% # <check_clayff_H2Odens.html check_clayff_H2Odens(atom,Box_dim)> % Check the approx. water density for a clayff system
% # <check_H2Odens.html check_H2Odens(atom,Box_dim)> % Computes the water density
% # <clayff_2004_atom.html clayff_2004_atom(atom,Box_dim,varargin)> % Assigns the original Clayff atom types by Cygan et al., 2004. Can also 'heal' edges 
% # <clayff_2004_param.html clayff_2004_param(Atom_label,varargin)> % Holds the ion and the original Clayff atomtype parameters
% # <clayff_atom_old.m  clayff_atom_old(atom,Box_dim,varargin)> % Assigns the Clayff atom types by MHolmboe. Can also 'heal' edges. This is an older version... 
% # <clayff_atom.html clayff_atom(atom,Box_dim,varargin)> % Assigns the Clayff atom types by MHolmboe. Can also 'heal' edges 
% # <clayff_param.html clayff_param(Atom_label,varargin)> % Holds the ion and Clayff atomtype parameters
% # <clayffmod_atom.html clayffmod_atom(atom,Box_dim,varargin)> % Assigns the modififed Clayff atom types. Can also 'heal' edges 
% # <tweak_charge_atom.html tweak_charge_atom(atom)> % This function tries to tweak the charge of the atom struct in case of rounding errors

%% INTERFACE from Heinz 2005, 2013, with atomtypes by MHolmboe
% # <charge_atom.html charge_atom(atom,Box_dim,ffname,watermodel,varargin)> % This function tries to charge the atom accorind to clayff or interface ff
% # <interface_atom.html interface_atom(atom,Box_dim,varargin)> % This function tries to assign all atoms according to the interface atom types (with modified atom names by MHolmboe), with some modifications for edges...
% # <interface_param.html interface_param(Atom_label,water_model)> %  This holds the extended INTERFACE ff parameters
% # <interface15_atom.html interface15_atom(atom,Box_dim,varargin)> % This function tries to assign all atoms according to the interface1.5 atom types (with modified atom names by MHolmboe), with some modifications for edges...
% # <interface15_param.html interface15_param(Atom_label,water_model)> %  This holds the extended INTERFACE1.5 ff parameters
% # <charge_interface_atom.html charge_interface_atom(atom,Box_dim,varargin)> % Sets the charge for Interface atomtypes 
% # <charge_interface15_atom.html charge_interface15_atom(atom,Box_dim,varargin)> % Sets the charge for Interface 1.5 atomtypes
% # <check_interface_charge.html check_interface_charge(atom)> %  This checks the charge of the INTERFACE atomtypes by Mholmboe
% # <check_interface15_charge.html check_interface15_charge(atom)> %  This checks the charge of the INTERFACE 1.5 atomtypes by Mholmboe
% # <tweak_charge_atom.html tweak_charge_atom(atom)> % This function tries to tweak the charge of the atom struct in case of rounding errors

%% Graphene oxide modeled with OPLS/aa stuff
% # <opls_go_atom.html opls_go_atom(atom,Box_dim,rmin,rlarge)> % This function tries to smear out the charge at around -OH and epoxides in GO
% # <oplsaa_go_param.html oplsaa_go_param(Atom_label,water_model)> % This custom function holds the extended oplsaa_aa ff for graphite oxide
% # <charge_opls_go_atom.html charge_opls_go_atom(atom,Box_dim,varargin)> % Sets the charge for some specific OPLS atomtypes

%% Writing of topology files for Lammps (Clayff) and Gromacs (Clayff/INTERFACE)
% # <write_atom_all.html write_atom_all(atom,Box_dim,filename_out,varargin)> % This function tries to write various files for you. Works best for systems designed for Clayff...
% # <write_atom_itp.html write_atom_itp(atom,Box_dim,filename_out,varargin)> % This script creates and prints a gromacs .itp file. Works best for clayff or interface ff with spc, spce or tip3p
% # <write_atom_lmp.html write_atom_lmp(atom,Box_dim,filename_out,varargin)> % This script creates and prints a lammps data file (.lj). Works best for Clayff(_2004) systems
% # <write_atom_oplsaa_go_itp.html write_atom_oplsaa_go_itp(atom,Box_dim,filename_out,varargin)> % This custom made script creates and prints a gromacs .itp file for 
% # <write_atom_psf.html write_atom_psf(atom,Box_dim,filename_out,varargin)> % This function writes an .psf file from the atom struct
