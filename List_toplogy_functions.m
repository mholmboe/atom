%% List of topology functions

%% Custom toplogy tools
%% Clayff, with atomtypes by MHolmboe
clayff_atom(atom,Box_dim,varargin) % Assigns the Clayff atom types by MHolmboe. Can also 'heal' edges
clayffmod_atom(atom,Box_dim,varargin) % Assigns the modififed Clayff atom types. Can also 'heal' edges 
clayff_param(Atom_label,varargin) % Holds the ion and Clayff atomtype parameters
charge_clayff_atom(atom,Box_dim,varargin) % Sets the charge for Clayff atomtypes 
check_clayff_charge(atom) % This checks the charge of the Clayff atomtypes by Mholmboe
check_clayff_H2Odens(atom,Box_dim) % Specific function that calculates the approx. water density in heterogeneous systems
mass_atom_clayff(atom,varargin) % This function fetches the atom weight from the clayff and interface ff's
tweak_charge_atom(atom) % This function tries to tweak the charge of the atom struct in case of rounding errors

%% INTERFACE from Heinz 2005, 2013, with atomtypes by MHolmboe
interface_atom(atom,Box_dim,varargin) % This function tries to assign all atoms according to the interface atom types (with modified atom names by MHolmboe), with some modifications for edges...
interface_param(Atom_label,water_model) %  This holds the extended INTERFACE ff parameters
interface15_atom(atom,Box_dim,varargin) % This function tries to assign all atoms according to the interface1.5 atom types (with modified atom names by MHolmboe), with some modifications for edges...
interface15_param(Atom_label,water_model) %  This holds the extended INTERFACE1.5 ff parameters
check_interface_charge(atom) %  This checks the charge of the INTERFACE atomtypes by Mholmboe
check_interface15_charge(atom) %  This checks the charge of the INTERFACE 1.5 atomtypes by Mholmboe
charge_interface_atom(atom,Box_dim,varargin) % Sets the charge for Interface atomtypes 
charge_interface15_atom(atom,Box_dim,varargin) % Sets the charge for Interface 1.5 atomtypes
mass_atom_clayff(atom,varargin) % This function fetches the atom weight from the clayff and interface ff's
tweak_charge_atom(atom) % This function tries to tweak the charge of the atom struct in case of rounding errors

%% Graphene oxide modeled with OPLS/aa stuff
opls_go_atom(atom,Box_dim,rmin,rlarge) % This function tries to smear out the charge at around -OH and epoxides in GO
oplsaa_go_param(Atom_label,water_model) % This custom function holds the extended oplsaa_aa ff for graphite oxide
charge_opls_go_atom(atom,Box_dim,varargin) % Sets the charge for some specific OPLS atomtypes

%% Writing of topology files for Lammps (Clayff) and Gromacs (Clayff/INTERFACE)
write_atom(atom,Box_dim,filename_out,varargin) % This function tries to write various files for you. Works best for systems designed for Clayff...
write_atom_lmp(atom,Box_dim,filename,varargin) % This script creates and prints a lammps data file (.lj). Works best for Clayff systems
write_atom_itp(atom,Box_dim,filename,varargin) % This script creates and prints a gromacs .itp file. Works best for clayff or interface ff with spc, spce or tip3p
write_atom_oplsaa_go_itp(atom,Box_dim,filename,varargin) % This custom made script creates and prints a gromacs .itp file for 
write_atom_psf(atom,Box_dim,filename_out,varargin) % This function writes an .psf file from the atom struct
