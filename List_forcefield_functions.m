%% List of topology functions
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%

%% Custom topology tools

%% minff, with atomtypes by MHolmboe
% # <minff_atom.html minff_atom(atom,Box_dim,varargin)> % Assign minff atom types with edge healing.
% # <charge_minff_atom.html charge_minf_atom(atom,Box_dim,varargin)> % Set the charge for MINFF atom types.

%% Clayff, with atomtypes by MHolmboe
% # <charge_atom.html charge_atom(atom,Box_dim,ffname,watermodel,varargin)> % Charge the atom according to Clayff or Interface force fields.
% # <charge_clayff_2004_atom.html charge_clayff_2004_atom(atom,Box_dim,varargin)> % Set the charge for the original Clayff atom types from the Cygan et al., 2004 paper.
% # <charge_clayff_atom.html charge_clayff_atom(atom,Box_dim,varargin)> % Set the charge for Clayff atom types.
% # <charge_opls_go_atom.html charge_opls_go_atom(atom,Box_dim,varargin)> % Set the charge for specific OPLS atom types.
% # <check_clayff_2004_charge.html check_clayff_2004_charge(atom)> % Check the charge of the original Clayff atom types.
% # <check_clayff_charge.html check_clayff_charge(atom)> % Check the charge of the Clayff atom types.
% # <check_clayff_H2Odens.html check_clayff_H2Odens(atom,Box_dim)> % Check the approximate water density for a Clayff system.
% # <check_H2Odens.html check_H2Odens(atom,Box_dim)> % Compute the water density.
% # <clayff_2004_atom.html clayff_2004_atom(atom,Box_dim,varargin)> % Assign the original Clayff atom types by Cygan et al., 2004, with edge healing.
% # <clayff_2004_param.html clayff_2004_param(Atom_label,varargin)> % Hold ion and Clayff atom type parameters for the original Clayff force field.
% # <clayff_atom.html clayff_atom(atom,Box_dim,varargin)> % Assign Clayff atom types with edge healing.
% # <clayff_param.html clayff_param(Atom_label,varargin)> % Hold ion and Clayff atom type parameters.
% # <clayff210_atom.html clayff210_atom(atom,Box_dim,varargin)> % Assign modified Clayff atom types, with edge healing.
% # <clayff211_atom.html clayff211_atom(atom,Box_dim,varargin)> % Assign modified Clayff atom types (faster version), with edge healing.
% # <tweak_charge_atom.html tweak_charge_atom(atom)> % Tweak the charge of the atom struct to correct rounding errors.

%% INTERFACE from Heinz 2005, 2013, with atomtypes by MHolmboe
% # <charge_atom.html charge_atom(atom,Box_dim,ffname,watermodel,varargin)> % Charge the atom according to Clayff or Interface force fields.
% # <interface_atom.html interface_atom(atom,Box_dim,varargin)> % Assign atoms according to the Interface atom types, with modifications for edges.
% # <interface_param.html interface_param(Atom_label,water_model)> % Hold extended INTERFACE force field parameters.
% # <interface15_atom.html interface15_atom(atom,Box_dim,varargin)> % Assign atoms according to the Interface 1.5 atom types, with modifications for edges.
% # <interface15_param.html interface15_param(Atom_label,water_model)> % Hold extended INTERFACE 1.5 force field parameters.
% # <charge_interface_atom.html charge_interface_atom(atom,Box_dim,varargin)> % Set the charge for Interface atom types.
% # <charge_interface15_atom.html charge_interface15_atom(atom,Box_dim,varargin)> % Set the charge for Interface 1.5 atom types.
% # <check_interface_charge.html check_interface_charge(atom)> % Check the charge of the INTERFACE atom types.
% # <check_interface15_charge.html check_interface15_charge(atom)> % Check the charge of the INTERFACE 1.5 atom types.
% # <interface15_silica_atom.html interface15_silica_atom(atom,Box_dim,varargin)> % Assign atom types for the Interface 1.5 force field, specific to silica.
% # <tweak_charge_atom.html tweak_charge_atom(atom)> % Tweak the charge of the atom struct to correct rounding errors.

%% Graphene oxide modeled with OPLS/aa
% # <opls_go_atom.html opls_go_atom(atom,Box_dim,rmin,rlarge)> % Smear out the charge around -OH and epoxide groups in graphene oxide.
% # <oplsaa_go_param.html oplsaa_go_param(Atom_label,water_model)> % Hold the extended OPLS-AA force field parameters for graphite oxide.
% # <charge_opls_go_atom.html charge_opls_go_atom(atom,Box_dim,varargin)> % Set the charge for specific OPLS atom types.

%% Writing topology files for Lammps (Clayff) and Gromacs (Clayff/INTERFACE)
% # <write_atom_itp.html write_atom_itp(atom,Box_dim,filename_out,varargin)> % Create and print a Gromacs .itp file for Clayff or Interface force fields.
% # <write_atom_lmp.html write_atom_lmp(atom,Box_dim,filename_out,varargin)> % Create and print a LAMMPS data file (.lj) for Clayff systems.
% # <write_atom_oplsaa_go_itp.html write_atom_oplsaa_go_itp(atom,Box_dim,filename_out,varargin)> % Create and print a Gromacs .itp file for OPLS-AA or GO systems.

%% Bonded and Nonbonded Parameters
% # <bonded_parameters.html bonded_parameters(atom,varargin)> % Define bonded parameters for atoms.
% # <nonbonded_parameters.html nonbonded_parameters(atom,varargin)> % Define nonbonded parameters for atoms.
% # <nonbonded_ff.html nonbonded_ff(atom,varargin)> % Defines nonbonded forcefield parameters.

%% Lennard-Jones and Coulomb Potentials
% # <buckinghamcoul.html buckinghamcoul(atom,Box_dim,varargin)> % Calculates interactions using the Buckingham potential and Coulombic forces.
% # <ljcoul_12_6.html ljcoul_12_6(atom,Box_dim,varargin)> % Handles Lennard-Jones 12-6 potential along with Coulombic interactions.
% # <ljcoul_C12C6C4.html ljcoul_C12C6C4(atom,Box_dim,varargin)> % Lennard-Jones and Coulomb potential with C12, C6, C4 terms.
% # <ljcoul_C12C6.html ljcoul_C12C6(atom,Box_dim,varargin)> % Lennard-Jones and Coulomb potential with C12, C6 terms.
% # <ljcoul.html ljcoul(atom,Box_dim,varargin)> % General Lennard-Jones and Coulomb interaction function.
% # <ljcoul_2x2x.html ljcoul_2x2x(atom,Box_dim,varargin)> % Lennard-Jones and Coulomb potential with 2x factors.
% # <ljcoul_2x.html ljcoul_2x(atom,Box_dim,varargin)> % Variation of Lennard-Jones Coulomb with 2x factors.

%% Objective Functions and Force Calculations
% # <buckinghamcoul_objective_func.html buckinghamcoul_objective_func(atom,Box_dim,varargin)> % Objective function for fitting Buckingham and Coulomb potentials.
% # <ljcoul_objective_func.html ljcoul_objective_func(atom,Box_dim,varargin)> % Objective function for Lennard-Jones and Coulomb potentials.
% # <ljcoul_2x_objective_func.html ljcoul_2x_objective_func(atom,Box_dim,varargin)> % Objective function for 2x Lennard-Jones and Coulomb potentials.
% # <ljcoul_force.html ljcoul_force(atom,Box_dim,varargin)> % Computes forces based on Lennard-Jones and Coulomb potentials.
% # <ljcoul_2x_force.html ljcoul_2x_force(atom,Box_dim,varargin)> % Force calculation involving 2x Lennard-Jones and Coulomb potentials.
% # <ljcoul_force_C12C6C4.html ljcoul_force_C12C6C4(atom,Box_dim,varargin)> % Force calculation for Lennard-Jones potential with C12, C6, C4 terms.
% # <ljcoul_force_objective_func.html ljcoul_force_objective_func(atom,Box_dim,varargin)> % Objective function for force calculations with Lennard-Jones and Coulomb potentials.

%% Automated Fitting Tools
% # <autofit_C6C8C10xljcoul.html autofit_C6C8C10xljcoul(atom,Box_dim,varargin)> % Automated fitting for Lennard-Jones Coulombic parameters with C6, C8, C10 terms.
% # <autofit_2xljcoul_func.html autofit_2xljcoul_func(atom,Box_dim,varargin)> % Function for automatic fitting of 2x Lennard-Jones and Coulomb parameters.
% # <autofit_2xljcoul_batch.html autofit_2xljcoul_batch(atom,Box_dim,varargin)> % Batch fitting for 2x Lennard-Jones and Coulomb potentials.
% # <autofit_2xljcoul.html autofit_2xljcoul(atom,Box_dim,varargin)> % Automated fitting of 2x Lennard-Jones and Coulomb parameters.
% # <autofit_geometric2LB.html autofit_geometric2LB(atom,Box_dim,varargin)> % Automated geometric fitting for Lennard-Jones potential with Born-Mayer interactions.
% # <autofit_buckcoul.html autofit_buckcoul(atom,Box_dim,varargin)> % Automated fitting for Buckingham and Coulomb potentials.
% # <autofit_force_2xljcoul.html autofit_force_2xljcoul(atom,Box_dim,varargin)> % Automated fitting of force parameters for 2x Lennard-Jones and Coulomb interactions.
% # <autofit_ljcoul.html autofit_ljcoul(atom,Box_dim,varargin)> % Automated fitting for Lennard-Jones and Coulombic interactions.

%% Miscellaneous
% # <print_top.html print_top(atom,Box_dim,varargin)> % Prints or generates topology-related data.
% # <import_ff_table.html import_ff_table(filename,varargin)> % Imports forcefield parameter tables.
% # <change_top.html change_top(atom,Box_dim,varargin)> % Modifies the topology file or its parameters.
% # <mass_atom_clayff.html mass_atom_clayff(atom)> % Defines or calculates masses for Clayff atoms.
% # <harmonic_bond.html harmonic_bond(atom,Box_dim,varargin)> % Defines or computes harmonic bond forces or parameters.
% # <smear_charge.html smear_charge(atom,Box_dim,varargin)> % Distributes charge across atoms, possibly using charge smearing techniques.
% # <check_interface_charge.html check_interface_charge(atom)> % Checks charges for INTERFACE atom types.
% # <check_interface15_charge.html check_interface15_charge(atom)> % Checks charges for INTERFACE 1.5 atom types.

%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
