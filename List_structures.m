%% List of available solvents
% .pdb or .gro files to be used with the function solvate_atom. Note that you could also use
% custom solvent boxes to solvate a simulation cell.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%

%% water
% # 864_spc.gro|.pdb % equilibrated SPC water box
% # 864_spce.gro|.pdb % equilibrated SPC/E water box
% # 864_tip3p.gro|.pdb % equilibrated TIP3P water box
% # 864_tip4p.gro|.pdb % equilibrated TIP4P water box
% # 864_tip5p.gro|.pdb % equilibrated TIP5P water box
% # 96spc_hex_ice_h.gro|.pdb % equilibrated SPC hex-ice water box
% # 96tip4p_hex_ice_h.gro|.pdb % equilibrated TIP4P hex-ice water box
% # 864_swm4_ndp.gro|.pdb % Polarizable water v1
% # 864_swm4_ndp_vds.gro|.pdb % Polarizable water v2

%% organics
% # 500xEtOH.gro|.pdb % equilibrated ethanol solvent box

%% mineral
% # 1xPyro_Lee_Guggenheim_1981_alfabeta90.pdb % Pyrophyllite mineral structure from Lee and Guggenheim, 1981
% # Pyrophyllite.pdb % Pyrophyllite mineral structure
% # 3WNaMMT.pdb % Montmorillonite (MMT) mineral structure
% # 3WNaMMT_4x6x2.pdb % Large Montmorillonite (MMT) mineral structure
% # - Hexagonal_pyro_mmt % Hexagonal pyrophyllite and montmorillonite mineral structures
% # - Hexagonal_talc_laponite % Hexagonal talc and laponite mineral structures

%% Conversion functions
% # <spc2tip4p.html spc2tip4p(filename)> % This function converts a .gro or .pdb file with spc water to some tip4p water
% # <spc2tip5p.html spc2tip5p(filename)> % This function converts a .gro or .pdb file with spc water to some tip5p water
% # <spce2tip4p.html spce2tip4p(filename)> % This function converts a .gro or .pdb file with spce water to some tip4p water
% # <tip3p2tip4p.html tip3p2tip4p(filename)> % This function converts a .gro file with tip3p water to tip4p water
