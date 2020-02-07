%% List of available solvents
% .pdb or .gro files to be used with the function solvate_atom. Note that you could also use
% custom solvent boxes to solvate a simulation cell.
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Water structures
% 864_spc.gro|.pdb % equilibrated SPC water box
% 864_spce.gro|.pdb % equilibrated SPC/E water box
% 864_tip3p.gro|.pdb % equilibrated TIP3P water box
% 864_tip4p.gro|.pdb % equilibrated TIP4P water box
% 864_tip5p.gro|.pdb % equilibrated TIP5P water box
% 96spc_hex_ice_h.gro|.pdb % equilibrated SPC hex-ice water box
% 96tip4p_hex_ice_h.gro|.pdb % equilibrated TIP4P hex-ice water box
% 864_swm4_ndp.gro|.pdb * Polarizable water v1
% 864_swm4_ndp_vds.gro|.pdb * Polarizable water v2

%% Other structures
% 864_spc.gro|.pdb % equilibrated SPC water box

%% Conversion functions
% # <spc2tip4p.html spc2tip4p(filename)> % This function converts a .gro or .pdb file with spc water to some tip4p water
% # <spc2tip5p.html spc2tip5p(filename)> % This function converts a .gro or .pdb file with spc water to some tip5p water
% # <spce2tip4p.html spce2tip4p(filename)> % This function converts a .gro or .pdb file with spc water to some tip4p water
% # <tip3p2tip4p.html tip3p2tip4p(filename)> % This function converts a .gro file with spc water to some tip4p water