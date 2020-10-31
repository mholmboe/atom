%% Import a .car file from Heinz's MODEL database in his INTERFACE 1.5 distribution
atom = import_atom_car('mont0_333_K_15_single_layer.car','no_counterions') % Remove the counterions, we will add them later when building the enitre system

% Is the charges correct?
% Before starting the simulation, make sure the topol.top file is correct,
% and note that you must change the moleculename in each .itp file, i.e.
% mon should be changed to MMT

