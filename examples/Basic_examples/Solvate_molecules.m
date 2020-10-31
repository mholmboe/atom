%% Examples demonstrating how to solvate with water or any custom solvent box
% (For a full list of water models and solvents 
% <List_structures.html List_structures>.

%% First set some convenient matlab settings
format compact; set(gcf,'Visible','on');

%% Pick a filename/structure to import 
filename_in='Pyrophyllite.pdb'; % default is 'Pyrophyllite.pdb'

%% Solvating molecules using <solvate_atom.html solvate_atom>
% The function <solvate_atom.html solvate_atom> can solvate either a
% box or an arbitrary orthogonal volume defined by the 
% <limits_variable.html limits> variable, by stacking pre-equilibrated 
% boxes in all directions until the box or specified volume is completeley 
% filled with _maxsol_ number of solvent molecules. If an existing 
% <atom_variable.html atom> struct is passed along (called solute_atom in 
% the examples below), any solvent molecules having atoms with a cutoff of 
% _r_ Ångström will be removed in order to remove atomic overlap. [] can be
% used to as to pass along an empty solute_atom struct. Note that
% the relative density of the solvent box to be stacked can be set with the
% variable _density_ (default 1, but 1.1 is usually also ok). Optionally a 
% string can be passed the specify the desired water model, like 'SPC' or 
% 'TIP4P' (see also <List_structures.html List_structures>), or 'custom' -
 % where the following two arguments must be a custom preequilibrated 
% solvent atom struct and the corresponding <Box_dim_variable.html Box_dim> 
% variable, see the last example.

%% First import some molecule and call it solute_atom
solute_atom = import_atom(filename_in); 
solute_atom = replicate_atom(solute_atom,Box_dim,[4 2 1]) % Replicate the 
% molecule just to get a bigger molecule

%% Set some variables
%
limits = [20 20 20] % The 1x3 or 1x6 <limits_variable.html limits> variable representing the volume to be solvated. Can be set to the <Box_dim_variable Box_dim>.
density = 1.1 % Relative density of solvent, one can use >1 to squeeze in extra...
r = 2 % Ångström, nearest solute - solvent distance 
maxsol = 100;

%% Run the function <solvate_atom.html solvate_atom>
SOL = solvate_atom(limits,density,r,maxsol) % Will solvate an empty box
SOL = solvate_atom(limits,density,r,maxsol,solute_atom) % Will solvate the solute_atom struct
SOL = solvate_atom(limits,density,r,'maxsol',solute_atom) % Will maximize the number of solvent molecules
SOL = solvate_atom(limits,density,r,'shell15',solute_atom) % Will add a 15 Ångström thick shell around the solute_atom
SOL = solvate_atom(limits,density,r,maxsol,solute_atom,'tip4p') % Will solvate the solute_atom struct with tip4p water
SOL = solvate_atom(limits,density,r,maxsol,solute_atom,'spc_ice') % Will solvate the solute_atom struct with a SPC ice structure

%%
% To solvate with a custom solvent, import some preequlibrated box and
% name its atom struct to _mysolvent_ and its _Box_dim mysolvent_Box_dim_.
% Note that this will overwrite any existing <Box_dim_variable.html Box_dim>
% variable, hence in practise you may want to import your solute molecule 
% after the custom solvent.
mysolvent=import_atom('500xEtOH.pdb');
mysolvent_Box_dim=Box_dim;
maxsol=30; % Since the ethanol molecule is larger than a water molecule, we need to decrease _maxsol_ to 30
SOL = solvate_atom(limits,density,r,maxsol,solute_atom,'custom',mysolvent,mysolvent_Box_dim) % Will solvate the solute_Atom with a custom solvent, like an ethanol box

%% Add the solvent to the initial molecule
System = update_atom({solute_atom SOL});
% plot_atom(System,Box_dim) % Now plot the final solvated system
% vmd(System,Box_dim) % If VMD is installed

