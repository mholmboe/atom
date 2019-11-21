%% Examples demonstrating how to center, translate, rotate, place a molecule
% (For a full list of functions that can move around molecules, go to 
% <List_build_functions.html List_build_functions> or the 
% <List_general_functions.html List_general_functions>)

%% First set some convenient matlab settings
format compact; set(gcf,'Visible','on');

%% Pick filenames to import and export
% Set some filenames
filename_in='Ethanol.pdb'; % default is 'Pyrophyllite.pdb'
filename_out='outEthanol'; % default is 'test'

%% First import some molecule
atom=import_atom(filename_in);

%%
% *Note* that you can issue plot_atom(atom,Box_dim) or vmd(atom,Box_dim)
% (if VMD is installed and Matlab knows the path to it) to see the molecule

%% Center a molecule
% Use <center_atom.html center_atom> to center the molecule in the middle
% of the box.
% *Note* that you could pass two additional arguments that can be used to 
% center only a specific resname (molecule name) and/or along specific 
% x|y|z dimensions. Look for instance into <center_atom.html center_atom> 
% and examples 2 and 3.

centered_atom=center_atom(atom,Box_dim);
% plot_atom(centered_atom,Box_dim) % Run command to see what changed

%% Translate a molecule
% Use <translate_atom.html translate_atom> to translate the molecule 
% somewhere.
% *Note* that you could pass one additional argument that can be used to 
% translate only a specific resname (molecule name). Look for instance into 
% <translate_atom.html translate_atom>  and examples 2 and 3.

translated_atom=translate_atom(atom,[0 5 10]);
% plot_atom(translated_atom,Box_dim) % Run command to see what changed

%% Place a molecule
% Use <place_atom.html place_atom> to place a molecule 
% somewhere
% *Note* that the <place_atom.html place_atom> function is dependent on
% the <COM_atom.html COM_atom> function and is therefore a bit slow for
% large molecules

placed_atom=place_atom(atom,[0 5 10]);
% plot_atom(placed_atom,Box_dim) % Run command to see what changed

%% Rotate a molecule
% Use <rotate_atom.html rotate_atom> to rotate the molecule. The second 
% argument can either be the string 'random' or a 1x3 array holding the 
% new rotation angles around the x,y,z axes.
% *Note* that the <rotate_atom.html rotate_atom> function is dependent on
% the <COM_atom.html COM_atom> function and is therefore a bit slow for
% large molecules

rotated_atom=rotate_atom(atom,Box_dim,'random');
% or to rotate with specified angles around the x,y,z axes.
rotated_atom=rotate_atom(atom,Box_dim,[0 90 180]); 
% plot_atom(rotated_atom,Box_dim) % Run command to see what changed

