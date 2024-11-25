%% Translate a molecule
% Let's import some molecuel and translate it with a 1x3 translation vector.

close all;
atom=import_atom('6x4MMT.gro');
atom=translate_atom(atom,[10 5 -5]); % [10 5 -5] is the translation vector
plot_atom(atom,Box_dim);
view(85,15);
%vmd(atom,Box_dim)