%% Wrap a molecule
% First let's import and translate some molecule...

close all;
atom=import_atom('6x4MMT.gro');
atom=translate_atom(atom,[10 5 -5]); % [10 5 -5] is the translation vector
show_atom(atom,Box_dim); view(85,15)
%vmd(atom,Box_dim)
%% 
% Let's wrap every atom back into the box!

close all;
atom=wrap_atom(atom,Box_dim);
show_atom(atom,Box_dim); view(85,15);
%vmd(atom,Box_dim)