%% Replicate a molecule
% First import a molecule...

atom=import_atom('1xMMT.gro');
plot_atom(atom,Box_dim); hold off;
%vmd(atom,Box_dim)
%% 
% Then replicate it like this...

atom_replicated=replicate_atom(atom,Box_dim,[6 4 1]); % This will also replicate Box_dim
plot_atom(atom_replicated,Box_dim); hold off;
%vmd(atom,Box_dim)
%% 
% There is also an option to replicate the molecule along the x/y/z directions 
% in a different order, like this

atom=import_atom('1xMMT.gro');
atom_replicated_yxz=replicate_atom(atom,Box_dim,[6 4 1],'yxz'); % This will also replicate Box_dim
plot_atom(atom_replicated_yxz,Box_dim);
%vmd(atom,Box_dim)