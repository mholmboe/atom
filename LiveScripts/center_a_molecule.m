%% Center a molecule
% Import a structure and center it

close all;
atom=import_atom('6x4MMT_x10y5z0.gro'); % A non-centered structure
atom=center_atom(atom,Box_dim); % This will also replicate Box_dim
show_atom(atom,Box_dim);
%vmd(atom,Box_dim)
%% 
% We can also center with respect to a certain resname group (here |'MMT'|)

close all;
atom=import_atom('6x4MMT_x10y5z0.gro');
atom=center_atom(atom,Box_dim,'MMT'); % Default 'all'
show_atom(atom,Box_dim);
%vmd(atom,Box_dim)
%% 
% Or we can center with respect a certain x, y or z direction (here in the |'z'| 
% direction)

close all;
atom=import_atom('6x4MMT_x10y5z0.gro');
atom=center_atom(atom,Box_dim,'MMT','z'); % Default 'xyz'
show_atom(atom,Box_dim);
%vmd(atom,Box_dim)