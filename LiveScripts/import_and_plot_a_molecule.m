%% Import and plotting a structure
%% Import a structure
% This is how we import any .gro | .pdb | .xyz file

atom=import_atom('96spc_hex_ice_h.gro');
%% 
% ..or you could try these commands specific for each filetype

% atom=import_atom_gro('filename.gro');
% atom=import_atom_pdb('filename.pdb');
% atom=import_atom_xyz('filename.xyz'); % any info regarding box size? If not add it manually to the Box_dim 1x3 (1x9 for triclinic) vector
%% Now let's plot the imported simulation box
% The atom scripts comes with a simple plot function called plot_atom, which 
% can display and rotate a simulation cell in 3D.

show_atom(atom,Box_dim,'halfvdw','box'); % 'axis' argument is optional
view(85,5); % change the default view
%% 
% However, you could also use VMD if the path in the function PATH2VMD() is 
% set correctly. Note that VMD will lock the matlab editor/parser as long a its 
% open, but can be closed by issuing 'exit' in the command window. Anyway,to run 
% VMD try invoking:

vmd(atom,Box_dim);
%% 
% or plot the file direclty by issuing something like:

% vmd('96spc_hex_ice_h.gro');