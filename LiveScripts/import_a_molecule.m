%% How-to import a .gro/.pdb/.xyz structure file and place it
% Here we will show how you could import a molecule and translate it somewhere 
% in the original simulation cell or in a new simulation cell having new dimensions.
% 
% 
% 
% Imports the structure file into matlabs variable space

% atom=import_atom(filename) 
atom=import_atom('1xMMT.gro');
%% 
% Import a structure and move/translate it by some x,y,z vector

% atom=import_atom(filename,[translation vector 1x3]) – Imports a structure and translates it
atom=import_atom('1xMMT.gro',[0 0 10]);

%% 
% 
% 
% Import a structure and center it in a new box and move it by some x,y,z vector, 
% i.e. imports a structure and centers it in the x/y plane at z=0, then translates 
% it with the translation vector 

% atom=import_atom(filename,[translation vector 1x3],[New_Box_dim 1x3]) 
atom=import_atom('1xMMT.gro',[0 0 10],[10 20 20]);