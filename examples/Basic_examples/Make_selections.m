%% Examples demonstrating how to make selections of filter the atom struct
% (For a full list of functions that select/filter the 
% <atom_variable.html atom> struct, go to 
% <List_build_functions.html List_build_functions> or the 
% <List_general_functions.html List_general_functions>)

%% First set some convenient matlab settings
format compact; set(gcf,'Visible','on');

%% Pick filenames to import and export
% Set some filenames
filename_in='Pyrophyllite.pdb'; % default example is 'Pyrophyllite.pdb'

%% Import some molecule
atom=import_atom(filename_in);
atom=replicate_atom(atom,Box_dim,[4 4 4]) % Replicate the molecule just to get a bigger system

%% Different types of string based selections
% The following selection commands can be applied both to the fields .type
% and .resname. In these examples we first create some index variable ind, 
% and then use this ind variable to extract all the matching atoms.

% Select atoms in the atom struct matching certain atomtypes using Matlabs 
% ismember() function. 
ind1 = ismember([atom.type],[{'Al' 'Si' 'H'}]) % gives a binary (1/0) logical array
atom1 = atom(ind1) % This creates a new atom1 struct with the filtered/selected atomtypes

% Similar for one atomtype using Matlabs strcmp() function
ind2 = strcmp([atom.type],'Al') % try also strncmp or strncmpi?
atom2 = atom(ind2) 

% Similar for one atomtype using Matlabs strncmpi() function, which is case
% insensitive and in this example only tries to match the 2 first
% characters
ind3 = find(strncmpi([atom.type],'al',2)) % Will find the indexes of 'Al' 'Alt'
atom3 = atom(ind3)

% Compare ind1, ind2 and ind3. Are they similar or different? What is the
% difference between ismember(), strcmp(), strncmpi()? Could we have used
% any other Matlab functions to do the same thing?

%% Different types of selections based on numeric values
% In these examples we skip the generation of index variables used above, 
% just because we can...

% Select atoms based on some coordinate values
pos_z_atom = atom([atom.z]>0) % finds all atoms with a positve z-coordinate
% or try this
specific_z_atom = atom([atom.z]>1|[atom.z]<-1) % Note the | character

% Select atoms based on their index
first100_atom = atom(1:100) % finds the first 100 atoms in the atom struct
first100_atom = atom([atom.index]<101) % also finds the first 100 atoms in the atom struct
molecule3_atom = atom([atom.molid]==3) % finds the atoms in the 3rd molecule (having .molid = 3)
