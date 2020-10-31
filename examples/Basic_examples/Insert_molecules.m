%% Examples demonstrating how to add new atoms or ions to an atom struct
% (For a full list of functions that can add atoms/ions/molecules, go to 
% <List_build_functions.html List_build_functions>

%% First set some convenient matlab settings
format compact; set(gcf,'Visible','on');

%% Create new atoms with <create_atom.html create_atom>
% This function creates new atoms/particles within a certain region defined
% by <limits_variable.html limits>. It can also add particles on a plane by setting Lx|Ly|Lz to 
% 0 or something small. We call the function by issuing something like:
%
% atom = create_atom(type,resname,limits,nmax,varargin) 
%
% *Input arguments*
%
% * {type} is particle/atomtype, like {'Si'}
% * {resname} is resname, like {'MIN'}
% * [limits] is a 1x3 or 1x6 array representing a volume, see also <limits_variable.html limits>
% * The number nmax is the max number of particles
% * Optional scale arguemnt (varargin{1}) is a how-many-diameters-between-the-particles-thingy

%%
% *Examples*
atom = create_atom('Na','Na',[10 20 30],10) % Adding 10 Na in a 10x20x30 Å box
atom = create_atom('Na','Na',[0 0 0 10 20 30],10) % same as above
atom = create_atom('Na','Na',[0 0 30 10 20 30],10) % Adding 10 Na on a plane at z=30Å in a 10x20x30 Å box
atom = create_atom('Na','Na',[10 20 30],10,2) % Double the minimum spacing between the Na atoms
atom = create_atom('Na','Na',[10 20 30],10,1.5,in_atom) % Making sure there is no atomic overlap with another in_atom struct

%% Insert new atoms or molecules with <insert_atom.html insert_atom>
% * This function inserts a whole molecule from a structure file or atom_in
% into a region defined by <limits> with an atom or molecule. varargin 
% can be used to assure that one atom type is at least some distance above 
% (in z) some other atom type
% 
% atom = insert_atom(atom_in,limits,rotate,r,nmax,varargin)
%
% *Input arguments*
%
% * atom_in is the atom struct to be inserted
% * [limits] is a 1x3 or 1x6 array representing a volume, see also 
% <limits_variable.html limits>
% * rotate argument can be a string like 'random', {'random'}, or be used
% to set some angles like [60 90 60].
% * The number maxsol is the max number of atoms/molecules to insert
% * Optional varargin{1} is any already existing atom struct with which we 
% want to avoid any atomic overlap with

%%
% *Examples*
atom_in = import_atom('Ethanol.pdb'); % Import the example Ethanol molecule
limits=[30 30 30]; % Set some volime defined by <limits_variable limits>
r=2; % rmin cutoff
maxsol=10; % max number of atoms/molecules to insert
atom1 = insert_atom(atom_in,limits,'rotate',r,maxsol)
atom2 = insert_atom(atom_in,limits,[10 20 30],r,maxsol,[])
atom3 = insert_atom(atom_in,limits,'rotate',r,2,atom_in,[1 2],0.3) % Here we make sure the first atomtype in a molecule has a z-coordinates 0.3 Å > than the second atomtype within the same molecule

%% Add ions with <ionize_atom.html ionize_atom>
% * This function adds ions within a certain region defined by <limits>
% * Can also add particles on a plane by setting Lx|Ly|Lz to 0 or something small
% * Compared to create_atom, this function can also add particles near a
% 'surface' or in the 'bulk', when an in_atom struct (representing a solid
% phase) is passed argument.
% * If slow, check out insert_atom or solvate_atom or grid_atom...
%
% *Input arguments*
%
% * {type} is particle/atomtype, like {'Na'}
% * {resname} is resname, like {'Na'}
% * [limits] is a 1x3 or 1x6 array representing a volume, see also <limits_variable.html limits>
% * The number nmax is the max number of particles
% * Optional scale argument (varargin{1}) is a how-many-diameters-between-the-particles-thingy
% * in_atom is any already existing atom struct (like a surface) with which 
% we want to avoid any atomic overlap with

% 
%%
% *Examples*
atom = ionize_atom('Na','Na',[10 20 30],10)
atom = ionize_atom('Na','Na',[10 20 30],10,2) % Nearest distance will be 2 * ionic radii
atom = ionize_atom('Na','Na',[10 20 30],10,2,in_atom) % Random placement without atomic overlap with in_atom
atom = ionize_atom('Na','Na',[10 20 30],10,2,in_atom,'surface') % Preferred placement at the 'surface' or 'bulk'
atom = ionize_atom('Na','Na',[10 20 30],10,2,in_atom,'surface'|'bulk','x'|'y'|'z'|20) % Preferred placement at the 'surface' or'bulk' within the x|y|z or [value] range
