%% Examples demonstrating how to plot the atom struct
% (See also the list
% <List_general_functions.html List_general_functions>)

%% First set some convenient matlab settings
format compact; set(gcf,'Visible','on');

%% Pick filenames to import and plot
% Set some filenames
filename_in='Pyrophyllite.pdb'; % default example is 'Pyrophyllite.pdb'

%% Import some molecule
atom=import_atom(filename_in);
% atom=replicate_atom(atom,Box_dim,[4 2 1]) % Replicate the molecule just to get a bigger system

%% Plot using <vmd.html vmd>
% The <vmd.html vmd> function is basically a redirect to the popular VMD
% molecular visualization software, that must be installed separately. Hence
% look into the <vmd.html vmd> function and set the path to VMD so Matlab 
% can find it. Typically on a Mac it looks something like this:
%
% PATH2VMD = '/Applications/VMD\ 1.9.2.app/Contents/MacOS/startup.command';
%
% The <vmd.html vmd> function can take the structure filename as input, or
% the <atom_variable.html atom> variable, and/or the <Box_dim_variable.html
% Box_dim> variable as input.
%%
% *Examples*
%
vmd(filename_in)
vmd(atom)
vmd(atom,Box_dim)

%% Plot using <plot_atom.html plot_atom>
% The <plot_atom.html plot_atom> function is a really simplistic viewer
% function that can show the <atom_variable.html atom> struct in a coordinate system.
% If first invoking the <bond_angle_atom.html bond_angle_atom> function
% generating the <Bond_index_variable.html Bond_index> variable, all bonds 
% (also across the PBC) can be plotted if supplying the
% <Bond_index_variable.html Bond_index> variable.

%%
% *Examples*
%
plot_atom(atom,Box_dim)
plot_atom(atom,Box_dim,2) % 2 is a scale factor for the plotted atoms
plot_atom(atom,Box_dim,1,Bond_index) % Bond_index from bond_angle_atom()
plot_atom(atom,Box_dim,1,[],'axis') % Will show an axis in the lower left corner

%% Plot using <plot_density_atom.html plot_density_atom>
% The <plot_density_atom.html plot_density_atom> function is a really 
% simplistic viewer function that can show the <atom_variable.html atom> 
% struct in a coordinate system, along with density profiles that has been 
% smoothed through a Gaussian convolution. Note that the system is wrapped
% before calculating the density profiles.

%%
% *Examples*
%
plot_density_atom(atom,Box_dim)
plot_density_atom(atom,Box_dim,2) % 2 is a scale factor for the plotted atoms
plot_density_atom(atom,Box_dim,1,3) % Here the second argument 3 scales the plotted density profiles


