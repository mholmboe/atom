%% Examples demonstrating how to import (i.e. read) and export (i.e. write) structure files, like .pdb|.gro|.xyz files
% (For a full list of importing and exporting functions, go to 
% <List_import_functions.html List_import_functions> or the 
% <List_export_functions.html List_export_functions>)

%% First set some convenient matlab settings
format compact; set(gcf,'Visible','on');

%% Pick filenames to import and export
% Set some filenames
filename_in='Pyrophyllite.pdb'; % default is 'Pyrophyllite.pdb'
filename_out='outPyrophyllite'; % default is 'outPyrophyllite'

%% How to read/import a structure file into Matlab
% Use the <import_atom.html import_atom> function as demonstrated below.
% This function automatically detects the file format, which should be
% either .pdb|.gro|.xyz

%%
% *Example* with <import_atom.html import_atom>
atom = import_atom(filename_in);
% plot_atom(atom,Box_dim) % Run command to see what changed

%%
% *Note* that you get information about the 
% composition and Box dimensions if such exist.

%   Found .pdb file
%   filename =
%       'Pyrophyllite.pdb'
%   .pdb file imported
%   composition = 
%       struct with fields:
% 
%      resnames: {'PYR'}
%     nresidues: 1
%        natoms: 40
% Atom_types =
%   1×4 cell array
%     {'Al'}    {'H'}    {'O'}    {'Si'}
% Atom_numbers =
%      4     4    24     8
% Box_dim =
%     5.1600    8.9658    9.1897         0         0         0         0   -1.6966         0

%%
% *Example* with <import_atom_pdb.html import_atom_pdb> (or <import_atom_gro.html import_atom_gro> or <import_atom_xyz.html import_atom_xyz>)
atom = import_atom_pdb(filename_in);

%%
% *Note* that with 
% <import_atom_pdb.html import_atom_pdb> or the
% <import_atom_gro.html import_atom_gro> or the 
% <import_atom_xyz.html import_atom_xyz>) functions you get less 
% information than with the more general <import_atom.html import_atom> 
% function.
%
% *Note* that when you import a structure file, you can pass one or two
% additional arguments after the filename to (1) translate, or (2)
% center and translate the structure to a new box. Look for instance at
% <import_atom_pdb.html import_atom_pdb> and examples 2 and 3.

%% How to write/export a structure file
% In order to export or write out a structure file, we need to specify the 
% format we want (like .pdb|.gro|.xyz) by calling the corresponding
% function  <write_atom_pdb.html write_atom_pdb> or the
% <write_atom_gro.html write_atom_gro> or the 
% <write_atom_xyz.html write_atom_xyz>. Apart from passing the
% <atom_variable.html atom> struct, we must also pass the
% <Box_dim_variable.html Box_dim> variable and an output filename.
%
% *Note* that when you export a .pdb file, the first cutoff is
% <rmaxshort_variable.html rmaxshort> for bonded H's and the 
% second cutoff is (<rmaxlong_variable.html rmaxlong>) for all other 
% bonds

write_atom_gro(atom,Box_dim,filename_out)
write_atom_pdb(atom,Box_dim,filename_out)
write_atom_pdb(atom,Box_dim,filename_out,1.25,2.25) % Will try to write the CONECT bond records
write_atom_xyz(atom,Box_dim,filename_out) % With the Box_dim, not standard  for .xyz files
write_atom_xyz(atom,filename_out) % Without the Box_dim, standard for .xyz  files


    
    