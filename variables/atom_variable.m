%% atom
% * The atom struct, the main variable in the atom MATLAB library.
%
%% Version
% 2.11
%

%% Introduction to the atom struct
% The <atom_variable.html atom> struct variable (as in Atomistic Topology Operations 
% in Matlab) is an indexed Matlab struct variable that stores molecule 
% attributes like atomtype names, moleculeID's, atomID's, coordinates
% and more as indexed fields in the atom struct. A brief list of the atom 
% struct fields are listed below.
% Normally all strings are cell strings and all numeric variables are
% double-precision floating-points. Unassigned values are set to 'NaN' or 0
% . *_atom_* is the default variable name for the atom struct, however 
% any string name (not starting with a numeric value) would also work, like 
% *_DNA_*, *_Ions_* or *_Water_* etc.
%
% Furthermore, different atom structs for instance representing a DNA 
% molecule, a set of counter-ions and solvating waters, can easy be 
% concatenated/appended to each other by issuing:
% 
% FullSystem = [DNA Ions Water]; % Caveat: all fields must be the same
%
% However, if one wants to update all the moleculeID's and atomID's, it is
% recommended to use the <update_atom.html update_atom> function like this:
%
% FullSystem = update_atom({DNA Ions Water}); % Updating the molid and
% index fields
%
% If invoking other functions like the
% <bond_angle_atom.html bond_angle_atom> or similar functions it can also 
% hold other optional fields carrying the bonding and angle info to the 
% nearest neighbours etc, or the mass, charge, bond valence and so on. 
% Since the atom  struct is indexed, atom(n) holds the attributes of the 
% n:th particle in the atom struct. Because of indexing, the particles in 
% the atom struct can easily be filtered and manipulated - based on a 
% selection of specific or a range of attributes, like atomtype names, 
% coordinates and so on. A few examples are shown below... see the 
% <List_all_functions.html List_all_functions> for all functions that
% operate on the <atom_variable.html atom> struct. Note that many of these functions
% also need the often accompaning <Box_dim_variable.html Box_dim> variable.
%
% *Default fields*
%
% # *molid* (molecule/residue ID number)
% # *resname* (molecule/residue name)
% # *type* (atomname/type)
% # *fftype* (secondary atomname/type)
% # *index* number (atomID number)
% # *neigh* (holds different subfields)
% # *bond* (holds different subfields)
% # *angle* (holds different subfields)
% # *x* (X-coordinate in Ångström)
% # *y* (Y-coordinate in Ångström)
% # *z* (Z-coordinate in Ångström)
% # *vx* (X-velocity in Ångström/ps)
% # *vy* (Y-velocity in Ångström/ps)
% # *vz* (Z-velocity in Ångström/ps)
%
% *Non default fields* that can be set after invoking other
% functions
%
% * *COM_x* (the correspondingn MolID's center-of-mass in x in Ångström)
% * *COM_y* (the correspondingn MolID's center-of-mass in y in Ångström)
% * *COM_z* (the correspondingn MolID's center-of-mass in z in Ångström)
% * *element* (element name)
% * *mass* (atomic mass in g/mol)
% * *charge* (atomic partial charge)
% * and so on...
%
%% atom struct examples
% Import some molecule from an .pdb or .gro file with
% atom=import_atom(filename), then try this some of these or similar
% commands relevant for your molecule
%
%% To find certain atomtypes... can you spot the differences between the different commands?
% # index=ismember([atom.type],[{'Al' 'Alt' 'Mgo'}]) % gives a binary (1/0) logical array
% # index=find(ismember([atom.type],[{'Al' 'Alt' 'Mgo'}])) % finds the values of matching indexes
% # index=strcmp([atom.type],'Al') % try also strncmp or strncmpi?
% # index=find(strncmpi([atom.type],'al',2) % Will find the indexes of 'Al' 'Alt'
% 
% Now create a new atom struct with the selected atomtypes 
%
% # new_atom=atom(index) % This creates a new_atom struct with the filtered/selected atomtypes

%% To select atoms/particles/sites based on their coordinates
% # positive_z_atom = atom([atom.z]>0) % finds all atoms with a positve z-coordinate
% # inbetween10n20_atom = atom([atom.z]>10&[atom.z]<20) % finds all atoms with a positve z-coordinate
% # first100_atom = atom(1:100) % finds the first 100 atoms in the atom struct
% # first100_atom = atom([atom.index]<101) % also finds the first 100 atoms in the atom struct

%% Replace some atomtype's or resname's
% # [atom(index).type]=deal({'Mg'})
% # [atom(index).resname]=deal({'MIN'})

%% Try this: 
atom=import_atom('Ethanol.pdb') % Will also output the Box_dim variable and some other stuff, see the variable explorer
atom=center_atom(atom,Box_dim) % Center the coordinates, i.e. [atom.x|.y|.z]
show_atom(atom,Box_dim,.1,1) % Plot the atom with the Box_dim

