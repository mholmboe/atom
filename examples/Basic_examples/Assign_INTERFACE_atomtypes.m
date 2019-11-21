%% Example on how to assign Interface FF atomtypes to a mineral atom struct
% See also the <Advanced_examples.html advanced examples>

%%
% Hendrik Heinz and co-workers distribute their forcefield along with a 
% large set of example structures and topology files that can be converted 
% for Lammps, Namd2 and Charmm. Here we will try to assign Interface FF 
% atomtypes to a mineral atom struct and also generate a Gromacs molecule 
% topology file, an so-called .itp file, along with a normal .pdb file 
% Hence start by downloading the 
% <https://bionanostructures.com/interface-md/ Interface FF package>.
%
% Note that the Interface FF is a fairly complex forcefield since 
% it contains many different atomtypes having well-defined bonds and angles.
% Because of this the functions presented below and the accuracy of the 
% output files they generate is highly uncertain, and you need to verify 
% them yourself by for instance reading up on the relevant Interface FF 
% publications. I take no responsibility...

%% First set some convenient matlab settings
format compact; set(gcf,'Visible','on');

%% Method 1: Convert a .car file from the official Interface FF distribution (and optionally write a some topology files)
% From the downloaded forcefield package, grab a .car file from the
% MODEL_DATABASE/CLAY_MINERALS directory, like the
% 'mont0_333_Na_15_cell.car' file which we can convert with the 
% <import_atom_car.html import_atom_car> function into a .pdb file as well 
% as two types of topology files, an .itp-file (for Gromacs) and a .psf 
% file (for Namd2). Note that this function only works
% with atomtypes also listed in the actual function and its dependencies. 
% If you find something is missing or if you cannot get it to work, 
% email michael.holmboe@umu.se.

atom = import_atom_car('mont0_333_Na_15_cell.car')

%%
% Some of the structures contain counterions, to remove them run:
atom = import_atom_car('mont0_333_Na_15_cell.car','no_counterions')

%%
% What output was generated? What can you do with it?

%% Method 2: Use the <interface_atom.html interface_atom> function to assign the Interface FF atomtypes
%

%% Set input and output filenames
filename_in='Pyrophyllite.pdb'; % default is 'Pyrophyllite.pdb'
filename_out='6x4x1_MMT_interface_2004.pdb'; % default is 'Pyrophyllite.pdb'

%% Import a unit cell structure and replicate it into a mineral layer
% Normally, when constructing a mineral particle and assigning the 
% Interface FF (Heinz et al, 2005,2013) atomtypes, one usually start of by 
% building the mineral particle from an X-ray determined unit cell 
% structure. In this example we will demonstrate how to generate a 
% montmorillonite layer particle from a pyrophyllite unit cell. Both 
% minerals are so-called 2:1 T-O-T sheet silicate minerals, with the 
% difference between the two isostructural minerals being (as you may know) 
% the fact that montmorillonite carries charge defects due to isomorphic 
% substitution. For montmorillonite this means substitution of Si4+ with 
% Al3+ in the two tetrahedral sheets, or octahedral Al3+ with Mg2+ or Fe2+ 
% in the octahedral sheet of the montmorillonite layer.

atom = import_atom(filename_in); % Imports a pyrophyllite unit cell
atom = replicate_atom(atom,Box_dim,[6 4 1]); % Replicate the structure by 6x4x1 into a Pyrophyllite clay layer

%% Perform isomorphic substitution using <substitute_atom.html substitute_atom>
% <substitute_atom.html substitute_atom> was written with centrosymmetric
% minerals in mind, therfore it works best if the layer is centered at z=0 
% in the xy-plane. It can handle both substitutions in the octahedral as
% well as in the tetrahedral sheets. It distributes the substituted sites
% randomly, except for the fact that one can choose a nearest distance
% between the substituted sites, which in this example is 5.5 Ångström. We
% do this to avoid oxygen atoms facing two substituted sites at the same
% time (Google Löwenstein's rule).

atom = substitute_atom(atom,Box_dim,6*4*2/3,'Al','Mgo',5.5) % Perform octahedral (only) substitutions on 2/3's of all Al sites. Here 5.5 is the minimum distance between the substituted Mg2+ sites
% atom = substitute_atom(atom,Box_dim,14,'Al','Mgo',5.5,2,'Si','Al',5.5) % 14 octahedral substitutions and 2 tetrahedral substitutions

%% Assign the Interface FF (Heinz, 2005) atomtypes to the montmorillonite atom struct
atom = interface_atom(atom,Box_dim,'interface') % Assign the Interface FF atom types to the atomstruct

%% Write a Interface FF .pdb file
write_atom_pdb(atom,Box_dim,filename_out); % Print the clay sheet to a .pdb file

%% Assign modified Interface FF atomtypes to the montmorillonite atom struct
% The original Interface FF does not contain all different atomtypes needed to
% model all sorts of clays/minerals, like for instance clay layers with 
% edges. Hence here the Interface FF forcefield is slightly modified...
% Primarily it uses different atom names (see conversion list below). 
% Secondly, it contains new oxygen atomtypes that do not exist in the 
% original Interface FF publication from Heinz et al 2005, which can be used 
% to model lets say clay edges.
% 
% Below is the list of atomtype names from the
% Heinz, 2005 paper and the ones modified here, having atomtype names that 
% simply make more sense to me. Note that I have also added a few to the 
% Heinz (Heinz, et al 2005) atomtypes, like oahe/oahhe/oshe etc.

% Interface FF from Heinz et al., 2005    = {'h*','ho','o*','oh','ob','obos','obts','obss', 'ohs', 'oas', 'oahhe','oahe', 'oshe','st','ao','at','mgo', 'mgh','cao','cah','feo','lio','Li','Na','K','Rb','Cs','Mg','Ca','Sr','Ba','F','Cl','Br','I'}';
% modified Interface FF (MHolmboe)        = {'Hw', 'H','Ow','Oh','O', 'Omg', 'Oalt','Odsub','Ohmg','Oalsi','Oalhh','Oalh','Osih','Si','Al','Alt','Mgo','Mgh','Cao','Cah','Feo','Lio','Li','Na','K','Rb','Cs','Mg','Ca','Sr','Ba','F','Cl','Br','I'}';

%% Assign the modified Interface FF atomtypes to the montmorillonite atom struct
atom = interface_atom(atom,Box_dim) % Assign the Interface FF (Heinz et al., 2005) atom types to the atomstruct

%% Heal and assign the modified Interface FF atomtypes to the montmorillonite atom struct
% In cases were atoms need healing, or in order to protonate edge groups,
% one can use a slighlt longer command like below. For more info look into
% the <interface_atom.html interface_atom> function and lines 46-77.
atom = interface_atom(atom,Box_dim,'interface','tip3p',[1:5]) % Assign the Interface FF atom types to the atomstruct

%% Write the new modified Interface FF .pdb file
write_atom_pdb(atom,Box_dim,strcat('mod_',filename_out)); % Print the clay layer to a .pdb file

