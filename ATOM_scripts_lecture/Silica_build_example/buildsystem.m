%% This script creates a structure file for clay+ions+alcohol. Normally one uses one or two structure files
%% with some clay layers, then one can add layers of ions and then 1-2 segments
%% of water, either on a grid or using a preequilibrated water body.
clear all;
format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename_out='preem.gro'; % Total system
Box_height=100;    % was 140 in version 1
rmin=2.3; % Nearest solute-solvent distance
nACN=1550; %3*724*10/14;% Number of ACN molecues
nSOL=170; %3*80*10/14;% Number of water molecues

%% Solute
% nMOL=4;

%% Electrolyte
nAcetate=12;
nNH4=28;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add solid slab
%% atom = import_atom(filename,translation_vector,Full_Box_dim)
Solid_slab=import_atom('silica.pdb'); % the .gro file creating after AssignClayff
Silica=center_atom(Solid_slab,Box_dim);Silica=translate_atom(Silica,[0 0 -Box_dim(3)/2]);
%% Set the full system box size
Full_Box_dim=Box_dim;Full_Box_dim(3)=Box_height;
System = wrap_atom(Silica,Full_Box_dim); % Do we want to wrap atoms into cell
% vmd(System,Full_Box_dim)

%% Add ions/molecules randomly wherever there is space (this is similar to solvate_atom())
% Ion = create_atom(Cation,Cation,limits,ceil(nCation/2),2,System);System = update_atom({System Ion}); % we are not using this function, insert_Atom better for molecules
Acetate_x1=import_atom('1xAcetate.gro');Acetate_Box_dim=Box_dim;
NH4_x1=import_atom('1xAmmonium.gro');NH4_Box_dim=Box_dim;
% MOL_x1=import_atom('mol.gro');MOL_Box_dim=Box_dim;[MOL_x1.resname]=deal({'MOL'});
%% Loop over lower and upper 1/2th of the box
for region=[1 2]
    limits=[0 0 (region-1)*1/2*Full_Box_dim(3) Full_Box_dim(1) Full_Box_dim(2) region*1/2*Full_Box_dim(3)]
    NH4=insert_atom(NH4_x1,limits,'rotate',rmin,ceil(nNH4/2),System);System = update_atom({System NH4});
end
for region=[1 2]
    limits=[0 0 (region-1)*1/2*Full_Box_dim(3) Full_Box_dim(1) Full_Box_dim(2) region*1/2*Full_Box_dim(3)]
    Acetate=insert_atom(Acetate_x1,limits,'rotate',rmin,ceil(nAcetate/2),System);System = update_atom({System Acetate});
end
% for region=[1 2]
%     limits=[0 0 (region-1)*1/2*Full_Box_dim(3) Full_Box_dim(1) Full_Box_dim(2) region*1/2*Full_Box_dim(3)]
%     MOL=insert_atom(MOL_x1,limits,'rotate',rmin/2,ceil(nMOL/2),System);System = update_atom({System MOL});
% end

%% Solvate the system
%% atom = solvate_atom(solute_atom,[1x3 or 1x6 Box vector],waterdensity,rmin,number of water or just 'max');
ACN_slab=import_atom('512xACN.gro');ACN_Box_dim=Box_dim;
ACN = solvate_atom([0 0 10 Full_Box_dim(1) Full_Box_dim(2) Full_Box_dim(3)-10],1.2,rmin,nACN,System,'custom',ACN_slab,ACN_Box_dim);System = update_atom({System ACN});
SOL = solvate_atom([0 0 15 Full_Box_dim(1) Full_Box_dim(2) Full_Box_dim(3)-15],1.1,rmin,nSOL,System,'tip4p');System = update_atom({System SOL}); % Large rmin, otherwise H2O is put into the silica slab
%% Put all atom structs together
% System = translate_atom(System,[0 0 Full_Box_dim(3)],'all'); % Do we want to translate the system
% System = center_atom(System,Full_Box_dim,{'SIL'},'z'); % Do we want to center the system
% System = wrap_atom(System,Full_Box_dim); % Do we want to wrap atoms into cell

%% Reorder the molecules in the order we want with respect to the resnames, so it matches the molecular topology file (topol.top in Gromacs)
System_reordered = reorder_atom(System,{'SIL' 'NH4' 'ACE' 'ACN' 'SOL'},'resname');

%% Write the system to a .gro structure file
write_atom_gro(System_reordered,Full_Box_dim,filename_out);

%% Plot the final structure in vmd vmd(System,Full_Box_dim) % Use VMD to
vmd(System_reordered,Full_Box_dim);
