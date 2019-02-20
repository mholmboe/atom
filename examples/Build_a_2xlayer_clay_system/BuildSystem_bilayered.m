%% This script creates a structure file for clay+ions+alcohol. Normally one uses one or two structure files
%% with some clay layers, then one can add layers of ions and then 1-2 segments
%% of water, either on a grid or using a preequilibrated water body.
clear all;
format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clayfiles={'MMT_clayff_1.gro' 'MMT_clayff_2.gro'}; % the .gro file creating after AssignClayff
filename_out='preem.gro'; % Total system
d001=45;    % 12.4 15.6 18.9 21.8 24.8 for 1W   2W   3W   4W   5W
rmin=2;
water_density=1.1;
nlayers=2;
nSOL=1440;%UCinX*UCinY*5*3; % 5 water molecules per unit cell and monolayer per MMT layer is reasonable 
Cation='Na';nCation=16; % Per interlayer
spacer=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add mineral lattices
%% atom = import_atom(filename,translation_vector,Full_Box_dim)
MMT1 = import_atom(clayfiles{1},[0 0 0]);
Full_Box_dim=Box_dim(1:3);Full_Box_dim(3)=nlayers*d001;
MMT2 = import_atom(clayfiles{2},[4*5.2/3 4*9/3 d001],Full_Box_dim);
System = update_atom({MMT1 MMT2});
% vmd(System,Full_Box_dim)

%% Add ions randomly wherever there is space (this function is similar to solvate_atom())
%% atom = create_atom(type,resname,limits,scale,maxion,in_atom)  
Ion1 = create_atom(Cation,Cation,[0 0 spacer Full_Box_dim(1) Full_Box_dim(2) spacer],nCation/2,1,System);
Ion2 = create_atom(Cation,Cation,[0 0 d001-spacer Full_Box_dim(1) Full_Box_dim(2) d001-spacer],nCation/2,1,System);
Ion3 = create_atom(Cation,Cation,[0 0 d001+spacer Full_Box_dim(1) Full_Box_dim(2) d001+spacer],nCation/2,1,System);
Ion4 = create_atom(Cation,Cation,[0 0 2*d001-spacer Full_Box_dim(1) Full_Box_dim(2) 2*d001-spacer],nCation/2,1,System);
System = update_atom({System Ion1 Ion2 Ion3 Ion4});
System = wrap_atom(System,Full_Box_dim); % Do we want to wrap atoms into cell
%vmd(System,Full_Box_dim);

%% Solvate the system1
%% atom = solvate_atom(solute_atom,[1x3 or 1x6 Box vector],waterdensity,rmin,number of water or just 'max');
SOL1 = solvate_atom([0 0 0 Full_Box_dim(1) Full_Box_dim(2) d001],water_density,rmin,nSOL,System);
SOL2 = solvate_atom([0 0 0+d001 Full_Box_dim(1) Full_Box_dim(2) 2*d001],water_density,rmin,nSOL,System);
System = update_atom({System SOL1 SOL2});
vmd(System,Full_Box_dim);

%% Put all atom structs together
% System = translate_atom(System,[0 0 -spacer/2],'all'); % Do we want to translate the system
% System = center_atom(System,Full_Box_dim,{'PYR'},'z'); % Do we want to center the system
System = wrap_atom(System,Full_Box_dim); % Do we want to wrap atoms into cell

%% Write the system to a .gro structure file
write_atom_gro(System,Full_Box_dim,filename_out);

%% Plot the final structure in vmd
vmd(System,Full_Box_dim) % Use VMD to plot the simulation box



