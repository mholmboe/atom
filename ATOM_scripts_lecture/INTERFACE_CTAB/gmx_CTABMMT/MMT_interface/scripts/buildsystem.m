%% This script creates a structure file for clay+ions+water. Normally one uses one or two structure files
%% with some clay layers, then one can add layers of ions and then 1-2 segments
%% of water, either on a grid or using a preequilibrated water body.
clear all;
format compact;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clayfiles={'interface_MMT_1.gro' 'interface_MMT_2.gro'}; % the .gro file creating after AssignClayff
molfiles={'1xCTAB_16_16Ang.pdb'}; % the .gro file creating after AssignClayff
filename_out='test_preem.gro'; % Total system
d001=60*1.1;    % 12.4 15.6 18.9 21.8 24.8 for 1W   2W   3W   4W   5W
nlayers=length(clayfiles);
Full_Box_dim=[31.1880 36.0600 nlayers*d001];
nSOL=2000; %UCinX*UCinY*5*3; % 5 water molecules per unit cell and monolayer per MMT layer is reasonable 
Cation='Na';nCation=8;
System=[]; % The final stom struct for the whole system, called ?System?
%% Define regions with an 1x6 array, like [xlo ylo zlo xhi yhi zhi];
v1=[0 0   0*Full_Box_dim(3) Full_Box_dim(1) Full_Box_dim(2) 1/4*Full_Box_dim(3)];
v2=[0 0 1/4*Full_Box_dim(3) Full_Box_dim(1) Full_Box_dim(2) 2/4*Full_Box_dim(3)];
v3=[0 0 2/4*Full_Box_dim(3) Full_Box_dim(1) Full_Box_dim(2) 3/4*Full_Box_dim(3)];
v4=[0 0 3/4*Full_Box_dim(3) Full_Box_dim(1) Full_Box_dim(2) 4/4*Full_Box_dim(3)];
v_lower=[0 0   0*Full_Box_dim(3) Full_Box_dim(1) Full_Box_dim(2) 1/2*Full_Box_dim(3)];
v_upper=[0 0 1/2*Full_Box_dim(3) Full_Box_dim(1) Full_Box_dim(2) 2/2*Full_Box_dim(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add mineral lattices
%% atom = import_atom(filename,translation_vector,Full_Box_dim)
MMT1 = import_atom(clayfiles{1},[0 0 0],Full_Box_dim); System = update_atom(MMT1);
MMT2 = import_atom(clayfiles{2},[5.2/3 9/3 d001],Full_Box_dim); System = update_atom({System MMT2});
%vmd(System,Full_Box_dim)

%% Add CTAB on a grid
CTAB = import_atom(molfiles,[0 0 18]); %System = update_atom({System CTAB});
CTAB1 = replicate_atom(CTAB,Box_dim,[2 2 1]);System = update_atom({System CTAB1});
CTAB2 = rotate_atom(CTAB1,Full_Box_dim,[180 0 0]);
CTAB2 = translate_atom(CTAB2,[0 0 30],'all');System = update_atom({System CTAB2});
CTAB3 = translate_atom(CTAB1,[0 0 d001],'all');System = update_atom({System CTAB3});
CTAB4 = translate_atom(CTAB2,[0 0 d001],'all');System = update_atom({System CTAB4});
% vmd(System,Full_Box_dim);

%% Or add any number of molecules to a region specified by the v1-v4 limits
% atom = insert_atom(atom,limits,rotate_vector or 'random',rmin,num,atom_solute,{type1 type2}(optional cell),difference in mean z between typ1 type2 (optional number))
% CTAB1 = insert_atom(import_atom(molfiles),v1,'random',2,8,System,{'N' 'C'},2);System = update_atom({System CTAB1});
% CTAB2 = insert_atom(import_atom(molfiles),v2,'random',2,8,System,{'C' 'N'},2);System = update_atom({System CTAB2});
% CTAB3 = insert_atom(import_atom(molfiles),v3,'random',2,8,System,{'C' 'N'},2);System = update_atom({System CTAB3});
% CTAB4 = insert_atom(import_atom(molfiles),v4,'random',2,8,System,{'N' 'C'},2);System = update_atom({System CTAB4});
% vmd(System,Full_Box_dim);

%% Add ions randomly wherever there is space (this function is similar to solvate_atom())
%% atom = create_atom(type,resname,limits,scale,maxion,in_atom)
Ion1 = create_atom(Cation,Cation,v1,4,nCation/2,System); System = update_atom({System Ion1});
Ion2 = create_atom(Cation,Cation,v2,4,nCation/2,System); System = update_atom({System Ion2});
Ion3 = create_atom(Cation,Cation,v3,4,nCation/2,System); System = update_atom({System Ion3});
Ion4 = create_atom(Cation,Cation,v4,4,nCation/2,System); System = update_atom({System Ion4});
%vmd(System,Full_Box_dim);

%% Solvate the system
System = wrap_atom(System,Full_Box_dim); % Wrap all solute atoms into the box before adding water
%% atom = solvate_atom(solute_atom,[1x3 or 1x6 Box vector],waterdensity,rmin,number of water or just 'max');
SOL1 = solvate_atom(v_lower,1.1,2,nSOL,System); System = update_atom({System SOL1});
SOL2 = solvate_atom(v_upper,1.1,2,nSOL,System); System = update_atom({System SOL2});
% vmd(System,Full_Box_dim);

%% A final touch
System = translate_atom(System,[0 0 d001/2],'all'); % Do we want to translate the system
System = wrap_atom(System,Full_Box_dim); % Do we want to wrap atoms into cell
System = center_atom(System,Full_Box_dim,'MMT','z'); % Do we want to center the system
System = charge_atom(System,Box_dim,'interface','tip3p') % Check the total clayff charge

%% Write the system to a .gro structure file
write_atom_gro(System,Full_Box_dim,filename_out);

%% Plot the final structure in vmd
vmd(System,Full_Box_dim) % Use VMD to plot the simulation box



