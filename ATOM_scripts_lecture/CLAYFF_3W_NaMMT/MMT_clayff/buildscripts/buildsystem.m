%% This script creates a structure file for clay+ions+water. Normally one uses one or two structure files
%% with some clay layers, then one can add layers of ions and then 1-2 segments
%% of water, either on a grid or using a preequilibrated water body.
clear all;
format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename_out='preem.gro'; % Total system
d001=18.9*1.1;    % 12.4 15.6 18.9 21.8 24.8 for 1W   2W   3W   4W   5W
nlayers=2;
Full_Box_dim=[30.9600   35.8633  nlayers*d001];
nSOL=360;%UCinX*UCinY*5*3; % 5 water molecules per unit cell and monolayer per MMT layer is reasonable 
Cation='Na';nCation=16;
System=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add mineral lattices
%% atom = import_atom(filename,translation_vector,Full_Box_dim)
MMT1 = import_atom('MMT_clayff_1.gro',[0 0 0],Full_Box_dim); System = update_atom({System MMT1});
MMT2 = import_atom('MMT_clayff_2.gro',[5.2/3 9/3 d001],Full_Box_dim); System = update_atom({System MMT2});
% vmd(System,Full_Box_dim)

%% Add ions using the copy_atom function
%% atom = copy_atom(atom,origtype,newtype,newresname,trans_vec,num);
% Ion1 = copy_type(MMT1,'Mgo',Cation,Cation,[0 0 d001/2],sum(ismember([MMT1.type],'Mgo'))); System = update_atom({System Ion1});
% Ion2 = copy_type(MMT2,'Mgo',Cation,Cation,[0 0 d001/2],sum(ismember([MMT2.type],'Mgo'))); System = update_atom({System Ion2});
% vmd([MMT1 MMT2 Ion1 Ion2],Full_Box_dim);

%% Add ions randomly wherever there is space (this function is similar to solvate_atom())
%% atom = create_atom(type,resname,limits,scale,maxion,in_atom)
Ion1 = create_atom(Cation,Cation,[0 0 0 Full_Box_dim(1) Full_Box_dim(2) Full_Box_dim(3)/2],nCation,3.5,System); System = update_atom({System Ion1});
Ion2 = create_atom(Cation,Cation,[0 0 Full_Box_dim(3)/2 Full_Box_dim(1) Full_Box_dim(2) Full_Box_dim(3)],nCation,3.5,System); System = update_atom({System Ion2});
%vmd(System,Full_Box_dim);

%% Solvate the system
System = wrap_atom(System,Full_Box_dim); % Wrap all solute atoms into the box before adding water
%% atom = solvate_atom([1x3 or 1x6 Box vector],waterdensity,rmin,number of water or just 'max',solute_atom);
SOL1 = solvate_atom([0 0 0 Full_Box_dim(1) Full_Box_dim(2) d001],1.0,1.9,nSOL,System); System = update_atom({System SOL1});
SOL2 = solvate_atom([0 0 d001 Full_Box_dim(1) Full_Box_dim(2) Full_Box_dim(3)],1.0,1.9,nSOL,System); System = update_atom({System SOL2});
% vmd(System,Full_Box_dim);

%% Put all atom structs together
% System = translate_atom(System,[0 0 d001/2],'all'); % Do we want to translate the system
% System = center_atom(System,Full_Box_dim,{'MMT'},'z'); % Do we want to center the system
System = wrap_atom(System,Full_Box_dim); % Do we want to wrap atoms into cell
System = charge_atom(System,Box_dim,'clayff','spc/e');

%% Write the system to a .gro structure file
write_atom_gro(System,Full_Box_dim,filename_out);

write_atom_lmp(System,Full_Box_dim,filename,1.25,1.25,'clayff','spce')

%% Plot the final structure with plot_atom() or show_atom() or vmd()
plot_atom(System,Full_Box_dim,.1)
% show_atom(System,Full_Box_dim)
% vmd(System,Full_Box_dim) % Use VMD to plot the simulation box

% write_atom_all(System,Full_Box_dim,'system',1.25,2.25,'clayff','spc/e');



