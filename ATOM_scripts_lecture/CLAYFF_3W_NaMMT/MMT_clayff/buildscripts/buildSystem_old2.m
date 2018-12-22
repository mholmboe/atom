%% This script creates a .xyz file for clay+ions+water. Normally one uses one or two .xyz files
%% with some clayff clay layers, then one can add 1-4 layers of ions and then 1-2 segments
%% of water, either on a grid or using a preequilibrated water body. If one wants to add
%% water first and then the ions, one have to add the water on a grid.

clear all;
format long;

%%% Your path to VMD and GMX, optional %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PATH2VMD = '/Applications/VMD\ 1.9.1.app/Contents/MacOS/startup.command'; % Add your own path to VMD here
PATH2GMX ='/usr/local/gromacs-5.02/bin/'; % Add your own path to GROMACS here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

files=[{'1_MMT.gro'} {'2_MMT.gro'}];
nlayers=length(files);
filename_out='preem.gro';

d001=18.9  % 12.4 15.6 18.9 21.8 24.8
nSOL=360; % 2750 3900 5100 6100 7100
nCatIon=16; nCatIonblk=0; nAnionblk=0; % per MMT layer

%% Add Clayff lattices
atom_Tot=[];
for i=1:nlayers
    import_atom_gro(char(files(i))); atom = translate_atom_func(atom,[i*5.1488/2 i*8.9979/2 -median([atom.z])+(i-1)*d001],'all'); atom_Tot = add2atom(XYZ_labels,XYZ_data,Box_dim,'MMT',atom_Tot);%atom_Tot=atom;
end

%% Add ions on grids
for i=1:nlayers
    ions2xyz('Na',nCatIon,'Cl',0,[Box_dim(1) Box_dim(2) (i-0.5)*d001],'xy');    atom_Tot = add2atom(XYZ_labels,XYZ_data,Box_dim,'Na',atom_Tot);
end
% ions2xyz('Na',nCatIon,'Cl',0,[Box_dim(1) Box_dim(2) 3/2*d001],'xy');    atom_Tot = add2atom(XYZ_labels,XYZ_data,Box_dim,'Na',atom_Tot);
% ions2xyz('Na',nCatIonblk/3,'Cl',0,[40 70 2*d001],'xz');                        atom_Tot = add2atom(XYZ_labels,XYZ_data,Box_dim,'Na',atom_Tot);
% ions2xyz('Na',nCatIonblk/3,'Cl',0,[40 90  2*d001],'xz');                       atom_Tot = add2atom(XYZ_labels,XYZ_data,Box_dim,'Na',atom_Tot);
% ions2xyz('Na',nCatIonblk/3,'Cl',0,[40 110  2*d001],'xz');                      atom_Tot = add2atom(XYZ_labels,XYZ_data,Box_dim,'Na',atom_Tot);
% ions2xyz('Cl',ceil(nAnionblk/2),'Cl',0,[40 100 2*d001],'xz');                 atom_Tot = add2atom(XYZ_labels,XYZ_data,Box_dim,'Cl',atom_Tot);
% ions2xyz('Cl',floor(nAnionblk/2),'Cl',0,[40 80  2*d001],'xz');                atom_Tot = add2atom(XYZ_labels,XYZ_data,Box_dim,'Cl',atom_Tot);
% ions2xyz('Na',4,'Cl',0,[40 90 2*d001],'xz');                        atom_Tot = add2atom(XYZ_labels,XYZ_data,Box_dim,'Na',atom_Tot);
% ions2xyz('Cl',4,'Cl',0,[40 95 2*d001],'xz');                        atom_Tot = add2atom(XYZ_labels,XYZ_data,Box_dim,'Cl',atom_Tot);
% ions2xyz('Ca',5,'Cl',0,[40 110  2*d001],'xz');                      atom_Tot = add2atom(XYZ_labels,XYZ_data,Box_dim,'Ca',atom_Tot);
% ions2xyz('Cl',4,'Cl',0,[40 100 2*d001],'xz');                 atom_Tot = add2atom(XYZ_labels,XYZ_data,Box_dim,'Cl',atom_Tot);
% ions2xyz('Cl',15,'Cl',0,[40 80  2*d001],'xz');                atom_Tot = add2atom(XYZ_labels,XYZ_data,Box_dim,'Cl',atom_Tot);

%% Set system box size
Box_dim = [Box_dim(1) Box_dim(2) length(files)*d001];
%% Wrap atoms into cell
atom_Tot = wrap_atom_func(atom_Tot,Box_dim);
%% Center the system
atom_Tot = center_atom_func(atom_Tot,Box_dim,'MMT','xyz');
%% Solvate the system
for i=1:nlayers
    sol2xyz(atom_Tot,[0 Box_dim(1) 0 Box_dim(2) (i-1)*d001 i*d001],nSOL,1.6,1,0); atom_Tot = add2atom(XYZ_labels,XYZ_data,Box_dim,'SOL',atom_Tot);
end
%% Check the charge after AssignClayff.m
Total_charge = check_clayff_charge(atom_Tot)
% %% Check the water density
H2Odens = check_clayff_H2Odens(atom_Tot,Box_dim)
%% Write the system
write_atom_gro(atom_Tot,Box_dim,filename_out);

%%%%%% Plot the structure %%%%%%%%%%%
system(strcat(char({PATH2VMD}),char(strcat({' '},{filename_out})))); % Use VMD to plot the simulation box   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

