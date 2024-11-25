%% solvate_atom.m
% * This function generates a certain region defined by <limits> with a
% * solvent structure of density <density>, or a solvent shell with
% * thickness shell10|15|20|25|30 Ångström around the any given solute
%
%% Function arguments
% * [limits] is a 1x3 or 1x6 volume variable
% * density (number) is the density
% * r is a number and the minimum distance between solvent-solute particles.
% * maxsol (number) maxsol is the max number of solvent molecules, or a
% * string 'maxsol' (allowing maximal solvation), or a string like shell10|15|20|25|30
% * indicating a solvatin shell around the solute.
% * solute_atom is the existing atom struct which needs solvation
% * Optional string can be the desired water model, like 'SPC' or 'TIP4P',
% or with 'custom' be used to solvate with a custom structure, like a
% ethanol slab
%
%% Dependencies
% * import_atom_gro
% * scale_atom
% * replicate_atom
% * translate_atom
% * merge_atom
% * slice_atom
% * update_atom
% * distance_matrix
% * cell_list_distance_matrix
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * SOL = solvate_atom(limits,density,r,maxsol) % Basic input arguments
% * SOL = solvate_atom(limits,density,r,'maxsol',solute_atom) % Will maximize the number of solvent molecules
% * SOL = solvate_atom(limits,density,r,maxsol,solute_atom) % Will account for existing solute sites
% * SOL = solvate_atom(limits,density,r,'shell15',solute_atom) % Will solvatize a 15Å shell around the sites in the solute_atom
% * SOL = solvate_atom(limits,density,r,maxsol,solute_atom,'tip4p') % Will use the tip4p water model
% * SOL = solvate_atom(limits,density,r,maxsol,solute_atom,'spc_ice') % Will use an hexagonal ice structure
% * SOL = solvate_atom(limits,density,r,maxsol,solute_atom,'custom',mysolvent,mysolvent_Box_dim) % mysolvent(_Box_dim) is an atom struct
%
function SOL = solvate_atom(limits,density,r,maxsol,varargin)

% Solvent shell thickness
if strcmpi(maxsol,'shell10')
    shellthickness=10;
elseif strcmpi(maxsol,'shell15')
    shellthickness=15;
elseif strcmpi(maxsol,'shell20')
    shellthickness=20;
elseif strcmpi(maxsol,'shell25')
    shellthickness=25;
elseif strcmpi(maxsol,'shell30')
    shellthickness=30;
elseif strcmpi(maxsol,'shell5')
    shellthickness=5;
elseif strcmpi(maxsol,'shell4')
    shellthickness=4;
elseif strcmpi(maxsol,'shell3')
    shellthickness=3;
elseif strcmpi(maxsol,'shell2')
    shellthickness=2;
elseif strncmpi(maxsol,'shell',5)
    shellthickness=10;
end

if numel(limits)==1
    limits(4:6)=limits(1);
    limits(1:3)=0;
    Lx=limits(4);
    Ly=limits(5);
    Lz=limits(6);
elseif numel(limits)==3
    Lx=limits(1);
    Ly=limits(2);
    Lz=limits(3);
    limits(4)=limits(1);
    limits(5)=limits(2);
    limits(6)=limits(3);
    limits(1:3)=0;
elseif numel(limits)==6
    Lx=limits(4)-limits(1);
    Ly=limits(5)-limits(2);
    Lz=limits(6)-limits(3);
elseif numel(limits)==9
    Lx=limits(1);
    Ly=limits(2);
    Lz=limits(3);
    limits(4)=limits(1);
    limits(5)=limits(2);
    limits(6)=limits(3);
    limits(1:3)=0;
end

% Old way of doing it...
% SOL=import_atom_gro(strcat('spc864_',num2str(density,'%3.2f'),'.gro')); % I was fastest... [2 2 1]*216.gro

if nargin == 4
    solute_atom=[];
else
    solute_atom=varargin{1};
end

if nargin>5
    watermodel=varargin(2);
    if iscell(watermodel)
        assignin('caller','watermodel',watermodel)
        watermodel=char(watermodel{1});
    end
    if strncmpi(watermodel,'spc_ice',7)
        SOL=import_atom_gro('96spc_hex_ice_h.gro');
        pause(1)
        disp('Adding hexagonal spc ice!!!')
    elseif strncmpi(watermodel,'tip4p_ice',9)
        SOL=import_atom_gro('96tip4p_hex_ice_h.gro');
        pause(1)
        disp('Adding hexagonal tip4p ice!!!')
    elseif strncmpi(watermodel,'spce',4)
        SOL=import_atom_gro('864_spce.gro');
        disp('Adding spc/e!!!')
    elseif strncmpi(watermodel,'spc',3)
        SOL=import_atom_gro('864_spc.gro');
        disp('Adding spc!!!')
    elseif strncmpi(watermodel,'tip3p',5)
        SOL=import_atom_gro('864_tip3p.gro');
        disp('Adding tip3p!!!')
    elseif strncmpi(watermodel,'tip4p',5)
        SOL=import_atom_gro('864_tip4p.gro');
        disp('Adding tip4p!!!')
    elseif strncmpi(watermodel,'tip5p',5)
        SOL=import_atom_gro('864_tip5p.gro');
        disp('Adding tip5p!!!')
    elseif strncmpi(watermodel,'swm4',4)
        SOL=import_atom_gro('864_swm4_ndp.gro');
        disp('Adding swm4_ndp!!!')
    elseif strncmpi(watermodel,'custom',5)
        disp('You are using your own solvent, you must be brave')
        SOL=varargin{3};
        if nargin>7
            Box_dim=varargin{4};
            if numel(Box_dim)>3
                disp('Custom solvent boxes must be orthogonal, i.e. numel(Box_dim)==3')
                pause
            end
        end
    end
else
    SOL=import_atom_gro('864_spc.gro');
    disp('Adding spc / spc/e!!!')
end

atomsperSOL=sum([SOL.molid]==1);

SOL=scale_atom(SOL,Box_dim,[1 1 1]./density,'all');

nx=ceil(Lx/Box_dim(1))
ny=ceil(Ly/Box_dim(2))
nz=ceil(Lz/Box_dim(3))
SOL=replicate_atom(SOL,Box_dim,[nx ny nz],'xzy');
disp('Replicated n times');
nx*ny*nz;

if (limits(1)+limits(2)+limits(3)) ~= 0
    disp('Translating the solvent box');
    SOL=translate_atom(SOL,[limits(1) limits(2) limits(3)],'all');
end

disp('nSOL before merge');
size(SOL,2)/atomsperSOL

if size(solute_atom,2) > 0
    if size(SOL,2) > 10000 || size(solute_atom,2) > 10000
        nSOL_block=size(SOL,2)/(nx*ny*nz);
        SOL_count=1;SOL_merged=[];count=1;
        while SOL_count<size(SOL,2)
            SOL_block= SOL(SOL_count:SOL_count+nSOL_block-1);
            SOL_block = merge_atom(solute_atom,limits(4:6)-.2,SOL_block,'molid','H',[r-.4 r]); % Can shell be implemented here instead?
            SOL_merged = [SOL_merged SOL_block];
            SOL_count=SOL_count+nSOL_block;
            disp('Number of solvent molecules...');
            count=count+1;
            size(SOL_merged,2)/atomsperSOL
        end
        SOL=SOL_merged;
    else
        SOL = merge_atom(solute_atom,limits(4:6)-.2,SOL,'molid','Hw',[r-.4 r]); % Can shell be implemented here instead?
    end
else
    SOL=slice_atom(SOL,[limits(1) limits(2) limits(3) limits(4:6)-.4],0);
end

SOL=update_atom(SOL);

% Randomize the order of the SOL molecules, and check for 'max' or 'shell' option
nSOL=size(SOL,2);
if iscellstr({maxsol}) == 1
    if strncmpi(maxsol,'max',3)
        maxsol=nSOL/atomsperSOL;
    elseif strncmpi(maxsol,'shell',5)
        disp('Will solvate a shell around the solute molecule')
        nSolute=size(solute_atom,2);
        if (size(SOL,2)+size(solute_atom,2))>50000
            dist_matrix = cell_list_dist_matrix_atom(update_atom({solute_atom SOL}),Box_dim);
            dist_matrix(dist_matrix==0)=2*shellthickness; % Set some high dummy value
        else
            dist_matrix = dist_matrix_atom(update_atom({solute_atom SOL}),Box_dim);
        end
        ind=[];
        for i=1:size(SOL,2) % Vectorize this!!!
            if min(dist_matrix(nSolute+i,1:nSolute))>shellthickness
                ind=[ind i];
            end
        end
        molid=unique([SOL(intersect([SOL.index],ind)).molid]);
        rm_ind = ismember([SOL.molid],molid);
        SOL(rm_ind)=[];
        nSOL=size(SOL,2);
        maxsol=size(SOL,2)/atomsperSOL;
        SOL=update_atom(SOL);
    end
end
rand_molid=randperm(nSOL/atomsperSOL);
if maxsol>size(rand_molid,2)
    disp('You have tried to add too many solvent molecules for the given region')
end
disp('Number of maximum solvent molecules possible')
size(rand_molid,2)
disp('Number of solvent molecules requested')
maxsol

rand_molid=rand_molid(1:maxsol);

rand_index=ismember([SOL.molid],rand_molid);
SOL=SOL(rand_index);

% Delete water molecules if not using the <maxsol> option
if iscellstr({maxsol}) == 0
    if atomsperSOL*maxsol > size(SOL,2)
        disp('Ooops, you asked for too much solvent...')
        maxsol
        disp('Maximum number of solvent molecules allowed without changing the density or rmin is:')
        size(SOL,2)/atomsperSOL
        SOL=SOL(1:atomsperSOL*maxsol);
    else
        % Use randperm instead... ?
        SOL=SOL(1:atomsperSOL*maxsol);
    end
end
SOL=update_atom(SOL);

disp('nSOL after merge')
size(SOL,2)/atomsperSOL

assignin('caller','SOL',SOL);
assignin('caller','limits',limits);

end
