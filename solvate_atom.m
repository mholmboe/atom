%% solvate_atom.m
% * This function generates a certain region defined by <limits> with a water structure of density <density>
% * Tested 15/07/2017
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * SOL = solvate_atom(limits,density,r,maxsol,solute_atom)
% * SOL = solvate_atom(limits,density,r,maxsol,solute_atom,'tip4p')
% * SOL = solvate_atom(limits,density,r,maxsol,solute_atom,'spc_ice')


function SOL = solvate_atom(limits,density,r,maxsol,solute_atom,varargin)

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
end

% Old way of doing it...
% SOL=import_atom_gro(strcat('spc864_',num2str(density,'%3.2f'),'.gro')); % I was fastest... [2 2 1]*216.gro

if nargin==6
    watermodel=varargin(1);
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
    elseif strncmpi(watermodel,'spc',3)
        SOL=import_atom_gro('864_spc.gro');
    elseif strncmpi(watermodel,'tip3p',5)
        SOL=import_atom_gro('864_tip3p.gro');
    elseif strncmpi(watermodel,'tip4p',5)
        SOL=import_atom_gro('864_tip4p.gro');
    elseif strncmpi(watermodel,'tip5p',5)
        SOL=import_atom_gro('864_tip5p.gro');
    end
else
    SOL=import_atom_gro('864_spc.gro');
end

atomsperSOL=sum([SOL.molid]==1);

SOL=scale_atom(SOL,Box_dim,[1 1 1]./density,'all')

nx=ceil(Lx/Box_dim(1))
ny=ceil(Ly/Box_dim(2))
nz=ceil(Lz/Box_dim(3))
SOL=replicate_atom(SOL,Box_dim,[nx ny nz],'xzy');
disp('Replicated n times');
nx*ny*nz;

if (limits(1)+limits(2)+limits(3)) ~= 0
    disp('Translating the water box');
    SOL=translate_atom(SOL,[limits(1) limits(2) limits(3)],'all');
end
% assignin('caller','SOL',SOL);
% pause
disp('nSOL before merge');
size(SOL,2)/atomsperSOL;

if size(solute_atom,2) > 0
    if size(SOL,2) > 10000 || size(solute_atom,2) > 10000
        nSOL_block=size(SOL,2)/(nx*ny*nz);
        SOL_count=1;SOL_merged=[];count=1;
        while SOL_count< size(SOL,2)
            SOL_block= SOL(SOL_count:SOL_count+nSOL_block-1);
            SOL_block = merge_atom(solute_atom,limits(4:6)-.4,SOL_block,'molid','Hw',[r-.4 r]);
            SOL_merged = [SOL_merged SOL_block];
            SOL_count=SOL_count+nSOL_block;
            disp('Number of water molecules...');
            count=count+1;
            size(SOL_merged,2)/atomsperSOL
        end
        SOL=SOL_merged;
    else
        SOL = merge_atom(solute_atom,limits(4:6)-.4,SOL,'molid','Hw',[r-.4 r]);
    end
else
    SOL = slice_atom(SOL,[limits(1) limits(2) limits(3) limits(4:6)-.4],0);
end
SOL=update_atom(SOL);



% If the whole atom struct contains water
if sum(strcmpi([SOL(1:atomsperSOL:end).type],'OW')) == size(SOL,2)/atomsperSOL
        % Randomize the order of for water molecules
        nSOL=size(SOL,2);
        ind_rand=randperm(nSOL);
        ind_sel=ismember(ind_rand,1:atomsperSOL:nSOL);
        OW_ind=ind_rand(ind_sel);
        if atomsperSOL == 3
            HW1_ind=OW_ind+1;
            HW2_ind=OW_ind+2;
            ind=reshape([OW_ind;HW1_ind;HW2_ind],[],1);
            SOL=SOL(ind);
        elseif atomsperSOL == 4
            HW1_ind=OW_ind+1;
            HW2_ind=OW_ind+2;
            MW_ind=OW_ind+3;
            ind=reshape([OW_ind;HW1_ind;HW2_ind;MW_ind],[],1);
            SOL=SOL(ind);
        elseif atomsperSOL == 5
            HW1_ind=OW_ind+1;
            HW2_ind=OW_ind+2;
            LP1_ind=OW_ind+3;
            LP2_ind=OW_ind+4;
            ind=reshape([OW_ind;HW1_ind;HW2_ind;LP1_ind;LP2_ind],[],1);
            SOL=SOL(ind);
        else
            disp('Something went wrong with the randomization process...')
        end
end

% Delete water molecules if not using the <maxsol> option
if iscellstr({maxsol}) == 0
    if atomsperSOL*maxsol > size(SOL,2)
        disp('Ooops, you asked for too much water...')
        disp('Maximum number of waters allowed without changing the water density or rmin is:')
        size(SOL,2)/atomsperSOL
        SOL=SOL(1:atomsperSOL*maxsol);
    else
        SOL=SOL(1:atomsperSOL*maxsol);
    end
end

disp('nSOL after merge')
size(SOL,2)/atomsperSOL

assignin('caller','SOL',SOL);

SOL=update_atom(SOL);

assignin('caller','limits',limits);



