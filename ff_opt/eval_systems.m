%% eval_systems.m
% * This special function is the objective function for lsqnonlin used in the
% forcefield optimization scheme called autofit_ff. This function is called by
% run_opt_ff_lsqnonlin.m and opt_ff.m
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples


function eval_systems(dirtype,varargin)

ParentSimDir = pwd;

%% Directory filename
filename=lower(dirtype);

fid = fopen('thermo_out.dat','at+');
fprintf(fid, '%s\n', ('%    GII_ref     GII_noH_ref    GII     GII_noH    Bond_diff    Angle_diff     mean_pressure   mean_potential  mean_posres     mean_temperature '));
fclose(fid);

fid = fopen('Box_out.dat','at+');
fprintf(fid, '%s\n', ('% Cell_ref   Cell_ave   Cell_std   diff_Cell_prcnt    mean_Box_X     mean_Box_Y  mean_Box_Z  diff_Box_X  diff_Box_Y   diff_Box_Z Cell '));
fclose(fid);

%% Initiate some variables
all_mean_posres=[];all_GII=[];
all_diff_BoxX=[];all_diff_BoxY=[];all_diff_BoxZ=[];

% Loop through each subfolder and extract the numeric part
subFolders = dir('OPT*');
numSimulation = []; % numSimulation=size(dir(strcat(dirtype,'*')),1);
for k = 1:length(subFolders)
    folderName = subFolders(k).name;

    % Use regular expression to extract the number from the folder name
    numStr = regexp(folderName, '\d+', 'match');

    if ~isempty(numStr)
        numSimulation(end+1) = str2double(numStr{1}); %#ok<*AGROW>
    end
end
numSimulation=sort(numSimulation);

% try
try
    for i=numSimulation %1:numSimulation

        tempdir=strcat(dirtype,num2str(i),'/');

        cd(tempdir);

        %% Set new charges
        load('system_info.mat');

        if size(ff,2)==1
            ff=ff.ff;
        end

        copyfile('topol.top','temp.top');

        try
            disp('Running gmx grompp')
            !bash run_gmx_grompp.sh
        catch
            disp('Something went wrong with gmx grompp ...')
            disp('numSimulation')
            i
        end
        if gpuDeviceCount("available") > 0 || strcmp(computer('arch'),'maca64') == 1
            %% Running in serial
            disp('Running gmx mdrun in serial')
            dirtypefiles=dir;
            if sum(strcmp({dirtypefiles.name},'opt.gro'))==0
                !bash run_gmx_mdrun.sh 'opt -ntomp 8 -rdd 0.5 -dds 0.90'
            end
            %% Running in serial
        end

        cd(ParentSimDir);
    end
catch
    disp('Did not finish gmx grompp for some reason...')
    cd(ParentSimDir);
end

if gpuDeviceCount("available") == 0 || ~strcmp(computer('arch'),'maca64') == 0
    %% -multidir simulations: requires gmx_mpi runs
    try
        if numSimulation==1 || ~strcmp(computer('arch'),'maca64') == 0
            tempdir=strcat(dirtype,num2str(i),'/');
            cd(tempdir);
            disp('Running gmx mdrun')
            !bash run_gmx_mdrun.sh 'opt -nt 4 -rdd 0.5 -dds 0.90'
            cd(ParentSimDir);
        else
            SimulationDirectories=[];
            for i=numSimulation %1:numSimulation
                SimulationDirectories=[SimulationDirectories ' ' strcat(dirtype,num2str(i))]; % Directories strin, starting with a ' ' space
            end
            disp('Running gmx_mpi mdrun, since we want to use the -multidir option!')
            eval(char(strcat('!bash ../gmx_scripts/run_gmx_mdrun_mpi.sh',{' '},num2str(numSimulation))));
        end
    catch
        disp('Did not finish gmx mdrun for some reason...')
        disp('numSimulation')
        i
        cd(ParentSimDir);
    end
    %% -multidir simulations: requires gmx_mpi runs
end

for i=numSimulation %1:numSimulation

    tempdir=strcat(dirtype,num2str(i),'/');

    cd(tempdir);

    try
        disp('Analyzing Box size!')
        !bash run_gmx_boxsize.sh # Note the different index order for zsh compared to bash

        Box_X=import_xvg('Box_X.xvg');
        analyze_start=floor(0.5*size(Box_X,1));
        init_Box_X=Box_X(1,2); % Since we want to compare with the Boxsize at t=0
        mean_Box_X=mean(Box_X(analyze_start:end,2));
        diff_Box_X=abs(mean(Box_X(:,2)-init_Box_X));
        all_diff_BoxX=[all_diff_BoxX diff_Box_X/mean_Box_X];

        Box_Y=import_xvg('Box_Y.xvg');
        init_Box_Y=Box_Y(1,2);
        mean_Box_Y=mean(Box_Y(analyze_start:end,2));
        diff_Box_Y=abs(mean(Box_Y(:,2)-init_Box_Y));
        all_diff_BoxY=[all_diff_BoxY diff_Box_Y/mean_Box_Y];

        Box_Z=import_xvg('Box_Z.xvg');
        init_Box_Z=Box_Z(1,2);
        mean_Box_Z=mean(Box_Z(analyze_start:end,2));
        diff_Box_Z=abs(mean(Box_Z(:,2)-init_Box_Z));
        all_diff_BoxZ=[all_diff_BoxZ diff_Box_Z/mean_Box_Z];
    catch
        disp('Could not extract the Box_X/Y/Z from the .edr file')
        pause(1)
        mean_Box_X=0;
        mean_Box_Y=0;
        mean_Box_Z=0;
        diff_Box_X=0;
        diff_Box_Y=0;
        diff_Box_Z=0;
        all_diff_BoxX=0;
        all_diff_BoxY=0;
        all_diff_BoxZ=0;
    end

    try
        !bash run_gmx_thermo.sh # Note the different index order for zsh compared to bash
        disp('Analyzing pressure!')

        mean_pressure=load('pressure.dat'); %mean(pressure(:,2));

        disp('Analyzing potential energy!')

        mean_potential=load('potential.dat'); %median(potential(:,2));

        disp('Analyzing position restraints!')

        mean_posres=load('posres.dat'); %median(posres(:,2));
        all_mean_posres=[all_mean_posres mean_posres];

        disp('Analyzing temperature!')

        mean_temperature=load('temperature.dat'); %median(temperature(:,2));

    catch
        disp('Could not extract the Pressure, Potential, Temperature.. from the .edr file')
        pause(1)
        load('prev_res.mat')
        mean_pressure=2*mean_pressure;
        mean_potential=2*mean_potential
        mean_posres=2*mean_posres
        all_mean_posres=[all_mean_posres 2*mean_posres];
        mean_temperature=mean_temperature;
    end

    %% Analyzing the trajectory
    !bash run_gmx_traj.sh

    if exist(fullfile(pwd,'prev_res.mat'),'file') ~= 2
        try
            copyfile('opt.gro','opt_t0.gro');
            % copyfile('temp.top',strcat('../first',num2str(i),'.top'));
            % copyfile('temp.top',strcat('first',num2str(i),'.top'));
            % copyfile('rmsf.xvg','rmsf_t0.xvg');
            % copyfile('rmsdev.xvg','rmsdev_t0.xvg');
            % copyfile('rmsdev.pdb','rmsdev_t0.pdb');
        catch
            disp('Could not save the inital simulation output!')
        end
    end
    save('prev_res.mat');
    cd(ParentSimDir);
end

indexes=[]; % Matlab struct
for i=numSimulation %1:numSimulation

    tempdir=strcat(dirtype,num2str(i),'/');

    cd(tempdir);

    disp('Analyzing bond distances/valences!')
    load('prev_res.mat');pause(0.1)
    load('system_info.mat');pause(0.1)
    indexes(i).ind=find(strncmpi([atom.resname],resname,2));
    preopt_prop = analyze_atom(element_atom(atom([indexes(i).ind])),Box_dim,0,target_Bond_index);
    preopt_prop = analyze_atom(element_atom(atom([indexes(i).ind])),Box_dim,0,target_Bond_index,[preopt_prop.oxstate]);
    preopt_GII=GII;
    preopt_GII_noH=GII_noH;
    preopt_Cell=Box_dim2Cell(Box_dim);

    if (max([atom([indexes(i).ind]).x])-min([atom([indexes(i).ind]).x]))>0.9*Box_dim(1) || ...
            (max([atom([indexes(i).ind]).y])-min([atom([indexes(i).ind]).y]))>0.9*Box_dim(2) || ...
            (max([atom([indexes(i).ind]).z])-min([atom([indexes(i).ind]).z]))>0.9*Box_dim(3)
        PBC=1
        [atom_ave,Box_dim_ave,Box_dim_std,Cell_ave,Cell_std]=write_ave_gro('conf','median');
    else
        PBC=0;
        [atom_ave,Box_dim_ave,Box_dim_std,Cell_ave,Cell_std]=write_ave_gro('conf','ave');
    end

    %% Since we are only looking at a fractional system, since else too many angles, bonds..
    atom_ave=atom_ave([indexes(i).ind]); % Also Box_dim_ave,Cell_ave,Box_dim_std,Cell_std...

    atom_ave = recalc_bond_atom(atom_ave,Box_dim_ave,target_Bond_index,target_Angle_index);%,target_Dihedral_index); % Will output new Bond_index, Angle_index, dihedral_index
    temp_Bond_index=Bond_index;temp_Angle_index=Angle_index;
    GII=0;
    temp_prop = analyze_atom(element_atom(atom_ave([indexes(i).ind])),Box_dim_ave,0,Bond_index,[preopt_prop.oxstate]);
    all_GII=[all_GII GII];

    dlmwrite('bond.ndx',{'[ B1B2 ]'},'delimiter','')
    dlmwrite('bond.ndx',target_Bond_index(:,1:2),'delimiter',' ','-append')

    dlmwrite('angle.ndx',{'[ A2A1A3 ]'},'delimiter','')
    dlmwrite('angle.ndx',target_Angle_index(:,1:3),'delimiter',' ','-append')


    %% To remove bonds and angles with H's
    ind_H=find(strncmpi([atom_ave.type],{'H'},1));
    [H_row,H_col]=ind2sub(size(temp_Bond_index),find(ismember(temp_Bond_index,ind_H)));
    temp_Bond_index(H_row,:)=[];
    target_Bond_index(H_row,:)=[];

    [H_row,H_col]=ind2sub(size(temp_Angle_index),find(ismember(temp_Angle_index,ind_H)));
    temp_Angle_index(H_row,:)=[];
    target_Angle_index(H_row,:)=[];

    Bond_diff=[0];Angle_diff=[0];RDF_diff=[0];CN_diff=[0];CoordNum=0;


    Bond_diff=double(temp_Bond_index(:,3))-double(target_Bond_index(:,3));
    Angle_diff=double(temp_Angle_index(:,4))-double(target_Angle_index(:,4));

    dlmwrite('../thermo_out.dat',[preopt_GII preopt_GII_noH GII GII_noH mean(abs(Bond_diff)) mean(abs(Angle_diff)) mean_pressure mean_potential mean_posres mean_temperature ],'precision','%1.6E','delimiter','\t','-append')
    dlmwrite('../Box_out.dat',[preopt_Cell Cell_ave Cell_std (Cell_ave-preopt_Cell)./Cell_ave mean_Box_X mean_Box_Y mean_Box_Z diff_Box_X/mean_Box_X diff_Box_Y/mean_Box_Y diff_Box_Z/mean_Box_Z],'precision','%1.6E','delimiter','\t','-append')

    save('prev_res.mat')

    !rm -f \#*
    !rm -f conf*.gro

    cd(ParentSimDir);

    save('prev_run.mat')

end


end
