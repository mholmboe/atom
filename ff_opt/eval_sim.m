%% eval_sim.m
% * This special function is the objective function for lsqnonlin used in the
% forcefield optimization scheme called autofit_ff. This function is called by
% run_opt_ff_lsqnonlin.m and opt_ff.m
%
%% Version
% 3.00
%
% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # f=@(x)eval_sim(x,scalefactors,dirtype,indexes);  % From run_opt_ff_lsqnonlin
% # res = eval_sim(fx,scalefactors,dirtype,indexes); % f1,f2,f3); % etc to run with extra parameters..

function res = eval_sim(param,scalefactors,dirtype,varargin)

ParentSimDir = pwd;

% Scale back all the parameters
param=param./scalefactors;

%% To also fit the charges or other extra parameters
fextra=nargin-3;
nparam=size(param,2);

disp('Will run with these parameters:')
param

%% PBC variable
PBC=1; % Assume periodic molecules at first, check later in system_info.mat

%% Directory filename
filename=lower(dirtype);

%% Initiate some variables
all_mean_posres=[];
all_diff_BoxX=[];all_diff_BoxY=[];all_diff_BoxZ=[];
skip_simulation=[];

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

try
    try
        for i=numSimulation %1:numSimulation

            tempdir=strcat(dirtype,num2str(i),'/');

            cd(tempdir);

            %% Set new charges
            load('system_info.mat');
            nOpt_labels=size(Opt_labels,2);

            if size(ff,2)==1
                ff=ff.ff;
            end

            if nparam>2*nOpt_labels
                disp('Setting new charges in the .itp topology file!')
                itp=import_itp('min.itp');
                temp_atom=atom;
                [temp_atom.type]=atom.fftype;
                if ~strncmp(Opt_labels,'H',1)
                    try
                        Hind=strcmp([ff.type],'H');
                        Hq=ff(Hind).charge;
                    catch
                        Hq=0.4;
                    end
                    temp_atom=charge_clayff_atom(temp_atom,Box_dim,[Opt_labels(~strncmp(Opt_labels,'O',1)) {'H'}],[param(end-(nparam-2*nOpt_labels)+1:end) Hq]);
                else
                    temp_atom=charge_clayff_atom(temp_atom,Box_dim,[Opt_labels(~strncmp(Opt_labels,'O',1))],[param(end-(nparam-2*nOpt_labels)+1:end)]);
                end
                [itp.atoms.charge]=[temp_atom.charge]';
                write_itp(itp,'min.itp');
            end

            print_top(param,'header.top','topol.top');
            change_top(param,'header.top');

            copyfile('temp.top',strcat('../last',num2str(i),'.top'));
            copyfile('temp.top',strcat('last',num2str(i),'.top'));

            dirtypefiles=dir;
            if sum(strcmp({dirtypefiles.name},strcat('prev_analysis_',num2str(i),'.mat')))==0
                iteration=0;
                Opt_ind=[];
                orig_param=param;
                prev_param=param;
                save('orig_param.mat','orig_param');
                save('prev_param.mat','prev_param');
                fid = fopen('thermo_out.dat','at+');
                fprintf(fid, '%s\n', ('% Sim   norm(res_temp)^2     GII     RDF_diff    CN_diff    CoordNum    Bond_diff    Angle_diff     mean_pressure   mean_potential  mean_posres     mean_temperature    param    mean_Box_X     mean_Box_Y  mean_Box_Z  diff_Box_X diff_Box_Y   diff_Box_Z Cell_ave Cell_std'));
                fclose(fid);
            else
                load('system_info.mat');
                load('prev_param.mat');
                Opt_labels_list=[Opt_labels;Opt_labels];
                Opt_labels_list=Opt_labels_list(:)';
                Opt_ind=find(prev_param-param);
                if numel(Opt_ind)>0
                    if sum(ismember(Atom_labels,Opt_labels_list(Opt_ind)))==0
                        skip_simulation=[skip_simulation i];
                        cd(ParentSimDir);
                        continue
                    end
                end
            end

            try
                disp('Running gmx grompp')
                % eval(char(strcat('!gmx grompp -f',{' '},filename,'.mdp -c pre',filename,'.gro -r pre',filename,'.gro -p temp.top  -n index.ndx -pp temp_gmx -po temp_gmx -maxwarn 3 -o',{' '},filename,'.tpr ')));% 2>/dev/null')));
                !bash run_gmx_grompp.sh
            catch
                disp('Something went wrong with gmx grompp ...')
                disp('numSimulation')
                i
                pause
            end
            if gpuDeviceCount("available") > 0 || strcmp(computer('arch'),'maca64') == 1
                %% Running in serial
                disp('Running gmx mdrun in serial')
                % eval(char(strcat('!gmx mdrun -nt 4 -deffnm',{' '},filename)));
                !bash run_gmx_mdrun.sh 'opt -ntomp 8 -rdd 0.5 -dds 0.90'
                %% Running in serial
            end
            cd(ParentSimDir);
        end
    catch
        disp('Did not finish gmx grompp for some reason...')
        cd(ParentSimDir);
    end

    if gpuDeviceCount("available") == 0 && strcmp(computer('arch'),'maca64') == 0
        %% -multidir simulations: requires gmx_mpi runs
        try
            if numSimulation==1 || ~strcmp(computer('arch'),'maca64') == 0
                tempdir=strcat(dirtype,num2str(i),'/');
                cd(tempdir);
                disp('Running gmx mdrun')
                % eval(char(strcat('!gmx mdrun -nt 4 -deffnm',{' '},filename)));
                !bash run_gmx_mdrun.sh 'opt -nt 4 -rdd 0.5 -dds 0.90'
                cd(ParentSimDir);
            else
                SimulationDirectories=[];
                for i=numSimulation
                    SimulationDirectories=[SimulationDirectories ' ' strcat(dirtype,num2str(i))]; % Directories strin, starting with a ' ' space
                end
                disp('Running gmx_mpi mdrun, since we want to use the -multidir option!')
                %           eval(char(strcat('!mpirun -np',{' '},num2str(numSimulation),{' '},'gmx_mpi mdrun -ntomp 1 -deffnm',{' '},filename,{' '},'-multidir',SimulationDirectories))); % -nt 6
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

        %         if ismember(i,skip_simulation)
        %             continue
        %         end

        try
            try
                atom=import_atom_gro('opt.gro');
            catch
                disp('Simulation did not finish properly..')
            end
            disp('Analyzing Box size!')
            !bash run_gmx_boxsize.sh # Note the different index order for zsh compared to bash
            %             % eval(char(strcat('!echo Box-X | gmx energy -f',{' '},filename,'.edr -b 0 -o Box_X.xvg 2>/dev/null')));
            %             !bash run_gmx_energy.sh Box-X
            Box_X=import_xvg('Box_X.xvg');
            analyze_start=floor(0.5*size(Box_X,1));
            init_Box_X=Box_X(1,2); % Since we want to compare with the Boxsize at t=0
            mean_Box_X=mean(Box_X(analyze_start:end,2));
            diff_Box_X=mean(Box_X(:,2)-init_Box_X);
            all_diff_BoxX=[all_diff_BoxX diff_Box_X/mean_Box_X];

            % eval(char(strcat('!echo Box-Y | gmx energy -f',{' '},filename,'.edr -o -b 0 Box_Y.xvg 2>/dev/null')));
            %             !bash run_gmx_energy.sh Box-Y
            Box_Y=import_xvg('Box_Y.xvg');
            init_Box_Y=Box_Y(1,2);
            mean_Box_Y=mean(Box_Y(analyze_start:end,2));
            diff_Box_Y=mean(Box_Y(:,2)-init_Box_Y);
            all_diff_BoxY=[all_diff_BoxY diff_Box_Y/mean_Box_Y];

            % eval(char(strcat('!echo Box-Z | gmx energy -f',{' '},filename,'.edr -o -b 0 Box_Z.xvg 2>/dev/null')));
            %             !bash run_gmx_energy.sh Box-Z
            Box_Z=import_xvg('Box_Z.xvg');
            init_Box_Z=Box_Z(1,2);
            mean_Box_Z=mean(Box_Z(analyze_start:end,2));
            diff_Box_Z=mean(Box_Z(:,2)-init_Box_Z);
            all_diff_BoxZ=[all_diff_BoxZ diff_Box_Z/mean_Box_Z];
            save(strcat('prev_res_',num2str(i),'.mat'));
        catch
            disp('Could not extract the Box_X/Y/Z from the .edr file')
            load(strcat('prev_res_',num2str(i),'.mat'))
            mean_Box_X=1.1*mean_Box_X;
            mean_Box_Y=1.1*mean_Box_Y;
            mean_Box_Z=1.1*mean_Box_Z;
            diff_Box_X=1.1*diff_Box_X;
            diff_Box_Y=1.1*diff_Box_Y;
            diff_Box_Z=1.1*diff_Box_Z;
            % all_diff_BoxX=all_diff_BoxX;
            % all_diff_BoxY=all_diff_BoxY;
            % all_diff_BoxZ=all_diff_BoxZ;
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
            save(strcat('prev_res_',num2str(i),'.mat'));
        catch
            disp('Could not extract the Pressure, Potential, Temperature.. from the .edr file')
            pause(1)
            load(strcat('prev_res_',num2str(i),'.mat'))
            mean_pressure=2*mean_pressure;
            mean_potential=2*mean_potential
            mean_posres=2*mean_posres
            all_mean_posres=[all_mean_posres 2*mean_posres];
            mean_temperature=mean_temperature;
        end

        %% Analyzing the trajectory
        !bash run_gmx_traj.sh

        if isempty(dir(fullfile(pwd,strcat('prev_analysis_*.mat'))))
            try
                copyfile('opt.gro','opt_t0.gro');
                copyfile('temp.top',strcat('../first',num2str(i),'.top'));
                copyfile('temp.top',strcat('first',num2str(i),'.top'));
                copyfile('rmsf.xvg','rmsf_t0.xvg');
                copyfile('rmsdev.xvg','rmsdev_t0.xvg');
                copyfile('rmsdev.pdb','rmsdev_t0.pdb');
            catch
                disp('Could not save the inital simulation output!')
            end
        end
        cd(ParentSimDir);
    end

    res=[];
    all_GII=[];
    indexes=[]; % Matlab struct
    for i=numSimulation %1:numSimulation

        tempdir=strcat(dirtype,num2str(i),'/');

        cd(tempdir);

        try
            disp('Analyzing bond distances/valences!')
            load(strcat('prev_res_',num2str(i),'.mat')); % pause(0.1)
            load('system_info.mat'); % pause(0.1)
            indexes(i).ind=find(strncmpi([atom.resname],resname,2));
            laststep='Run_simulation';
            thermo=readlines('thermo_out.dat');
            iteration=size(thermo,1)-1;

            %     [Y,I] = sort(Bond_index(:,3));
            %     target_Bond_index=Bond_index(I,:);
            %     [Y,I] = sort(Angle_index(:,4));
            %     target_Angle_index=Angle_index(I,:);
            %     [Y,I] = sort(Dihedral_index(:,5));
            %             %     target_Dihedral_index=Dihedral_index(I,:);
            %             target_Bond_index=Bond_index;
            %             target_Angle_index=Angle_index;
            % target_Dihedral_index=Dihedral_index;

            if ~ismember(i,skip_simulation)

                if (max([atom([indexes(i).ind]).x])-min([atom([indexes(i).ind]).x]))>0.9*Box_dim(1) || ...
                        (max([atom([indexes(i).ind]).y])-min([atom([indexes(i).ind]).y]))>0.9*Box_dim(2) || ...
                        (max([atom([indexes(i).ind]).z])-min([atom([indexes(i).ind]).z]))>0.9*Box_dim(3)
                    PBC=1;
                    [atom_ave,Box_dim_ave,Box_dim_std,Cell_ave,Cell_std]=write_ave_gro('conf','median');
                else
                    PBC=0;
                    [atom_ave,Box_dim_ave,Box_dim_std,Cell_ave,Cell_std]=write_ave_gro('conf','ave');
                end

                %% Since we are only looking at a fractional system, since else too many angles, bonds..
                atom_ave=atom_ave([indexes(i).ind]); % Also Box_dim_ave,Cell_ave,Box_dim_std,Cell_std...

                atom_ave = recalc_bond_atom(atom_ave,Box_dim_ave,target_Bond_index,target_Angle_index);%,target_Dihedral_index); % Will output new Bond_index, Angle_index, dihedral_index
                temp_Bond_index=Bond_index;temp_Angle_index=Angle_index;
                GII=0;GII_noH=0;
                %             temp_prop = analyze_atom(element_atom(atom_ave([indexes(i).ind])),Box_dim_ave,0,target_Bond_index);
                temp_prop = analyze_atom(element_atom(atom_ave([indexes(i).ind])),Box_dim_ave,0,Bond_index,[ValenceProperties.oxstate]);
                write_atom_gro(atom_ave,Box_dim_ave,strcat('ave_struct_',num2str(iteration),'.gro'));
                replace_string('Created in Matlab',strcat('GII=',num2str(GII),'; GII_noH=',num2str(GII_noH)),strcat('ave_struct_',num2str(iteration),'.gro'));
            else
                load(strcat('prev_analysis_',num2str(i),'.mat'),'atom_ave','Box_dim_ave','Cell_ave','Cell_std','temp_Bond_index','temp_Angle_index','GII','GII_noH');
                copyfile(strcat('ave_struct_',num2str(iteration-1),'.gro'),strcat('ave_struct_',num2str(iteration),'.gro'));
            end

            
            % % In case of noPBC
            % %     if PBC==0
            % dlmwrite('bond.ndx',{'[ B1B2 ]'},'delimiter','')
            % dlmwrite('bond.ndx',target_Bond_index(:,1:2),'delimiter',' ','-append')
            % %         eval(char(strcat('!echo 0 | gmx distance -f temp.xtc -s opt.tpr -n bond.ndx -rmpbc -pbc yes -len 0.15 -tol 1 -oav -oall -oxyz -oh -oallstat -b 5 2>/dev/null')));
            % %         diststat=import_xvg('diststat.xvg');
            % %         Bond_index=target_Bond_index(:,1:3);
            % %         Bond_index(:,3)=10*diststat(:,2);
            % %         Bond_index(:,4)=10*diststat(:,3);
            % %
            % dlmwrite('angle.ndx',{'[ A2A1A3 ]'},'delimiter','')
            % dlmwrite('angle.ndx',target_Angle_index(:,1:3),'delimiter',' ','-append')
            % %         eval(char(strcat('!echo 0 | gmx angle -f temp.xtc -type angle -b 5 -n angle.ndx -od -ov -all 2>/dev/null')));
            % %         angaver=import_xvg('angaver.xvg');
            % %         angaver(:,1:2)=[];
            % %         Angle_index=target_Angle_index(:,1:3);
            % %         Angle_index(:,4)=mean(angaver)';
            % %     end

            %% To analyze the RDF and CN data
            try
                if numel(rdf_target)>0
                    !bash run_gmx_rdf.sh 'ref' 'sel'
                    rdf=import_xvg('rdf_opt.xvg');rdf=rdf(:,1:2); % Relevant RDF data should be found in column 1:2
                    cn=import_xvg('cn_opt.xvg');cn=cn(:,1:2); % Relevant CN data should be found in column 1:2
                end
            catch
                disp('Did not succeed with the RDF/CN analysis!!!')
                pause(1)
            end


            %% To remove bonds and angles with H's
            ind_H=find(strncmpi([atom_ave.type],{'H'},1));
            [H_row,H_col]=ind2sub(size(target_Bond_index),find(ismember(target_Bond_index,ind_H)));
            target_Bond_index(H_row,:)=[];
            if ~ismember(i,skip_simulation)
                temp_Bond_index(H_row,:)=[];
            end
            [H_row,H_col]=ind2sub(size(target_Angle_index),find(ismember(target_Angle_index,ind_H)));
            target_Angle_index(H_row,:)=[];
            if ~ismember(i,skip_simulation)
                temp_Angle_index(H_row,:)=[];
            end

            Bond_diff=[0];Angle_diff=[0];RDF_diff=[0];CN_diff=[0];CoordNum=0;
            if numel(rdf_target)>0

                rdf((rdf(:,1)<0.125),2)=0;
                rdf_target((rdf_target(:,1)<0.125),2)=0;

                rdf_temp=interp1(rdf(:,1),rdf(:,2),rdf_target(:,1));
                cn_temp=interp1(cn(:,1),cn(:,2),cn_target(:,1));

                rdf_temp=rdf_temp/max(rdf_temp);
                rdf_target(:,2)=rdf_target(:,2)/max(rdf_target(:,2));

                RDF_diff=(rdf_temp-rdf_target(:,2));
                CN_diff=(cn_temp-cn_target(:,2));

                try
                    CoordNum=eval_cn(rdf,cn);
                catch
                    disp('Did not compute CoordNum')
                end

                if PBC==0
                    res_temp=double(RDF_diff); %;CN_diff/100];
                elseif PBC==1
                    res_temp=double(1/10*(abs(diff_Box_X)+abs(diff_Box_Y)+abs(diff_Box_Z)) + RDF_diff); %;CN_diff/100]);
                end
            else

                Bond_diff=double(temp_Bond_index(:,3))-double(target_Bond_index(:,3));
                Angle_diff=double(temp_Angle_index(:,4))-double(target_Angle_index(:,4));

                if PBC==0
                    res_temp=double([Bond_diff;Angle_diff/100]);
                elseif PBC==1
                    % if exist('GII')
                    %     res_temp=double(2*GII + 2*GII_noH + 1/10*(diff_Box_X+diff_Box_Y+2*diff_Box_Z) + [Bond_diff;Angle_diff/100]);
                    % else

                    % res_temp=double(1/10*(abs(diff_Box_X)+abs(diff_Box_Y)+abs(diff_Box_Z)) + [Bond_diff;Angle_diff/100]);
                    res_temp = double(...
                        1/10*( abs(Cell_ave(1)-target_Cell(1)) + abs(Cell_ave(2)-target_Cell(2)) + abs(Cell_ave(3)-target_Cell(3)) )...
                        + 1/10*( abs(Cell_ave(4)-target_Cell(4)) + abs(Cell_ave(5)-target_Cell(5)) + abs(Cell_ave(6)-target_Cell(6)) )...
                        + [Bond_diff;Angle_diff]);
                    % end
                end
            end

            dlmwrite('res_temp.dat',size(res_temp) ,'delimiter',' ','-append')

            %% Make sure we have the same res_temp size due to unknown occasional bug
            if iteration==1
                res_temp_size=size(res_temp,1);
                orig_res_temp=res_temp;
                save('res_temp_size.mat','res_temp_size','orig_res_temp');
            end
            load('res_temp_size.mat')
            res_temp=interp1(1:size(res_temp,1),res_temp,1:res_temp_size)';
            if find(isnan(res_temp))>0
                res_temp=2*orig_res_temp;
            end

            dlmwrite('thermo_out.dat',[iteration norm(res_temp)^2 GII mean(abs(RDF_diff)) mean(abs(CN_diff)) CoordNum mean(abs(Bond_diff)) mean(abs(Angle_diff)) mean_pressure mean_potential mean_posres mean_temperature param mean_Box_X mean_Box_Y mean_Box_Z diff_Box_X/mean_Box_X diff_Box_Y/mean_Box_Y diff_Box_Z/mean_Box_Z Cell_ave Cell_std],'precision','%1.6E','delimiter','\t','-append')
            %         dlmwrite('progress_param.dat',[norm(res)^2 param GII penalty mean_posres/10000 50*diff_Box_X 50*diff_Box_Y 50*diff_Box_Y mean_potential mean_temperature mean_Box_X mean_Box_Y mean_Box_Z],'precision','%1.6E','delimiter','\t','-append')
            if ~ismember(i,skip_simulation)
                save(strcat('prev_analysis_',num2str(i),'.mat'))
            end
            %             save(strcat(_',num2str(i),'.mat'))
            pause(0.1)
            res=[res;res_temp];

            res_size=size(res)

            %             !rm -f conf*.gro
            !rm -f \#*

        catch
            load(strcat('prev_analysis_',num2str(i),'.mat'),'res_temp');
            res_temp=2*res_temp;
            res=[res;res_temp];
            disp('Could not analyze this particular simulation for some reason')
            i
            pause(1)
        end

        prev_param=param;
        save('prev_param.mat','prev_param');

        cd(ParentSimDir);

        save('prev_run.mat')

        all_GII=[all_GII GII];

    end

    ParentSimDirFiles=dir;
    if sum(strcmp({ParentSimDirFiles.name},'progress_param.dat'))==0
        fid = fopen('progress_param.dat','at+');
        fprintf(fid, '%s\n', ('% sim  norm(res)^2  penalty  mean_GII parameters  mean_diff_BoxX   mean_diff_BoxY   mean_diff_BoxZ'));
        fclose(fid);
    end

    penalty=mean(all_mean_posres)/(5000*length(numSimulation));
    %     res=res+penalty;
    %     res=smooth(res,0.02,'rloess')';
    %     res=smooth(rmsd(:,2)+rmsf+penalty',0.02,'rloess');%,11);
    %     res=rmsd(:,2);%,11);

    dlmwrite('progress_param.dat',[iteration norm(res)^2 penalty mean(all_GII) param mean(all_diff_BoxX) mean(all_diff_BoxY) mean(all_diff_BoxZ)],'precision','%1.6E','delimiter','\t','-append')


catch
    disp('The optimization run did not finish for some reason!!!')
    pause(5)
    cd(ParentSimDir);
    load('prev_run.mat')
    penalty=2*penalty;
    res=2*res;
    dlmwrite('progress_param.dat','Crashed','delimiter','\t','-append')
    dlmwrite('res.dat','Crashed','-append')
end

dlmwrite('res.dat',size(res),'delimiter',' ','-append','precision','%12i')
% res=interp1(1:size(res,1),res,1:50000*numSimulation,'nearest',0)';
% dlmwrite('res.dat',size(res),'delimiter',' ','-append','precision','%12i')

save('prev_run.mat')

end
