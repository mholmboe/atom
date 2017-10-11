%% import_traj.m
% * This function imports an structure and a xtc, trr, dcd, xyz or gro trajectory file.
% * Does it work?
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * atom = import_traj('conf.gro',traj.xtc)
% * atom = import_traj('conf.gro',traj.trr)
% * atom = import_traj('conf.gro',traj.dcd)
% * atom = import_traj('traj.xyz')
% * atom = import_traj('traj.gro')


function atom = import_traj(filenameconf,filenametraj,varargin)
%

% Its always nice to import a structure file first
atom = import_atom(filenameconf);

if regexp(filenametraj,'.xtc') > 1
    disp('Found .xtc file, will use the excellent Gro2Mat scripts');
    disp('See Journal of Computational Chemistry,Volume 35, Issue 20')
    
    disp('Note that in maxTol, startFrame, endFrame values can be set,')
    disp('look into Gro2Mat parseTrj function for help')
    disp('Note that lines 55-56 in parseTrj can give unneccesary problems')
    import_xtc(filenameconf,filenametraj);
    delete('./#*'); delete('./temp*');
    assignin('caller','trj',trj);
%     assignin('caller','traj',traj);
elseif regexp(filenametraj,'.trr') > 1
    disp('Found .trr file, will use the Evans excellent readGmx2Matlab and trr2matlab functions');
    disp('See http://se.mathworks.com/matlabcentral/fileexchange/33312-convert-gromacs-v-4-5-trajectory-files-into-matlab-matrix')
    
    disp('Note that in princple both coords, velocities, forces could be imported, and that stride can be used')
    disp('Look into Evans readGmx2Matlab and trr2matlab functions for options')
    import_trr(filenametraj);
    assignin('caller','trj',trj);
elseif regexp(filenametraj,'.dcd') > 1
    disp('Found .dcd file, will use the Justin Gullingsrud excellent matdcd script');
    disp('See http://www.ks.uiuc.edu/Development/MDTools/matdcd for some documentation')
    disp(' ')
    disp(' ')
    traj = readdcd(filenametraj,1:nAtoms);
elseif regexp(filenametraj,'.gro') > 1
    % Note that varagin can be used for importing n frames or using stride
    disp('Found .gro trajectory');
    atom = import_gro_traj(filenametraj);
elseif regexp(filenametraj,'.xyz') > 1
    % Note that varagin can be used for importing n frames or using stride
    disp('Found .xyz trajectory');
    atom = import_xyz_traj(filenametraj);
end

assignin('caller','atom',atom)
assignin('caller','traj',traj);
assignin('caller','Box_dim',Box_dim)
assignin('caller','XYZ_labels',XYZ_labels);
assignin('caller','XYZ_data',XYZ_data);


