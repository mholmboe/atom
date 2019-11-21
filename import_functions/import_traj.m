%% import_traj.m
% * This general function imports an structure and a xtc, trr, dcd, xyz or gro trajectory file.
% * For more options, look into each repsective import function, to use
% * stride and other thingies.
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_traj('conf.gro','traj.xtc')
% # atom = import_traj('conf.gro','traj.trr')
% # atom = import_traj('conf.gro','traj.dcd')
% # atom = import_traj('traj.pdb')
% # atom = import_traj('traj.xyz')
% # atom = import_traj('traj.gro')
%
function [atom,traj] = import_traj(filename,varargin)

if nargin==1
    filenameconf=filename;
    filenametraj=filename;
elseif nargin>1
    filenameconf=filename;
    filenametraj=varargin{1};
end

stride=1;
if nargin>2
    stride=varargin{2};
end

if regexp(filenametraj,'.xtc') > 1
    disp('This function imports a structure file and a xtc file')
    disp('and relies on the mxdrfile by Jon Kapla')
    disp('http://kaplajon.github.io/mxdrfile/')
    atom = import_xtc(filenameconf,filenametraj,stride);
    
    %     disp('Alternatively, use the Gro2Mat package with the function:')
    %     disp('import_xtcv2(filenameconf,filenametraj)')
    %     disp('Look into import_traj.m and in the Journal of Computational Chemistry,')
    %     disp('Volume 35, Issue 20')
    %     disp('Found .xtc file, will use the excellent Gro2Mat scripts');
    %     disp('Note that in maxTol, startFrame, endFrame values can be set,')
    %     disp('look into Gro2Mat parseTrj function for help')
    %     disp('Note that lines 55-56 in parseTrj can give unneccesary problems')
    %     import_xtcv2(filenameconf,filenametraj);
    %     delete('./#*'); delete('./temp*');
    %     assignin('caller','trj',trj);
    
elseif regexp(filenametraj,'.trr') > 1
    disp('This function imports a structure file and a trr file')
    disp('and relies on the mxdrfile by Jon Kapla')
    disp('http://kaplajon.github.io/mxdrfile/')
    atom = import_trr(filenameconf,filenametraj,stride);
    
    % disp('Alternatively, use Evans readGmx2Matlab and trr2matlab functions')
    % disp('See http://se.mathworks.com/matlabcentral/fileexchange/33312-convert-gromacs-v-4-5-trajectory-files-into-matlab-matrix')
    % disp('and in the the given import_trrv2.m function')
    % disp('Found .trr file, will use the Evans excellent readGmx2Matlab and trr2matlab functions');
    % disp('Note that in princple both coords, velocities, forces could be imported, and that stride can be used')
    % import_trrv2(filenametraj);
    % assignin('caller','trj',trj);
elseif regexp(filenametraj,'.dcd') > 1
    disp('Found .dcd file, will use the Justin Gullingsrud excellent matdcd script');
    disp('See http://www.ks.uiuc.edu/Development/MDTools/matdcd for some documentation')
    disp(' ')
    disp(' ')
    traj = readdcd(filenametraj,1:nAtoms);
elseif regexp(filenametraj,'.pdb') > 1
    % Note that varargin can be used for importing n frames or using stride
    disp('Found .pdb trajectory');
    atom = import_mc_pdb_traj(filenametraj);
elseif regexp(filenametraj,'.gro') > 1
    % Note that varargin can be used for importing n frames or using stride
    disp('Found .gro trajectory');
    atom = import_gro_traj(filenametraj);
elseif regexp(filenametraj,'.xyz') > 1
    % Note that varargin can be used for importing n frames or using stride
    disp('Found .xyz trajectory');
    atom = import_xyz_traj(filenametraj);
end

if ~exist('atom','var')
    disp('Its always nice to import a structure file first')
    atom = import_atom(filenameconf);
end

assignin('caller','atom',atom)
assignin('caller','traj',traj);
assignin('caller','Box_dim',Box_dim)
assignin('caller','XYZ_labels',XYZ_labels);
assignin('caller','XYZ_data',XYZ_data);


