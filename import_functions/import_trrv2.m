%% import_trr.m
% * This function imports an trr file using Evan's scripts, see below
% * This function uses trr2matlab + readGmx2Matlab to read gromacs .trr trajectory files
% * into the variable 'traj'. Note that trr2matlab can output both coordinates, velocities and forces
% * if they exist in the .trr file
% * Please report problems/bugs to michael.holmboe@umu.se
% * Can varargin be used for stride?

%% Examples
% * traj = import_trr(traj.trr)
% * traj = import_trr(traj.trr,varargin)


function traj = import_trrv2(filenametrr,varargin)
%%

disp('Found .trr file, will use the Evans excellent readGmx2Matlab and trr2matlab functions');
disp('See http://se.mathworks.com/matlabcentral/fileexchange/33312-convert-gromacs-v-4-5-trajectory-files-into-matlab-matrix')
disp(' ')
disp(' ')
disp('Note that in princple both coords, velocities, forces could be imported, and that stride can be used')
disp('Look into Evans readGmx2Matlab and trr2matlab functions for options')


writeX=0;writeV=0;writeF=0;
% First Evan's script does this... see his comments below for options
if nargin==1
    trr2matlab(filenametrr,'x');
else
    evalstr=eval("strcat('trr2matlab(filenametrr')");
    if sum(strcmp(varargin, 'x')) > 0
        x='''x''';
        evalstr=eval("strcat(evalstr,'comma',x)");
        writeX = true;
    end
    if sum(strcmp(varargin, 'v')) > 0
        v='''v''';
        evalstr=eval("strcat(evalstr,'comma',v)");
        writeV = true;
    end
    if sum(strcmp(varargin, 'f')) > 0
        f='''f''';
        evalstr=eval("strcat(evalstr,'comma',f)");
        writeF = true;
    end
    evalstr=strcat(evalstr,')');
    evalstr=strrep(evalstr,'comma',',');
    eval(evalstr);
end

% For trr2matlab.m, published on Matlab fileexchange
% by Evan Arthur, University of Michigan, October 2011
%
% trr2matlab.m
% Inputs:
%     - path to trr file (required, must be first input)
%     - path to output file (optional)
%         if none given, a default name is chosen (such as 'xdata.binary')
%     - 'x' (optional)
%         outputs xyz atomic coordinates
%     - 'v' (optional)
%         outputs xyz of instantaneous velocities
%     - 'f' (optional)
%         outputs xyz of atomic forces
%     - 'single' or 'double' (optional)
%         selects whether to read file as single or double precision
%         This program automatically detects this. Use this option as
%         an override.
%     - [integer] (optional)
%         output frequency on display window of the frame currently being
%         read. 1 = every frame's statistics are output. 1000 = every 1000
%         frames statistics are output

% Outputs:
%     - xyz data
%         output either by default or if 'x' option is given
%         default name is 'xdata.binary'
%     - velocity data
%         output either by default or if 'v' option is given
%         default name is 'vdata.binary'
%     - force data
%         output either by default or if 'f' option is given
%         default name is 'fdata.binary'

% Example inputs and outputs:
%     trr2matlab ('traj.trr')
%           outputs all atomic coordinates, velocities, and forces as files
%           'xdata.binary', 'vdata.binary', and 'fdata.binary'
%     trr2matlab ('traj.trr', 'x', 'f')
%           outputs all atomic coordinates and forces as files
%           'xdata.binary' and 'fdata.binary' (velocity data is not output)
%     trr2matlab ('traj.trr', 'x')
%           outputs only atomic coordinates as file 'xdata.binary'
%           (velocity and force data are not output)
%     trr2matlab ('traj.trr', 'f', 'proteinA')
%           outputs only atomic forces as file 'proteinA_xdata.binary'
%           (velocity and coordinates data are not output)
%     trr2matlab ('traj.trr', 'x', 'single')
%           outputs only atomic forces as file 'xdata.binary'
%           (velocity and force data are not output)
%           forces single precision mode
%     trr2matlab ('traj.trr', 'x', 10)
%           outputs only atomic forces as file 'xdata.binary'
%           (velocity and force data are not output)
%           outputs statistics of information in window every 10 frames

% notes on Single/Double Precision:
%     This program detects the precision of the trr file automatically. If
%     the detection fails, write in the input 'single' or 'double'. The
%     program outputs garbage or fails spectacularly when this is not done
%     properly.
%
%     Since single precision exceeds the margin of error for most analyses,
%     this program only outputs data as single-precision numbers. If I get
%     requests for the excrutiatingly accurate double-precision output, I
%     will put in an option for it.

if (writeX+writeV+writeF)==0
    [trj] = readGmx2Matlab('xdata.binary');
    
    nAtoms=trj.num_atoms;
    nFrames=trj.num_frames;
    timestep=trj.time_step;
    
    xdata=(reshape(trj.trajectory(:,1,:)*10,nAtoms,nFrames)).';
    ydata=(reshape(trj.trajectory(:,2,:)*10,nAtoms,nFrames)).';
    zdata=(reshape(trj.trajectory(:,3,:)*10,nAtoms,nFrames)).';
    traj=zeros(size(xdata,1),size(xdata,2)*3);
    traj(:,1:3:end)=xdata;
    traj(:,2:3:end)=ydata;
    traj(:,3:3:end)=zdata;
    delete('./#*'); delete('./temp*');
    assignin('caller','traj',traj)
    
    delete('./#*'); delete('./temp*');
end

if writeF>0
    [trj] = readGmx2Matlab('fdata.binary');
    nAtoms=trj.num_atoms;
    nFrames=trj.num_frames;
    timestep=trj.time_step;
    fxdata=(reshape(trj.trajectory(:,1,:)*10,nAtoms,nFrames)).';
    fydata=(reshape(trj.trajectory(:,2,:)*10,nAtoms,nFrames)).';
    fzdata=(reshape(trj.trajectory(:,3,:)*10,nAtoms,nFrames)).';
    ftraj=zeros(size(fxdata,1),size(fxdata,2)*3);
    ftraj(:,1:3:end)=fxdata;
    ftraj(:,2:3:end)=fydata;
    ftraj(:,3:3:end)=fzdata;
    delete('./#*'); delete('./temp*');
    assignin('caller','ftraj',ftraj);
end

if writeV>0
    [trj] = readGmx2Matlab('vdata.binary');
    
    nAtoms=trj.num_atoms;
    nFrames=trj.num_frames;
    timestep=trj.time_step;
    
    vxdata=(reshape(trj.trajectory(:,1,:)*10,nAtoms,nFrames)).';
    vydata=(reshape(trj.trajectory(:,2,:)*10,nAtoms,nFrames)).';
    vzdata=(reshape(trj.trajectory(:,3,:)*10,nAtoms,nFrames)).';
    vtraj=zeros(size(vxdata,1),size(vxdata,2)*3);
    vtraj(:,1:3:end)=vxdata;
    vtraj(:,2:3:end)=vydata;
    vtraj(:,3:3:end)=vzdata;
    delete('./#*'); delete('./temp*');
    assignin('caller','vtraj',vtraj);
end

if writeX>0
    [trj] = readGmx2Matlab('xdata.binary');
    
    nAtoms=trj.num_atoms;
    nFrames=trj.num_frames;
    timestep=trj.time_step;
    
    xdata=(reshape(trj.trajectory(:,1,:)*10,nAtoms,nFrames)).';
    ydata=(reshape(trj.trajectory(:,2,:)*10,nAtoms,nFrames)).';
    zdata=(reshape(trj.trajectory(:,3,:)*10,nAtoms,nFrames)).';
    traj=zeros(size(xdata,1),size(xdata,2)*3);
    traj(:,1:3:end)=xdata;
    traj(:,2:3:end)=ydata;
    traj(:,3:3:end)=zdata;
    
    delete('./#*'); delete('./temp*');
    assignin('caller','traj',traj);
end

clear trj
trj.nAtoms=nAtoms;
trj.nFrames=nFrames;
trj.timestep=timestep;
assignin('caller','trj',trj);

% Assign the Box_dim matrix to the calling workspace
assignin('caller','Box_dim',Box_dim);

