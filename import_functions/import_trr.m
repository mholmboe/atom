%% import_trr.m
% * This function imports a structure file and a trr file
% and relies on the 'mxdrfile' by Jon Kapla, see below
% * http://kaplajon.github.io/mxdrfile/
% * Hence you need to download his mxdrfile package and add it to your
% Matlab path, see the installation inctructions on Jon's website.
% * This particular function basically just stores the trajectory with 
% coordinates as columns and timesteps as rows, and writes the box 
% dimensions to Box_dim
%
%% Version
% 2.03
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_trr('conf.gro',traj.xtc)
% # atom = import_trr('conf.pdb',traj.xtc,10)
% # atom = import_trr('conf.pdb',traj.xtc,10,10000)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script imports an trr trajectory file
% trr files with libxdrfile from gmx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of mxdrfile.
%
% Copyright Â© 2014 Jon Kapla. All Rights Reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% 1. Redistributions of source code must retain the above copyright
%    notice,this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%
% 3. The name of the author may not be used to endorse or promote products
%    derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
% DIRECT,INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Examples
% * atom = import_trr('conf.gro',traj.trr)
% * atom = import_trr('conf.pdb',traj.trr,10)

function atom = import_trr(filenameconf,filenametrr,varargin)

% If you need to set your local Gromacs PATH (but you could also call PATH2GMX)
% PATH=getenv('PATH');
% setenv('PATH', [getenv('PATH'),':','/usr/local/gromacs-2016.2/bin']);

if regexp(filenameconf,'.gro') > 1
    disp('Found .gro file');
    atom=import_atom_gro(filenameconf);
elseif regexp(filenameconf,'.pdb') > 1
    disp('Found .pdb file');
    filenamepdb=filenameconf;
    atom=import_atom_pdb(filenamepdb);
elseif regexp(filenameconf,'.xyz') > 1
    disp('Found .xyz file');
    disp('Does it carry any Box size info?');
    disp('Its better to use a .ggro or .pdb file...');
    filenamepdb=filenameconf;
    atom=import_atom_pdb(filenamepdb);
else
    disp('Unknown format of the structure file,');
    disp('it should be a .gro or .pdb file');
end

if nargin>2
    stride=varargin{1};
else
    stride=1;
end

if nargin>3
    MaxFrames=varargin{2};
else
    MaxFrames=10000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loadmxdrfile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,rTraj]=inittraj(filenametrr,'r');
nAtoms=rTraj.natoms;
Box_dim=zeros(MaxFrames,9);
traj=zeros(MaxFrames,3*nAtoms);
vtraj=zeros(MaxFrames,3*nAtoms);
ftraj=zeros(MaxFrames,3*nAtoms);
trajstep=zeros(MaxFrames,1);
trajtime=zeros(MaxFrames,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rframe=0;frame=1;
while true % Frameloop
    [rstatus,rtraj]=read_trr(rTraj);
    if(not(rstatus))
        rframe=rframe+stride;
    else
        break
    end
    
    trajstep(frame,1)=rtraj.step;
    trajtime(frame,1)=rtraj.time;

    traj(frame,1:3:end)=rtraj.x.value(1,:)';
    traj(frame,2:3:end)=rtraj.x.value(2,:)';
    traj(frame,3:3:end)=rtraj.x.value(3,:)';
    
    vtraj(frame,1:3:end)=rtraj.v.value(1,:)';
    vtraj(frame,2:3:end)=rtraj.v.value(2,:)';
    vtraj(frame,3:3:end)=rtraj.v.value(3,:)';
    
    ftraj(frame,1:3:end)=rtraj.f.value(1,:)';
    ftraj(frame,2:3:end)=rtraj.f.value(2,:)';
    ftraj(frame,3:3:end)=rtraj.f.value(3,:)';
    
    Box_dim(frame,1)=rtraj.box.value(1,1);
    Box_dim(frame,2)=rtraj.box.value(2,2);
    Box_dim(frame,3)=rtraj.box.value(3,3);
    
    Box_dim(frame,4)=rtraj.box.value(2,1);
    Box_dim(frame,5)=rtraj.box.value(3,1);
    
    Box_dim(frame,6)=rtraj.box.value(1,2); % #6 xy
    
    Box_dim(frame,7)=rtraj.box.value(3,2);
         
    Box_dim(frame,8)=rtraj.box.value(1,3); % #8 xz
    Box_dim(frame,9)=rtraj.box.value(2,3); % #9 yz
    
    frame=frame+1;
    
    if mod(rframe,100)==0
        frame-1
    end
    
end
[status,rTraj]=closetraj(rTraj);

traj(frame:end,:)=[];
Box_dim(frame:end,:)=[];
if sum(sum(Box_dim(:,4:9))) == 0
    Box_dim(:,4:end)=[];
end

Box_dim=10*Box_dim; % nm to Å

if sum(sum(traj)) ~= 0
    traj=10*traj; % nm to Å
    assignin('caller','traj',traj)
end

if sum(sum(vtraj)) ~= 0
    assignin('caller','vtraj',vtraj)
end

if sum(sum(ftraj)) ~= 0
    assignin('caller','ftraj',ftraj)
end

disp('Note the units, length in Å, v and f uses Gromacs nm units')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delete('./#*'); delete('./temp*');
assignin('caller','atom',atom)
assignin('caller','Box_dim',Box_dim)
assignin('caller','trajstep',trajstep)
assignin('caller','trajtime',trajtime)
assignin('caller','XYZ_labels',XYZ_labels)
assignin('caller','XYZ_data',XYZ_data)
