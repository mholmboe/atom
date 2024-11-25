%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function tries to extract the z component of the velocieties
%% in a Gromacs .trr file, using libxdrfile from gmx. Hacked by MHolmboe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This file is part of mxdrfile.
%%
%% Copyright © 2014 Jon Kapla. All Rights Reserved.
%%
%% Redistribution and use in source and binary forms, with or without
%% modification, are permitted provided that the following conditions are
%% met:
%%
%% 1. Redistributions of source code must retain the above copyright
%%    notice,this list of conditions and the following disclaimer.
%%
%% 2. Redistributions in binary form must reproduce the above copyright
%%    notice, this list of conditions and the following disclaimer in the
%%    documentation and/or other materials provided with the distribution.
%%
%% 3. The name of the author may not be used to endorse or promote products
%%    derived from this software without specific prior written permission.
%%
%% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS "AS IS"
%% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
%% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
%% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
%% DIRECT,INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
%% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
%% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
%% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
%% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%% POSSIBILITY OF SUCH DAMAGE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This scripts imports a gromacs trajectory and rotates by the angles alfa, beta, gamma

%% If you want to add any gromacs commands, remember to set the Gromacs path
%% (here '/usr/local/gromacs-2016.2/bin') for bash within Matlab right:
%% >>PATH = getenv('PATH')
%% >>setenv('PATH', [getenv('PATH'),':','/usr/local/gromacs-2016.2/bin'])
%% then you could try
%% system('gmx editconf -h')

function extract_vz_from_trr(trajname,outtrajname,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loadmxdrfile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3
    stride=1;
else
    stride=varargin{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,rTraj]=inittraj(trajname,'r');
[~,wTraj]=inittraj(outtrajname,'w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pathstr,name,ext] = fileparts(trajname);
if(strcmp(ext,'.trr'))
    frame=0;
    while true % Frameloop
        [rstatus,traj]=read_trr(rTraj);
        if(not(rstatus))
            frame=frame+1;
            %                 disp('Frame'),disp(frame)
        else
            break
        end

        if ~mod(frame,stride) % Test of stride
            %             if ~mod(frame,stride) % Test of stride
            if ~mod(frame,100*stride)
                disp('Frame'),disp(frame)
            end
            %% Do something with the coordinates
            %             V_data(:,1)=traj.v.value(1,:)';
            %             V_data(:,2)=traj.v.value(2,:)';
            V_data(:,3)=traj.v.value(3,:)';

            traj.v.value(1,:)=zeros(1,size(V_data(:,3),1));
            traj.v.value(2,:)=zeros(1,size(V_data(:,3),1));
            traj.v.value(3,:)=V_data(:,3)';

            %% Write newcoords to a new trr file
            wstatus=write_trr(wTraj, traj);

        end % end of test of stride
    end
    [status,rTraj]=closetraj(rTraj);
    [status,wTraj]=closetraj(wTraj);

else
    disp('Could not find .trr file!!!')
end

end

