%% import_pdb_traj.m
% * This function imports a .pdb trajectory
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_pdb_traj(filename)
% # atom = import_pdb_traj(filename,maxFrames)
% # atom = import_pdb_traj(filename,maxFrames,stride)
%
function atom = import_pdb_traj(filename,varargin)

inputfile = fopen(filename, 'r');
C = textscan(inputfile, '%s', 'Delimiter', '\n');
fclose(inputfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% System dimensions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    %% CRYST1 section
    CRYST1_rows = strfind(C{1}, 'CRYST1');
    CRYST1_rows = find(~cellfun('isempty', CRYST1_rows));
    CRYST1_data=cell(numel(CRYST1_rows),7);
    for i=1:numel(CRYST1_rows)
        CRYST1_data(i,:)=strsplit(char(C{1,1}(CRYST1_rows(i))));
    end
    CRYST1_data=str2double(CRYST1_data(:,2:end));
    
    a=CRYST1_data(:,1);
    b=CRYST1_data(:,2);
    c=CRYST1_data(:,3);
    alfa=CRYST1_data(:,4);
    beta=CRYST1_data(:,5);
    gamma=CRYST1_data(:,6);
    lx = a;
    xy = b.*cos(deg2rad(gamma));xy=round(xy,4);
    ly = (b.^2-xy.^2).^.5;
    xz = c.*cos(deg2rad(beta));xz=round(xz,4);
    yz = (b.*c.*cos(deg2rad(alfa))-xy.*xz)./ly;yz=round(yz,4);
    lz = (c.^2-xz.^2-yz.^2).^0.5;
    
    Box_dim=[lx ly lz zeros(numel(lx),1) zeros(numel(lx),1) xy zeros(numel(lx),1) xz yz];
catch
    disp('Could not fetch the CRYST1 Box_dim data')
    CRYST1_rows=[];
    CRYST1_data=[];
    Box_dim=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    %% MODEL section
    MODEL_rows = strfind(C{1}, 'MODEL');
    MODEL_rows = find(~cellfun('isempty', MODEL_rows));
    MODEL_data=cell(numel(MODEL_rows),2);
    for i=1:numel(MODEL_rows)
        MODEL_data(i,:)=strsplit(char(C{1,1}(MODEL_rows(i))));
    end
    MODEL_data=str2double(MODEL_data(:,2));
    
catch
    disp('Could not fetch the MODEL lines')
    MODEL_rows=[];
    MODEL_data=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    %% TER section
    TER_rows = strfind(C{1}, 'TER');
    TER_rows = find(~cellfun('isempty', TER_rows));
    TER_data=cell(numel(TER_rows),2);
    for i=1:numel(TER_rows)
        TER_data(i,:)=strsplit(char(C{1,1}(TER_rows(i))));
    end
    TER_data=str2double(TER_data(:,2));
catch
    disp('Could not fetch the TER lines')
    TER_rows=[];
    TER_data=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    %% ENDMDL section
    ENDMDL_rows = strfind(C{1}, 'ENDMDL');
    ENDMDL_rows = find(~cellfun('isempty', ENDMDL_rows));
    ENDMDL_data=cell(numel(ENDMDL_rows),2);
    for i=1:numel(ENDMDL_rows)
        ENDMDL_data(i,:)=strsplit(char(C{1,1}(ENDMDL_rows(i))));
    end
    ENDMDL_data=str2double(ENDMDL_data(:,2));
catch
    disp('Could not fetch the ENDMDL lines')
    ENDMDL_rows=[];
    ENDMDL_data=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA section
%% Init
if numel(CRYST1_rows)>0
    if numel(MODEL_rows)>0
        if CRYST1_rows(1)<MODEL_rows(1)
            DATA_init_rows=1+MODEL_rows;
        else
            DATA_init_rows=1+CRYST1_rows;
        end
    end
end
%% End
if numel(ENDMDL_rows)>0
    if numel(TER_rows)>0
        if ENDMDL_rows(1)<TER_rows(1)
            DATA_end_rows=-1+TER_rows;
        else
            DATA_end_rows=-1+ENDMDL_rows;
        end
    else
        DATA_end_rows=-1+ENDMDL_rows;
    end
elseif numel(TER_rows)>0
    DATA_end_rows=-1+TER_rows;
end

if numel(unique(DATA_end_rows-DATA_init_rows+1))>1
    [nAtoms,nAtoms_max_ind]=max(DATA_end_rows-DATA_init_rows+1);
else
    nAtoms=DATA_end_rows(1)-DATA_init_rows(1)+1;
    nAtoms_max_ind=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin>1
    Frames=1:varargin{1};
    if nargin>2
        Frames=1:varargin{2}:max(Frames);
    end
else
    Frames=1:MODEL_data(end);
end

j=1;atom=[];
for i=DATA_init_rows(nAtoms_max_ind):DATA_end_rows(nAtoms_max_ind)
    line = char(C{1,1}(i));
    atom(j).molid = str2double(line(23:26));
    atom(j).resname = {strtrim(line(18:20))};
    %         atom(j).type = {strtrim(line(13:16))}; % Changed to 17 for better
    %         compatiblity
    %         atom(j).fftype = {strtrim(line(13:16))}; % Changed to 17 for better
    %         compatiblity
    atom(j).type = {strtrim(line(13:17))};
    atom(j).fftype = {strtrim(line(13:17))};
    atom(j).index = mod(j,100000); %str2double(line(7:11));
    atom(j).neigh.type = {};
    atom(j).neigh.index = [0;0;0;0;0;0];
    atom(j).neigh.dist = [0;0;0;0;0;0];
    atom(j).bond.type = [0;0;0;0;0;0];
    atom(j).bond.index = [0;0;0;0;0;0];
    atom(j).angle.type = [0;0;0;0;0;0];
    atom(j).angle.index = [0;0;0;0;0;0];
    atom(j).x = str2double(line(31:38));
    atom(j).y = str2double(line(39:46));
    atom(j).z = str2double(line(47:54));
    atom(j).vx = NaN;
    atom(j).vy = NaN;
    atom(j).vz = NaN;
    
    try
        occupancy(j,1)=str2double(line(55:60));
        tempfactor(j,1)=str2double(line(61:66));
    catch
    end
    j = j + 1;
end

if isnan([atom(1).molid])
    [atom.molid]=deal(1);
end

traj=zeros(size(Frames,2),max(nAtoms)*3);
for i=1:length(Frames)
    data=char(C{1,1}(DATA_init_rows(i):DATA_end_rows(i)));
    temp_nAtoms=size(data,1);
    traj(i,1:3:3*temp_nAtoms)=str2num(data(:,31:38));
    traj(i,2:3:3*temp_nAtoms)=str2num(data(:,39:46));
    traj(i,3:3:3*temp_nAtoms)=str2num(data(:,47:54));
    if mod(i,10)==0
        i
    end
end

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
XYZ_labels=[atom.type]';

assignin('caller','atom',atom);
assignin('caller','traj',traj);
assignin('caller','nAtoms',max(nAtoms));
assignin('caller','Box_dim',Box_dim);
assignin('caller','XYZ_labels',XYZ_labels);
assignin('caller','XYZ_data',XYZ_data);

if Frames(end) == size(traj,1)
    sprintf('.pdb traj file imported %d frames', Frames(end))
else
    sprintf('Number of frames is wrong somehow')
end


