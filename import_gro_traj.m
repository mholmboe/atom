%% import_gro_traj.m
% * This function imports a .gro trajectory
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * atom = import_gro_traj(filename)
% * atom = import_gro_traj(filename,maxFrames)
% * atom = import_gro_traj(filename,maxFrames,stride)

function atom = import_gro_traj(filename,varargin)


% Get the number of atoms and frames
fileID = fopen(filename,'r');
Line1 = {fgets(fileID)};
Line2 = {fgets(fileID)};
Title=strsplit(char(Line2));
nAtoms=str2double(Line2);
LineEnd = textscan(fileID, '%s',inf,'delimiter', '\n');
fclose(fileID);

disp('Try this if nFrames gives you problems')
disp('nFrames=(size(LineEnd{1,1},1)-1)/(nAtoms+2)-1;')
nFrames=(size(LineEnd{1,1},1)-1)/(nAtoms+2);

% Set the number of frames to a lower value, or use stride
if nargin>1;
    Frames=[1:varargin{1}];
    if nargin>2;
        Frames=[1:varargin{2}:max(Frames)];
    end
else
    Frames=1:nFrames;
end

%
Box_dim=[];atom=[];
traj=zeros(length(Frames),nAtoms*3);
velo=zeros(length(Frames),nAtoms*3);
fileID = fopen(filename,'r');
for t=1:length(Frames)
    
    frewind(fileID)
    startRow = 1+nAtoms*(Frames(t)-1)+3*(Frames(t)-1);
    
    Line2 = textscan(fileID, '%s%[^\n\r]', 1,'HeaderLines',startRow);
    
    formatSpec = '%5s%5s%5s%5.0f%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, nAtoms, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines',1, 'ReturnOnError', false);
    
    assignin('caller','dataArray',dataArray)
    traj(t,1:3:end)=[dataArray{:,5}]*10;
    traj(t,2:3:end)=[dataArray{:,6}]*10;
    traj(t,3:3:end)=[dataArray{:,7}]*10;
    velo(t,1:3:end)=[dataArray{:,8}]*10;
    velo(t,2:3:end)=[dataArray{:,9}]*10;
    velo(t,3:3:end)=[dataArray{:,10}]*10;
    
    Box_string = textscan(fileID, '%s',1,'delimiter', '\n','HeaderLines', 1);
    Box_dim(t,:)=str2double(strsplit(char(Box_string{1,1})))*10;
    
    if t==1;
        
        MolID = str2double(dataArray{:,1}); % Converts to double
        Resname = strtrim(dataArray{:,2});
        atomtype = strtrim(dataArray{:,3});
        ind=num2cell(mod(1:nAtoms,100000));
        
        XYZ_labels = dataArray{:,3};
        XYZ_data = [dataArray{:,5} dataArray{:,6} dataArray{:,7}];
        
        nmol=1;first_in=[1];last_in=[];
        for i=1:nAtoms;
            if i > 1 && MolID(i) ~= MolID(i-1)
                nmol=nmol+1;
                atom(i).molid=nmol;
                first_in(atom(i).molid,1)=i;last_in(atom(i).molid-1,1)=i-1;
            elseif i > 1
                atom(i).molid=atom(i-1).molid;
            elseif i == 1
                atom(i).molid=1;
            end
            %atom(i).index=mod(i,100000);
        end
        %[atom.molid]=deal(MolID{:});
        %         [atom.resname]=deal(Resname{:});
        %         [atom.type]=deal(atomtype{:});
        %         [atom.fftype]=deal(atomtype{:});
        %         [atom.index]=deal(ind{:});
        
        nmol=1;first_in=[1];last_in=[];
        for i=1:nAtoms;
            if i > 1 && MolID(i) ~= MolID(i-1)
                nmol=nmol+1;
                atom(i).molid=nmol;
                first_in(atom(i).molid,1)=i;last_in(atom(i).molid-1,1)=i-1;
            elseif i > 1
                atom(i).molid=atom(i-1).molid;
            elseif i == 1
                atom(i).molid=1;
            end
            atom(i).resname=Resname(i);
            atom(i).type=atomtype(i);
            atom(i).fftype=atomtype(i);
            atom(i).index=mod(i,100000);
            atom(i).neigh.type  = {};
            atom(i).neigh.index  = zeros(6,1);
            atom(i).neigh.dist  = zeros(6,1);
            atom(i).bond.type  = zeros(6,1);
            atom(i).bond.index  = zeros(6,1);
            atom(i).angle.type  = zeros(6,1);
            atom(i).angle.index  = zeros(6,1);
            atom(i).x=traj(1,1+3*(i-1));
            atom(i).y=traj(1,2+3*(i-1));
            atom(i).z=traj(1,3+3*(i-1));
            atom(i).vx=velo(1,1+3*(i-1));
            atom(i).vy=velo(1,2+3*(i-1));
            atom(i).vz=velo(1,3+3*(i-1));
        end
        
    end
    if mod(i,10)==0;
        i
    end
end

%     atom(i).fx=X_coord(i)/Box_dim(1);
%     atom(i).fy=Y_coord(i)/Box_dim(2);
%     atom(i).fz=Z_coord(i)/Box_dim(3);

fclose(fileID);

atom = resname_atom(atom);

assignin('caller','MolID',MolID);
assignin('caller','atom',atom);
assignin('caller','traj',traj);
assignin('caller','velo',velo);
assignin('caller','nAtoms',nAtoms);
assignin('caller','Box_dim',Box_dim);
assignin('caller','XYZ_data',XYZ_data);
assignin('caller','XYZ_labels',XYZ_labels);

disp('.gro file imported')

if Frames(end) == size(traj,1)
    sprintf(' %d configurations/frames found in imported .gro file', Frames(end))
else
    sprintf('Number of frames is wrong somehow')
end

