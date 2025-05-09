%% import_xyz_traj.m
% * This function imports a .xyz trajectory
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # [atom,traj] = import_xyz_traj('traj.xyz')
% # [atom,traj] = import_xyz_traj('traj.xyz',1000,10)

function [atom,traj] = import_xyz_traj(filename,varargin)
%%

fileID = fopen(filename,'r');
Line1 = {fgets(fileID)};
Line2 = {fgets(fileID)};
nAtoms=str2double(Line1);
Box_cell=strsplit(char(Line2));
LineEnd = textscan(fileID, '%s',inf);
nFrames=sum(ismember(LineEnd{1,1},{num2str(nAtoms)}))+1;

fclose(fileID);

idx = strfind(char(Box_cell(1)),'#');
Box_dim=[];
if numel(idx) >0
    for i=1:length(Box_cell)
        [num, status] = str2num(char(Box_cell(i)));
        j=1;
        if status==1
            Box_dim=[Box_dim num];
            j=j+1;
        end
    end
end

Boxinfo='maybe';
if nargin>1
    Frames=[1:varargin{1}];
    if nargin>2
        Frames=[1:varargin{2}:max(Frames)];
    end
else
    if numel(Box_dim)==0
        Frames=1:nFrames;
        Boxinfo='no';
    else
        Frames=1:nFrames;
    end
end
Frames;
traj=zeros(length(Frames),nAtoms*3);
Box_dim=zeros(length(Frames),3);
Comments={};
fileID = fopen(filename,'r');

if strcmp(Boxinfo,'no')
    for t=1:length(Frames)
        
        frewind(fileID)
        startRow = 1+nAtoms*(Frames(t)-1)+2*(Frames(t)-1);

        Line2 = textscan(fileID, '%[^\n\r]', 1,'HeaderLines',startRow);
        
        Comments(t,1)=Line2;

        
        formatSpec = '%s%f%f%f%[^\n\r]';
        dataArray = textscan(fileID, formatSpec, nAtoms, 'Delimiter',{'\t',' '}, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines',0, 'ReturnOnError', false);
        
        traj(t,1:3:end)=[dataArray{:,2}];
        traj(t,2:3:end)=[dataArray{:,3}];
        traj(t,3:3:end)=[dataArray{:,4}];
        
        % Extra stuff - will print a .xyz for every frame. Note you could set the Box_dim
        XYZ_labels = dataArray{:,1};
        XYZ_data = [dataArray{:,2} dataArray{:,3} dataArray{:,4}];
        % write_xyz(XYZ_labels,XYZ_data,[],strcat('out_',num2str(t),'.xyz'),Comments(t,1));
        
        if t==1
            XYZ_labels = dataArray{:,1};
            XYZ_data = [dataArray{:,2} dataArray{:,3} dataArray{:,4}];
        end
       
        Box_dim(t,:) = [(max([dataArray{:,2}])-min([dataArray{:,2}]))*1.001    (max([dataArray{:,3}])-min([dataArray{:,3}]))*1.001  (max([dataArray{:,4}])-min([dataArray{:,4}]))*1.001];
        
        if mod(t,10)==0
            t
        end
    end
    assignin('caller','Comments',Comments);
    disp('Guessing the box dimensions to be .1% larger than the max-min distances')
else
    for t=1:length(Frames)
        
        frewind(fileID)
                startRow = 1+nAtoms*(Frames(t)-1)+2*(Frames(t)-1)
        
%         Line1 = {fgets(fileID)}
        Line2 = textscan(fileID, '%s%[^\n\r]', 1,'HeaderLines',startRow);
        
%         Line1 = {fgets(fileID)};
%         Line2 = {fgets(fileID)};
%         nAtoms=str2double(Line1);
        Box_cell=strsplit(char(Line2));
        
        idx = strfind(char(Box_cell(1)),'#');
        Box_dim=[];
        if numel(idx) >0
            for i=1:length(Box_cell)
                [num, status] = str2num(char(Box_cell(i)));
                j=1;
                if status==1
                    Box_dim=[Box_dim num];
                    j=j+1;
                end
            end
        end
        
        formatSpec = '%s%f%f%f%[^\n\r]';
        dataArray = textscan(fileID, formatSpec, nAtoms, 'Delimiter',{'\t',' '}, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines',0, 'ReturnOnError', false);
        
        traj(t,1:3:end)=[dataArray{:,2}];
        traj(t,2:3:end)=[dataArray{:,3}];
        traj(t,3:3:end)=[dataArray{:,4}];
        
%         % Extra stuff - will print a .xyf for every frame. Note you could set the Box_dim
%         XYZ_labels = dataArray{:,1};
%         XYZ_data = [dataArray{:,2} dataArray{:,3} dataArray{:,4}];
%         write_xyz(XYZ_labels,XYZ_data,[15 15 15],strcat('out_',num2str(t),'.xyz'));

        
        if t==1
            XYZ_labels = dataArray{:,1};
            XYZ_data = [dataArray{:,2} dataArray{:,3} dataArray{:,4}];
        end
        
        if mod(t,10)==0
            t
        end
        
    end
end
fclose(fileID);

if sum(traj(end,:)==0)
    traj(end,:)=[];
end

for i=1:nAtoms
    atom(i).resname={'MOL'};
    atom(i).molid=1;
    atom(i).type    = XYZ_labels(i,1);
    atom(i).fftype  = XYZ_labels(i,1);
    atom(i).charge      = 0;
    atom(i).index=mod(i,100000);
    atom(i).neigh.type  = {};
    atom(i).neigh.index  = zeros(6,1);
    atom(i).neigh.dist  = zeros(6,1);
    atom(i).bond.type  = zeros(6,1);
    atom(i).bond.index  = zeros(6,1);
    atom(i).angle.type  = zeros(6,1);
    atom(i).angle.index  = zeros(6,1);
    atom(i).x = XYZ_data(i,1);
    atom(i).y = XYZ_data(i,2);
    atom(i).z = XYZ_data(i,3);
    atom(i).vx=0;
    atom(i).vy=0;
    atom(i).vz=0;
end

atom = resname_atom(atom);

% assignin('caller','atom',atom);
% assignin('caller','traj',traj);
assignin('caller','nAtoms',nAtoms);
assignin('caller','Box_dim',Box_dim);
assignin('caller','XYZ_labels',XYZ_labels);
assignin('caller','XYZ_data',XYZ_data);


if Frames(end) == size(traj,1)
    sprintf('.xyz traj file imported %d frames', Frames(end))
else
    sprintf('Check the number of frames!')
end


