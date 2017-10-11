%% import_atom_multiple_gro.m
% * This function import multiple .gro files. Does it work? Untested

function frame = import_atom_multiple_gro(filename,varargin)
%%

C = textread(filename, '%s','delimiter', '\n');
Title=C{1}; nAtoms=str2num(cell2mat(C(2)));

if nargin>1;
    nFrames=[1:varargin{1}];%floor(size(C,1)/(nAtoms+3))];
    if nargin>2;
        nFrames=[1:varargin{2}:max(nFrames)];%floor(size(C,1)/(nAtoms+3))];
    end
else
    nFrames=1:size(C,1)/(nAtoms+3);
end

frame=[];
traj=zeros(length(nFrames),nAtoms*3);
for t=1:length(nFrames)
    
    startRow = 3+nAtoms*(nFrames(t)-1)+3*(nFrames(t)-1);
    endRow = 2+(nFrames(t))*nAtoms+3*(nFrames(t)-1);
    Box_dim=str2double(strsplit(C{endRow+1}))*10;
    
    fileID = fopen(filename,'r');
    %% Read columns of data as strings:
    formatSpec = '%5s%5s%5s%5.0f%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines', startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % atom=struct('molid',strtrim(dataArray{:,1}));
    %% Have to use the deal() function for the others...
    
    nAtoms=size(dataArray{:,5}(:),1);
    MolID = str2double(dataArray{:,1}); % Converts to double
    Resname = strtrim(dataArray{:,2});
    atomtype = strtrim(dataArray{:,3});
    ind=num2cell(mod(1:nAtoms,100000));
    
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
    [atom.resname]=deal(Resname{:});
    [atom.type]=deal(atomtype{:});
    [atom.fftype]=deal(atomtype{:});
    [atom.index]=deal(ind{:});
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    X_coord=num2cell(dataArray{:,5}*10);
    Y_coord=num2cell(dataArray{:,6}*10);
    Z_coord=num2cell(dataArray{:,7}*10);
    X_velo=num2cell(dataArray{:,8}*10);
    Y_velo=num2cell(dataArray{:,9}*10);
    Z_velo=num2cell(dataArray{:,10}*10);
    
    [atom.x]=deal(X_coord{:});
    [atom.y]=deal(Y_coord{:});
    [atom.z]=deal(Z_coord{:});
    [atom.vx]=deal(X_velo{:});
    [atom.vy]=deal(Y_velo{:});
    [atom.vz]=deal(Z_velo{:});
    
%     startRow = 3+nAtoms*(nFrames(t)-1)+3*(nFrames(t)-1);
%     endRow = 2+(nFrames(t))*nAtoms+3*(nFrames(t)-1);
%     
%     data=rawtraj(startRow:endRow,:);
% %     Box_dim=str2double(strsplit(C{endRow+1}))*10;
%     
    traj(t,1:3:end)=[dataArray{:,5}]*10;
    traj(t,2:3:end)=[dataArray{:,6}]*10;
    traj(t,3:3:end)=[dataArray{:,7}]*10;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
    XYZ_labels=[atom.type]';
    
    frame(t).atom=atom;
    frame(t).Box_dim=Box_dim;
    
end

assignin('caller','XYZ_labels',XYZ_labels)
assignin('caller','XYZ_data',XYZ_data)
assignin('caller','frame',frame)
assignin('caller','traj',traj)
assignin('caller','nAtoms',nAtoms)
assignin('caller','Box_dim',Box_dim)
assignin('caller','MolID',MolID)

if nFrames(end) == size(frame,2)
    sprintf('.gro file imported %d frames', nFrames(end))
else
    sprintf('Number of frames is wrong somehow')
end


