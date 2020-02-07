%% create_grid_atom.m
% * This old function puts ions on a grid plane and adds it to an atom struct. 
% * Use create_atom() instead.
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * atom = grid2atom(atom_label,12,[20 15 2],'xy',[5 5 0])

function atom = create_grid_atom(atom_label,nM,limits,dim,varargin)
 
% if length(limits) == 3;

if strcmp(dim,'xy')
    lx=limits(1);
    ly=limits(2);
    lz=limits(3);
elseif strcmp(dim,'xz')
    lx=limits(1);
    ly=limits(3);
    lz=limits(2);
elseif strcmp(dim,'yz')
    lx=limits(3);
    ly=limits(1);
    lz=limits(2);
end

% elseif length(limits) == 6;
%
%     if strcmp(dim,'xy');
%         lx=limits(4)-limits(1);
%         ly=limits(5)-limits(2);
%         lz=limits(6)-limits(3);
%     elseif strcmp(dim,'xz');
%         lx=limits(4)-limits(1);
%         ly=limits(6)-limits(3);
%         lz=limits(5)-limits(2);
%     elseif strcmp(dim,'yz');
%         lx=limits(6)-limits(3);
%         ly=limits(4)-limits(1);
%         lz=limits(5)-limits(2);
%     end
%
% end

if nM > 0;
    
    X_distance = 20;
    A = (lx)*(ly);
    
    % elseif nM > 4 && nM < 10;
    %     X_distance = lx/3;
    % elseif nM > 9 && nM < 17;
    %     X_distance = lx/4;
    % elseif lx > 0 && ly > 0 && nM >= 11;
    
    % X_distance = floor((Area^0.5-2*Area^0.5/nM)/(nM^0.5))
    % end
    
    X_distance = (A/(nM^0.5+1)^2)^0.5
    
    XYZ_data=zeros(nM,3);
    Y_distance = X_distance;
    xCoord = [X_distance-2:X_distance:lx-2];
    Y_distance = ly/ceil(nM/size(xCoord,2));
    yCoord = [Y_distance-2:Y_distance:ly-2];
    
    i=1;
    while size(xCoord,2)*size(yCoord,2) < nM
        X_distance = X_distance - 1;
        Y_distance = Y_distance - 1;
        xCoord = [X_distance-2:X_distance:lx-2];
        yCoord = [Y_distance-2:Y_distance:ly-2];
        i = 1 +1;
    end
    
    M_count = 1;
    i=1;
    while (i <= length(yCoord)) && (M_count <= nM);
        j=1;
        while j <= length(xCoord) && (M_count <= nM);
            XYZ_data(j+(i-1)*length(xCoord),:) = [xCoord(j),yCoord(i),lz];
            XYZ_labels(j+(i-1)*length(xCoord),1) = {atom_label};
            j = j + 1;
            M_count = M_count + 1;
        end
        i = i + 1;
    end
    
    %     %
    %     nRepIon=0;
    %     for i=1:size(XYZ_data,1);
    %         %  if XYZ_data(i,2) > 100 && nRepIon < (nR) || XYZ_data(i,2) < 20 && nRepIon < (nR);
    %         if mod(i,2) && nRepIon < (nR);
    %             XYZ_labels(i,1) = {Replacement};
    %             nRepIon=nRepIon+1;
    %         end
    %     end
    %
    %     nRepIon;
    
else
    
    XYZ_data=[];
    XYZ_labels=[];
    
end

if strcmp(dim,'xy')
    XYZ_data=XYZ_data;
elseif strcmp(dim,'xz')
    XYZ_data=[XYZ_data(:,1) XYZ_data(:,3) XYZ_data(:,2)];
elseif strcmp(dim,'yz')
    XYZ_data=[XYZ_data(:,3) XYZ_data(:,1) XYZ_data(:,2)];
end

% if length(limits) == 6;
%     if strcmp(dim,'xy');
%         XYZ_data(:,1)=XYZ_data(:,1)+limits(1);
%         XYZ_data(:,2)=XYZ_data(:,2)+limits(2);
%     elseif strcmp(dim,'xz');
%         XYZ_data(:,1)=XYZ_data(:,1)+limits(1);
%         XYZ_data(:,3)=XYZ_data(:,3)+limits(3);
%     elseif strcmp(dim,'yz');
%         XYZ_data(:,2)=XYZ_data(:,3)+limits(2);
%         XYZ_data(:,3)=XYZ_data(:,3)+limits(3);
%     end
% end

assignin('base','XYZ_data',XYZ_data)
assignin('base','XYZ_labels',XYZ_labels)

atom = add2atom(XYZ_labels,XYZ_data,atom_label,[]);

if nargin == 5
    trans_vec=cell2mat(varargin(1));
    atom = translate_atom(atom,trans_vec,'all');
end

disp('atom grid created!')

% XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
% XYZ_labels=[atom.type]';




