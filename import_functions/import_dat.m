%% import_dat.m
% * This function imports the text data files that has non-data on the 
% first line
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% #  Data = import_dat('energy.dat')
% #  Data = import_dat('energy.dat','plot')
%
function Data = import_dat(filename,varargin)

% Check if .dat is given in the filename
if regexp(filename,'.dat') ~= false
    filename = filename;
else
    filename = strcat(filename,'.dat');
end

% Open and read the .dat file
inputfile = fopen(filename, 'r');
C = textscan(inputfile, '%s', 'Delimiter', '\n');
fclose(inputfile);

nRows = size(C{1,1},1);
nColumns=size(strsplit(char(C{1,1}(end-1,:))),2);

% Search a specific string and find all rows containing matches
xaxis = strfind(C{1}, '@    xaxis  label ');
xaxislabel_row = find(~cellfun('isempty', xaxis));
xaxislabel=C{1,1}(xaxislabel_row);
xaxislabel = regexprep(xaxislabel,'@    xaxis  label "','');
xaxislabel = regexprep(xaxislabel,'"','');

yaxis = strfind(C{1}, '@    yaxis  label ');
yaxislabel_row = find(~cellfun('isempty', yaxis));
yaxislabel=C{1,1}(yaxislabel_row);
yaxislabel = regexprep(yaxislabel,'@    yaxis  label "','');
yaxislabel = regexprep(yaxislabel,'"','');
yaxislabel=strrep(yaxislabel,'\S-','-^');
yaxislabel=strrep(yaxislabel,'-^','^-^');
yaxislabel=strrep(yaxislabel,'\S','^');
yaxislabel=strrep(yaxislabel,'\N','');

yaxis_all_legends=[];
for i=1:nColumns-1
    yaxis = strfind(C{1},strcat('@ s',num2str(i-1),' legend '));
    yaxis_legend_row = find(~cellfun('isempty', yaxis));
    yaxis_legend=C{1,1}(yaxis_legend_row);
    yaxis_legend=strrep(yaxis_legend,'-','_');
    yaxis_legend = regexprep(yaxis_legend,strcat('@ s',num2str(i-1),' legend '),'');
    yaxis_all_legends = [yaxis_all_legends strtrim(regexprep(yaxis_legend,'"',''))];
end

% Option 1 (the faster way of doing it...)
D=C;

yaxis_all_legends=C{1,1}([1]);
yaxis_all_legends=strsplit(char(yaxis_all_legends));

nondatarows1 = strfind(C{1}, '#');
nondatarows1 = find(~cellfun('isempty',nondatarows1));
nondatarows2 = strfind(C{1}, '@');
nondatarows2 = find(~cellfun('isempty',nondatarows2));
nondatarows3 = strfind(C{1}, '&');
nondatarows3 = find(~cellfun('isempty',nondatarows3));
D{1,1}([1])=[];

D = regexp(D{1,1}, '\s+', 'split');
D = vertcat(D{:});
Data = cellfun(@str2double,D);
Data=[[1:size(Data,1)]' Data];

try
    
    for i=1:nColumns
        if i==1
            assignin('caller','x',Data(:,i));
        else
            assignin('caller',strcat('y',num2str(i-1)),Data(:,i));
            yaxis_all_legends(i-1)=strtrim(yaxis_all_legends(i-1));
            yaxis_all_legends(i-1)=strrep(yaxis_all_legends(i-1),' ','_');
            yaxis_all_legends(i-1)=strrep(yaxis_all_legends(i-1),'_',' ');
            yaxis_all_legends(i-1)=strrep(yaxis_all_legends(i-1),'.','');
            yaxis_all_legends(i-1)=strrep(yaxis_all_legends(i-1),'{}','');
            yaxis_all_legends(i-1)=strrep(yaxis_all_legends(i-1),'\','');
            yaxis_all_legends(i-1)=strrep(yaxis_all_legends(i-1),'/','');
            yaxis_all_legends(i-1)=strrep(yaxis_all_legends(i-1),'=','_');
            yaxis_all_legends(i-1)=strrep(yaxis_all_legends(i-1),',','_');
            yaxis_all_legends(i-1)=strrep(yaxis_all_legends(i-1),'(','_');
            yaxis_all_legends(i-1)=strrep(yaxis_all_legends(i-1),')','_');
            assignin('caller',char(yaxis_all_legends(i-1)),Data(:,i));
        end
    end
catch
    try
        yaxis_all_legends=strrep(yaxis_all_legends,'_',' ');
    catch
    end
    disp('Could not read the column labels, but you')
    disp('still have your data in the variable Data!')
end

yaxis_all_legends

% If we quickly want to plot the data, add a second argument when calling the function
hold on
if nargin>1
    if nargin>2
        plotColumns=varargin{2};
    else
        plotColumns=1:nColumns;
    end
    for i=1:nColumns
        if ismember(i,plotColumns)
            plot(Data(:,1),Data(:,i+1),'LineWidth',2)
            xlabel(xaxislabel,'FontSize',24);
            ylabel(yaxislabel,'FontSize',24);
            legend(yaxis_all_legends(i))
        end
    end
    set(gcf,'color','w');%,'units','normalized','position',[0,0,.4,.6]);
    %     set(gca, 'FontName', 'Arial','FontSize',22,'TickDir','out','Ytick',min(Data(:,2:end)):ceil(max(Data(:,2:end))/10)*10/5:ceil(max(Data(:,2:end))/10)*10)
    set(gca,'LineWidth',2,'FontName', 'Arial','FontSize',22,'TickDir','out','Xtick',floor(Data(1,1)/10)*10:ceil(Data(end,1)/10)*10/5:ceil(Data(end,1)/10)*10)
    
end
hold off

assignin('caller','Data',Data);
assignin('caller','xaxislabel',xaxislabel);
assignin('caller','yaxislabel',yaxislabel);
assignin('caller','yaxis_all_legends',yaxis_all_legends);


