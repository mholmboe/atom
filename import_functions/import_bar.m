%% import_bar.m
% * This function imports the type of .bar files that the MD package
% Gromacs uses for some of its text-based data output
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% #  Data = import_bar('energy.bar')
% #  Data = import_bar('energy.bar','plot')
%
function Data = import_bar(filename,varargin)

% Open and read the input file
inputfile = fopen(filename, 'r');
C = textscan(inputfile, '%s', 'Delimiter', '\n');
fclose(inputfile);

nRows = size(C{1,1},1);
nColumns=size(strsplit(char(C{1,1}(end-1,:))),2);

% Search a specific string and find all rows containing matches
results = strfind(C{1}, 'point');
results_rows = find(~cellfun('isempty', results));
results_set=[1;diff(results_rows)];
results_data=C{1,1}(results_rows);
results_data = regexprep(results_data,'point','');
results_data = regexprep(results_data,'DG','');
results_data = regexprep(results_data,',','');
results_data = regexprep(results_data,'+/-','');

Data=[];
for i=1:size(results_data,1)
    temp=regexp(results_data{i,:}, '\s+', 'split');
    temp_data = cellfun(@str2double,temp);
    if results_set(i)==1
        Data(i,1:2)=temp_data(end-1:end);
    else
        results_set(i:end)=2;
        Data(i-size(Data,1),3:4)=temp_data(end-1:end);
    end
end

if results_set(end)>1
    Data(:,3:4)=flipud(Data(:,3:4));
    Data(:,3)=-Data(:,3);
end

% Search a specific string and find all rows containing matches
total = strfind(C{1}, 'total');
total_rows = find(~cellfun('isempty', total));
total_set=[1;diff(total_rows)];
total_data=C{1,1}(total_rows);
total_data = regexprep(total_data,'total','');
total_data = regexprep(total_data,'DG','');
total_data = regexprep(total_data,',','');
total_data = regexprep(total_data,'+/-','');

Total_Data=[];Data=[Data;Data(end,:)];
for i=1:size(total_data,1)
    temp=regexp(total_data{i,:}, '\s+', 'split');
    temp_data = cellfun(@str2double,temp);
    j=2*(i-1)+1;
    Data(end,j:j+1)=temp_data(end-1:end);
end



% If we quickly want to plot the data, add a second argument when calling the function
if nargin>1
    hold on
    plot(Data(1:end-1,1),'LineWidth',2)
    if size(Data,2)>2
        plot(Data(1:end-1,3),'LineWidth',2)
    end
    set(gcf,'color','w');%,'units','normalized','position',[0,0,.4,.6]);
    %     set(gca, 'FontName', 'Arial','FontSize',22,'TickDir','out','Ytick',min(Data(:,2:end)):ceil(max(Data(:,2:end))/10)*10/5:ceil(max(Data(:,2:end))/10)*10)
    set(gca,'LineWidth',2,'FontName', 'Arial','FontSize',22,'TickDir','out')%,'Xtick',floor(Data(1,1)/10)*10:ceil(Data(end,1)/10)*10/5:ceil(Data(end,1)/10)*10)
    xlabel('vdw - q','FontSize',24);
    ylabel('kJ / mol','FontSize',24);
    xlim('auto')
end

assignin('caller','Data',Data);


