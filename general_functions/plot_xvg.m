%% plot_xvg.m
% * This function imports and plots the type of .xvg files that the MD 
% package Gromacs uses for some of its text-based data output. Note that 
% the function can output the Data and optionally plot only certain data 
% the data columns given as a second argument.
% Note also that the x and y ranges can also be set, see examples below.
%
%% Version
% 2.082
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% #  Data = plot_xvg('energy.xvg')
% #  Data = plot_xvg('energy.xvg',[1 5 3]) % Will plot the 5th and 3rd column vs. the 1st one.
% #  Data = plot_xvg('energy.xvg',[1 5 3],[0 200]) % Will also set the X-axis to be  between 0 and 200
% #  Data = plot_xvg('energy.xvg',[1 5 3],[0 200],[0 1.5]) % Will also set the Y-axis to be between 0 and 1.5
%
function Data = plot_xvg(filename,varargin)

% Check if the filename contains .xvg
if regexp(filename,'.xvg') ~= false
    filename = filename;
else
    filename = strcat(filename,'.xvg');
end

% Import the .xvg file, assigning the variables:
% Data, xaxislabel, yaxislabel and yaxis_all_legends
Data=import_xvg(filename);

% In case you specified which columns to plot
if nargin>1
    columns=varargin{1};
    if max(columns)>size(Data,2)
        disp('You wanted to plot this many columns...')
        size(columns,2)-1
        disp('...but only this many was found, except the first x-column...')
        size(Data,2)
        disp('Plotting all data columns...')
    elseif numel(columns)<1
        disp('Plotting all data columns...')
    else
        Data=Data(:,columns);
        try
        yaxis_all_legends=yaxis_all_legends(columns(2:end)-1);
        catch
            disp('Did not get the yaxis_all_legends')
        end
    end
    
end

% In case you did set the X-range to be plotted
if nargin>2
    xlimits=varargin{2};
    indlo=find(Data(:,1)<xlimits(1));
    indhi=find(Data(:,1)>xlimits(2));
    Data([indlo; indhi],:)=[];
end

% Create figure
% figure1 = figure('Color',[1 1 1]);
%plot1=plot(Data(:,1),Data(:,2:end),'LineWidth',1.5);
plot(Data(:,1),Data(:,2:end),'LineWidth',1.5);
set(gcf,'color','w');%,'units','normalized','position',[0,0,.4,.6]);
%     set(gca, 'FontName', 'Arial','FontSize',22,'TickDir','out','Ytick',min(Data(:,2:end)):ceil(max(Data(:,2:end))/10)*10/5:ceil(max(Data(:,2:end))/10)*10)
set(gca,'LineWidth',2,'FontName', 'Arial','FontSize',22,'TickDir','out');%,'Xtick',...
%floor(Data(1,1)/10)*10:ceil(Data(end,1)/10)*10/5:ceil(Data(end,1)/10)*10)

xlabel(xaxislabel,'FontSize',24);
ylabel(yaxislabel,'FontSize',24);

% In case you did set the X-range to be plotted
xlim('auto')
if nargin>3
    ylim(varargin{3});
else
    ylim('auto')
end

% Add a plot legend
legend(yaxis_all_legends,'Location','best')

% assignin('caller','plot1',plot1);
% assignin('caller','figure1',figure1);
assignin('caller','Data',Data);
assignin('caller','xaxislabel',xaxislabel);
assignin('caller','yaxislabel',yaxislabel);
assignin('caller','yaxis_all_legends',yaxis_all_legends);
