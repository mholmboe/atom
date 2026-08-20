function DataAll = plot_all_xvg(filenamePattern,varargin)
% plot_xvg  Import and plot one or more .xvg files (with wildcard support)
%
% Usage:
%   Data = plot_xvg('energy.xvg', ... )
%   DataAll = plot_xvg('energy*.xvg', [1 5 3], [0 200], [0 1.5], 2)
%
% Inputs:
%   filenamePattern   String, may include '*' wildcard (e.g. 'ene*.xvg')
%   varargin{1}       Columns to plot (e.g. [1 5 3])
%   varargin{2}       X-axis limits [xmin xmax]
%   varargin{3}       Y-axis limits [ymin ymax]
%   varargin{4}       Smoothing span (integer for Gaussian smoothing)
%
% Output:
%   DataAll           Cell array (1×Nfiles) of imported & processed data matrices

% Find matching files
fileList = dir(filenamePattern);
 
% % Define an “except” pattern: match *.xvg but not *_test.xvg
% pat = wildcardPattern('*.xvg','Except','*_test.xvg');
% fileList = dir(pat);

if isempty(fileList)
    error('No files match pattern "%s".', filenamePattern);
end

nFiles = numel(fileList);
DataAll = cell(1,nFiles);

% Loop over each file
for k = 1:nFiles
    fname = fileList(k).name;
    fprintf('Importing and plotting %s...\n', fname);

    % Import
    Data = import_xvg(fname);  %#ok<NASGU>

    % Optional: select columns
    if numel(varargin) >= 1 && ~isempty(varargin{1})
        cols = varargin{1};
        if max(cols) <= size(Data,2)
            Data = Data(:, cols);
        else
            warning('Requested columns exceed data size; using all columns.');
        end
    end

    % Optional: apply x-limits (filter out points outside)
    if numel(varargin) >= 2 && ~isempty(varargin{2})
        xl = varargin{2};
        idx = Data(:,1) < xl(1) | Data(:,1) > xl(2);
        Data(idx,:) = [];
    end

    % Optional: apply smoothing
    if numel(varargin) >= 4 && ~isempty(varargin{4})
        span = varargin{4};
        for c = 2:size(Data,2)
            Data(:,c) = smooth(Data(:,c), span);
        end
    end

    % Store processed data
    if k==1
        DataAll = Data;
    else
        DataAll = [DataAll Data(:,2:end)];
    end
    % Plot in its own figure
    hold on;
    % figure('Name', fname, 'Color', [1 1 1]);
    plot(Data(:,1), Data(:,2:end), 'LineWidth', 1.5);
    set(gca, 'LineWidth', 2, 'FontName', 'Arial', 'FontSize', 16, 'TickDir', 'out');
    xlabel(xaxislabel, 'FontSize', 18); %# assumes import_xvg creates this var
    ylabel(yaxislabel, 'FontSize', 18);
    if numel(varargin) >= 3 && ~isempty(varargin{3})
        ylim(varargin{3});
    end
    legend(yaxis_all_legends, 'Interpreter', 'none', 'Location', 'best');

    drawnow;
end

% Assign to caller for convenience
% assignin('caller','DataAll',DataAll);
end
