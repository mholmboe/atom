function write_xvg(xdata,ydata,inputfile,varargin)
%% write_xvg.m
%% This function allows you to write a xvg file based on a previous file,
%% in so that the new file has the same header and meta data, hence the
%% only thing that changes is the actual data. Good for averaging data.
%%
%% Please report bugs to michael.holmboe@umu.se

% Check if filename is cell
if iscell(inputfile)
    inputfile=char(inputfile);
end

if size(xdata,1) > size(xdata,2)
    xdata=xdata';
end

if size(ydata,1) > size(ydata,2)
    ydata=ydata';
end

%% Check if .xvg is given in the filename
if regexp(inputfile,'.xvg') ~= false
    inputfile = inputfile;
else
    inputfile = strcat(inputfile,'.xvg');
end

if nargin > 3
    outputfile=varargin{1};

    % Check if filename is cell
    if iscell(outputfile)
        outputfile=char(outputfile);
    end

    if regexp(char(outputfile),'.xvg') ~= false
        outputfile=outputfile;
    else
        outputfile=strcat(outputfile,'.xvg');
    end
else
    outputfile=regexprep(inputfile,'.xvg','');
    outputfile=strcat(outputfile,'_new.xvg');
end

%% Open and read the .xvg file
fid = fopen(inputfile,'r');
C = textscan(fid, '%s', 'Delimiter','\n');
fclose(fid);

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
nondatarows1 = strfind(C{1}, '#');
nondatarows1 = find(~cellfun('isempty',nondatarows1));
nondatarows2 = strfind(C{1}, '@');
nondatarows2 = find(~cellfun('isempty',nondatarows2));
nondatarows3 = strfind(C{1}, '&');
nondatarows3 = find(~cellfun('isempty',nondatarows3));
D{1,1}([nondatarows1;nondatarows2;nondatarows3])=[];

D = regexp(D{1,1}, '\s+', 'split');
D = vertcat(D{:});
Data = cellfun(@str2double,D);

% End option 1

try
    for i=1:nColumns
        if i==1
            assignin('caller','x',Data(:,i));
        else
            assignin('caller',strcat('y',num2str(i-1)),Data(:,i));
            yaxis_all_legends(i-1)=strtrim(yaxis_all_legends(i-1));
            yaxis_all_legends(i-1)=strrep(yaxis_all_legends(i-1),' ','_');
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
    disp('Could not read the column labels, but you')
    disp('still have your data in the variable Data!')
end

assignin('caller','Data',Data);
assignin('caller','xaxislabel',xaxislabel);
assignin('caller','yaxislabel',yaxislabel);
assignin('caller','yaxis_all_legends',yaxis_all_legends);
row=max([nondatarows1;nondatarows2;nondatarows3])+1;

C{1,1}(row:end)=[];
% C{1,1}(row:end)=1;

fid = fopen(outputfile, 'wt');
[nrows,ncols] = size(C{1,1});
for i = 1:nrows
    fprintf(fid,'%s\r\n',char(C{1,1}(i)));
end

A=mat2cell([xdata' ydata'],[max(size(xdata))],[1+min(size(ydata))]);
for i = 1:max(size(xdata))
    fprintf(fid, '%8.5f %8.5f\r\n', A{1,1}(i,:));
end
fclose(fid);

