function Data = import_xvg(filename)


inputfile = fopen(filename, 'r');
C = textscan(inputfile, '%s', 'Delimiter', '\n');

nRows = size(C{1,1},1);
nColumns=size(strsplit(char(C{1,1}(end-1,:))),2);

%// Search a specific string and find all rows containing matches
xaxis = strfind(C{1}, '@    xaxis  label ');
xaxis_label_row = find(~cellfun('isempty', xaxis));
xaxis_label=C{1,1}(xaxis_label_row);
xaxis_label = regexprep(xaxis_label,'@    xaxis  label "','');
xaxis_label = regexprep(xaxis_label,'"','');

yaxis = strfind(C{1}, '@    yaxis  label ');
yaxis_label_row = find(~cellfun('isempty', yaxis));
yaxis_label=C{1,1}(yaxis_label_row);
yaxis_label = regexprep(yaxis_label,'@    yaxis  label "','');
yaxis_label = regexprep(yaxis_label,'"','');

yaxis_all_legends=[];
for i=1:nColumns-1
    yaxis = strfind(C{1},strcat('@ s',num2str(i-1),' legend '));
    yaxis_legend_row = find(~cellfun('isempty', yaxis));
    yaxis_legend=C{1,1}(yaxis_legend_row);
    yaxis_legend=strrep(yaxis_legend,'-','_');
    yaxis_legend = regexprep(yaxis_legend,strcat('@ s',num2str(i-1),' legend '),'');
    yaxis_all_legends = [yaxis_all_legends strtrim(regexprep(yaxis_legend,'"',''))];
end

%% Option 1 (the faster way of doing it...)
D=C;
nondatarows1 = strfind(C{1}, '#');
nondatarows1 = find(~cellfun('isempty',nondatarows1));
nondatarows2 = strfind(C{1}, '@');
nondatarows2 = find(~cellfun('isempty',nondatarows2));
D{1,1}([nondatarows1;nondatarows2])=[];

D = regexp(D{1,1}, '\s+', 'split');
D = vertcat(D{:});
Data = cellfun(@str2double,D);

%% End option 1

% %% Option 1
% tic
% frewind(inputfile);
% i=1; j=0; allRows=nRows;
% Data=cell(nRows,1);NewData=cell(nRows,nColumns);
% while i<allRows+1;
%     tline = fgetl(inputfile);
%     if length(findstr('#',tline)) > 0 || length(findstr('@',tline)) > 0;
%         allRows=allRows-1;
%         i=i-1;
%         j=j+1;
%     else
%         Data{i,1}=strtrim(tline);
%         NewData(i,:)=strsplit(Data{i,1});
%     end
%     i=i+1;
% end
% Data = cellfun(@str2double,NewData);
% Data=Data(1:allRows,:);
% toc
% %% End option 1
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
    disp('Could not read the coulumn labels, but you')
    disp('still have your data in the Data variable!')
end


assignin('caller','Data',Data);
assignin('caller','yaxis_all_legends',yaxis_all_legends);

fclose(inputfile);



