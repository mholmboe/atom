function change_top(param,varargin)

if nargin>1
    topfilename=varargin{1};
else
    topfilename='topol.top';
end

%% Open and read the topology file
inputfile = fopen(topfilename, 'r');
C = textscan(inputfile, '%s', 'Delimiter', '\n');
fclose(inputfile);
nRows = size(C{1,1},1);
nColumns=size(strsplit(char(C{1,1}(end-1,:))),2);
ffparams_rows = strfind(C{1}, '#define VAR'); % all the #define VARX ...
ffparams_rows = find(~cellfun('isempty', ffparams_rows));
orig_ffparams=C{1,1}(ffparams_rows);

new_ffparams={};
for i=1:numel(param)
    define_string=strcat('#define VAR',num2str(i));
    define_string=strcat(define_string,{' '}',num2str(param(i),6));
    new_ffparams(i,1)=define_string;
end

copyfile(topfilename,'temp.top')

for i=1:numel(param)
    %     replace_string(orig_ffparams(i),new_ffparams(i),'temp.top','temp2.top')
    replace_row(orig_ffparams(i),new_ffparams(i),'temp.top','temp2.top')
    copyfile('temp2.top','temp.top')
    
end

if nargin>2
    copyfile('temp.top','last.top');
end

!rm temp2.top

end