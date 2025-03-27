%% replace_row.m
% * This special function is used to replace rows in text files.
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # replace_row('old string','new string','text.dat','newtext.dat')
%

function replace_row(old_row,new_row,filename,outfilename)

if iscell(old_row)
    old_row=char(old_row);
end

if iscell(new_row)
    new_row=char(new_row);
end

% Open and read the topology file
fid = fopen(filename, 'r');
C = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
% nRows = size(C{1,1},1);
% nColumns=size(strsplit(char(C{1,1}(end-1,:))),2);
ffparams_rows = strfind(C{1}, old_row); % all the #define VARX ...
ffparams_rows = find(~cellfun('isempty', ffparams_rows));
% ffparams=C{1,1}(ffparams_rows);
C{1,1}(ffparams_rows(1))={new_row};

fid = fopen(outfilename,'w');
CT = C{1,1}.';
fprintf(fid,'%s\n', CT{:});
fclose(fid);

end

