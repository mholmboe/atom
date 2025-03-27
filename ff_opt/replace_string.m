%% replace_string.m
% * This special function is used to replace strings in text files.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # replace_string('old string','new string','text.dat')
% # replace_string('old string','new string','text.dat','newtext.dat')

function replace_string(oldstring,newstring,filename,varargin)

if isnumeric(oldstring)
    oldstring=char(oldstring);
end

if isnumeric(newstring)
    newstring=char(newstring);
end

if iscell(oldstring)
    oldstring=char(oldstring);
end
   
if iscell(newstring)
    newstring=char(newstring);
end
    
fid  = fopen(filename,'r');
f=fread(fid,'*char')';
fclose(fid);

f = strrep(f,oldstring,newstring);

if nargin==3
    outfilename=filename;
else
    outfilename=varargin{1};
end

if strcmp(outfilename,filename)
    bckp=strcat(filename,'_bckp')
    movefile(filename,bckp)
end

fid  = fopen(outfilename,'w');
fprintf(fid,'%s',f);
fclose(fid);

delete(bckp);