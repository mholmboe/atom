function print_top(param,varargin)

if nargin>1
    filename=varargin{1};
else
    filename='temp.top';
end

fid = fopen(filename, 'wt'); % open a text file that we can write into
fprintf(fid, '%s % s\r\n',';','Topology file written by MATLAB');
fprintf(fid, '\r\n');
for i=1:numel(param)
    % #define VAR1 XX
    fprintf(fid, '%s \r\n',char(strcat('#define VAR',num2str(i),{' '},'XX',{'   '})));
end
fprintf(fid, '\r\n');
fclose(fid);

if nargin>2
    includefile=varargin{2};
    fid  = fopen(includefile,'r');
    f=fread(fid,'*char')';
    fclose(fid);    
    fid  = fopen(filename,'a');
    fprintf(fid,'%s',f);
    fclose(fid);
end
