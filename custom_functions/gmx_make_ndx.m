%% gmx_make_ndx.m
% * This function can help you print custom gromacs .ndx files, you
% just have to find a way to find the indexes id you want to print out
% with the group name {groupname}. The printed index file (default
% index.ndx) can be appended with new groups. See also <Make_selections.html Make_selections>
% how to make custom selections and find the corresponding atom indexes, MolID's etc.
%
% * Example on how to select all molID's of water (having [atom.types],
% the site names 'Ow') having z < 10, run this:
%
% ind=find(strcmpi([atom.type],'Ow')&[atom.z]<10); % find all Ow atoms with z coordinates < 10
%
% id=[atom(ind).molid]; % get the molID's of the Ow atoms from the previous command
%
% gmx_make_ndx(id,'SOL_low_z','index.ndx') % Write out it indexes (here called id) under the group name SOL_low_z to a file called index.ndx
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # gmx_make_ndx(id,groupname,filename) % Basic input arguments
% # gmx_make_ndx([101 102 103],'SOL') % molID's, groupname, will output a file called index.ndx
% # gmx_make_ndx([101 102 103],'SOL','SOL_ind.ndx') % molID's, groupname, filename

function gmx_make_ndx(id,groupname,varargin)

if nargin==1
    disp('You did not supply enough input arguments!')
    disp('See the examples in the function.')
elseif nargin==2
    filename='index.ndx'
elseif nargin >2
    filename=varargin{1};
end

if iscell(groupname)
    groupname=char(groupname);
end

if iscell(filename)
    filename=char(filename);
end

if regexp(filename,'.ndx') ~= false
    filename = filename;
else
    filename = strcat(filename,'.ndx');
end

% Format the ind vector to have 15 entries per row
ext_ind=zeros(1,15*ceil(length(id)/15));
ext_ind(1:length(id))=id;
ext_ind=reshape(ext_ind,15,[])';

% Print the index file
fid = fopen(filename, 'a+');
fprintf(fid, '%s %s %s\r\n','[',groupname,']');
for i = 1:size(ext_ind,1)
    row=ext_ind(i,:);
    row(row==0)=[];
    if max(row)<1000
        fprintf(fid, '%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i\r\n', row);
    elseif max(row)<10000
        fprintf(fid, '%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i\r\n', row);
    elseif max(row)<100000
        fprintf(fid, '%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i\r\n', row);
    end
end
fprintf(fid, '\r\n');
fprintf(fid, '\r\n');
fclose(fid);

end