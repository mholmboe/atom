%% write_xyz.m
% * This function writes an xyz file from the XYZ_labels and XYZ_data
% matrix
%
%% Version
% 2.11
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # write_xyz(XYZ_labels,XYZ_data,filename_out) % Basic input arguments
% # write_xyz(XYZ_labels,XYZ_data,Box_dim,filename_out) % Basic input arguments, with Box_dim as a comment
% # write_xyz(XYZ_labels,XYZ_data,[],filename_out,Comment) % Basic input arguments, with a comment string
%
function write_xyz(XYZ_labels,XYZ_data,varargin)

if nargin==3
    filename_out=varargin{1};
    Comment='Written in MATLAB by MHolmboe'
elseif nargin==4
    Box_dim=varargin{1};
    filename_out=varargin{2};
elseif nargin>4
    
    filename_out=varargin{2};
    Comment=varargin{3};
end


if regexp(filename_out,'.xyz') ~= false
    filename_out = filename_out;
else
    filename_out = strcat(filename_out,'.xyz');
end

nAtoms=size(XYZ_labels,1)
fid = fopen(filename_out, 'wt');
fprintf(fid, '%-5i\r\n',nAtoms);

if exist('Box_dim','var')
    if numel(Box_dim)==1
        Box_dim(1)=Box_dim(1);
        Box_dim(2)=Box_dim(1);
        Box_dim(3)=Box_dim(1);
    end
    
    if length(Box_dim)==3
        fprintf(fid, '# %10.5f%10.5f%10.5f\r\n',Box_dim);
    elseif length(Box_dim)==6
        fprintf(fid, '# %10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\r\n',Box_dim);
    elseif length(Box_dim)==9
        fprintf(fid, '# %10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\r\n',Box_dim);
    end
else
    fprintf(fid, '%s\r\n',char(Comment{1}));
end

for i = 1:nAtoms
    Atom_section(1:4) = [XYZ_labels(i), XYZ_data(i,1), XYZ_data(i,2), XYZ_data(i,3)];
    fprintf(fid, '%-5s%10.5f%10.5f%10.5f\r\n', Atom_section{1:4});
end

fprintf(fid, '\r\n');

fclose(fid);
disp('.xyz structure file written')
