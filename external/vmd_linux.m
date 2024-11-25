%% vmd_linux.m
% * This function plots the atom struct in VMD on linux. If passing a
% filename as a function argument you can likley use any type of fileformat
% VMD can handle.
%
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # vmd_linux(atom)
% # vmd_linux(atom,Box_dim)
% # vmd_linux('filename.gro')
% # vmd_linux('filename.pdb')
% # vmd_linux('filename.pdb','filename.xtc')

function vmd_linux(varargin)

% You need to set your own VMD path here!!!
PATH2VMD = '/usr/local/bin/vmd';


if nargin==1
    filename=varargin{1};
    if ~isstruct(filename)
        system(strcat(char({PATH2VMD()}),char(strcat({' '},{filename}))));
    else
        atom=filename;
        write_atom_gro(atom,[0 0 0],'temp_out.gro');
        system(strcat(char({PATH2VMD()}),char(strcat({' '},{'temp_out.gro'}))));
    end
else
    atom=varargin{1};
    if size(varargin{2},1)==1
    Box_dim=varargin{2};
    if numel(Box_dim)==1
        Box_dim(1)=Box_dim(1);
        Box_dim(2)=Box_dim(1);
        Box_dim(3)=Box_dim(1);
    end
    write_atom_gro(atom,Box_dim,'temp_out.gro');
    system(strcat(char({PATH2VMD()}),char(strcat({' '},{'temp_out.gro'}))));
    else
        % Assuming second argument is a traj file
        filenametraj=varargin{2};
        system(strcat(char({PATH2VMD()}),char(strcat({' '},{filename},{' '},{filenametraj}))));
    end
end

