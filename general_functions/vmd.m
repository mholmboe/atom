%% vmd.m
% * This function plots the atom struct in VMD. If passing a filename as a 
% function argument you can likley use any type of fileformat VMD can
% handle.
% * Make sure to edit your own PATH2VMD function!!!
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # vmd(atom)
% # vmd(atom,Box_dim)
% # vmd('filename.gro')
% # vmd('filename.pdb')
% # vmd('filename.pdb','filename.xtc')

function vmd(varargin)

% You need to set your own VMD path here!!!
PATH2VMD = '/Applications/VMD\ 1.9.4a38.app/Contents/MacOS/startup.command';
           '/Applications/VMD 1.9.4a38.app/Contents/MacOS'

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

