%% The GROMACS path on your computer
% Add your own path to a mpi version of GROMACS here.
% Note that you might need a '\' to skip spaces, or double quotations as
% in "'PATH'"
%
%% Version
% 2.11
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # PATH2GMX        % Can be used to invoke Gromacs utilities, see for instance the gmx() function
% # PATH2GMX('add') % adds Gromacs to the MATLAB path

function PATH2GMX = PATH2GMX(varargin)

PATH2GMX ='/usr/local/gromacs-2021/bin'; % Note, this is a version compiled wth mpi

if nargin>0
    
    % Does this work for the 'base' ?
    PATH=getenv('PATH');
    if regexp(PATH,strcat(':',PATH2GMX))
        PATH=strrep(PATH,strcat(':',PATH2GMX),'');
        setenv('PATH',PATH);
    end
    setenv('PATH', [getenv('PATH'),':',PATH2GMX]);
    % Now you could try
    % system('gmx editconf -h')
    
end

end
