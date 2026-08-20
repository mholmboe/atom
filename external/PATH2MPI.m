%% The OPEN-MPI path on your computer
% Add your own path to Open-mpi here. Note that you might need a '\' to 
% skip spaces, or double quotations as in "'PATH'"
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # PATH2MPI

function PATH2MPI = PATH2MPI(varargin)

PATH2MPI ='/opt/homebrew/bin/';

if nargin == 1
    % Does this work for the 'base' ?
    PATH=getenv('PATH');
    if regexp(PATH,strcat(':',PATH2MPI))
        PATH=strrep(PATH,strcat(':',PATH2MPI),'');
        setenv('PATH',PATH);
    end
    setenv('PATH', [getenv('PATH'),':',PATH2MPI]);
    % Now you could try
    % system('mpirun')
end

end