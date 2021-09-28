%% The OPEN-MPI path on your computer
% Add your own path to Open-mpi here. Note that you might need a '\' to skip spaces

function PATH2MPI = PATH2MPI(varargin)

PATH2MPI ='/usr/local/bin/';

if nargin == 1
    % Does this work for the 'base' ?
    PATH=getenv('PATH');
    PATH=strrep(PATH,'PATH2MPI','');
    setenv('PATH', [getenv('PATH'),':',PATH2MPI]);
    % Now you could try
    % system('mpirun')
end

end