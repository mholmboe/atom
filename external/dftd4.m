%% To run DFTD4 from MATLAB
% Add your own path to DFTD4 here, if you have installed it separately
% Note that you might need a '\' to skip spaces, or double quotations as
% in "'PATH'"
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # dftd4('POSCAR -c 0 -2')
% 

function dftd4(varargin)

if nargin>0
    
    eval(char(strcat({'!/Users/miho0052/opt/miniconda3/bin/dftd4'},{' '},varargin{1})))
    
else
    
    eval('!/Users/miho0052/opt/miniconda3/bin/dftd4 -h')
    
end