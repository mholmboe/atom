%% To run DFTD3 from MATLAB
% Add your own path to DFTD3 here, if you have installed it separately
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
% # dftd3('POSCAR -pbc -zero')

function dftd3(varargin)

if nargin>0

    eval(char(strcat({'!/Users/miho0052/opt/anaconda3/bin/dftd3'},{' '},varargin{1})))
    
else
    
    eval('!/Users/miho0052/opt/anaconda3/bin/dftd3 -h')
    
end