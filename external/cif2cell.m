%% To run CIF2CELL from MATLAB
% Add your own path to CIF2CELL here, if you have installed it separately
% Note that you might need a '\' to skip spaces, or double quotations as
% in "'PATH'"
%

function cif2cell(varargin)

if nargin>0

    eval(char(strcat({'!/Users/miho0052/opt/miniconda3/bin/cif2cell'},{' '},varargin{1})))
    
else
    
    eval('!/Users/miho0052/opt/miniconda3/bin/cif2cell -h')
    
end