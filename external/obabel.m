%% To run OBABEL from MATLAB
% Add your own path to OBABEL here, if you have installed it separately
% Note that you might need a '\' to skip spaces, or double quotations as
% in "'PATH'"

function obabel(varargin)

if nargin>0
    
%     eval(char(strcat({'!/usr/local/bin/obabel'},{' '},varargin{1})))
    eval(char(strcat({'!/Users/miho0052/opt/miniconda3/bin/obabel'},{' '},varargin{1})))
    
else
    
%     eval('!/usr/local/bin/obabel -H')
eval('!/Users/miho0052/opt/miniconda3/bin/obabel -H')
    
end