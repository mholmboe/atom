%% The PYTHON path on your computer
% Add your own path to PYTHON (conda3) here.
% Note that you might need a '\' to skip spaces, or double quotations as
% in "'PATH'"

function python(varargin)

if nargin>0
    
    eval(char(strcat({'!/Users/miho0052/opt/miniconda3/bin/python'},{' '},varargin{1})))
    
else
    
    eval(char(strcat({'!/Users/miho0052/opt/miniconda3/bin/python'})))
    
end