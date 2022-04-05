%% The VMD path on your computer
% Add your own path to VMD here, and/or directly in the vmd() function.
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
% # PATH2VMD

function PATH2VMD = PATH2VMD(varargin)

% Example on a Mac computer
PATH2VMD = '/Applications/VMD\ 1.9.2.app/Contents/MacOS/startup.command';

% Example on a Win computer
% PATH2VMD = '"C:\Program Files (x86)\University of Illinois\VMD\vmd.exe"';

if nargin>0
    PATH=getenv('PATH');
    if regexp(PATH,strcat(':',PATH2VMD))
        PATH=strrep(PATH,strcat(':',PATH2VMD),'');
        setenv('PATH',PATH);
    end
    setenv('PATH', [getenv('PATH'),':',PATH2VMD]);
end

end

