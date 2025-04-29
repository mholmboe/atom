%% The Octave path on your computer
% Add your own path to Octave here.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # PATH2OCTAVE        % Can be used to invoke Octave utilities
% # PATH2OCTAVE('add') % adds Octave to the MATLAB path

function PATH2OCTAVE = PATH2OCTAVE(varargin)

PATH2OCTAVE ='/opt/homebrew/Cellar/octave/9.3.0/bin/octave-cli-9.3.0';

if nargin>0
    
    PATH=getenv('PATH');
    if regexp(PATH,strcat(':',PATH2GMX))
        PATH=strrep(PATH,strcat(':',PATH2GMX),'');
        setenv('PATH',PATH);
    end
    setenv('PATH', [getenv('PATH'),':',PATH2GMX]);
    
end

end
