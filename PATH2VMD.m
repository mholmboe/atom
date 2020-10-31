%% The VMD path on your computer
% Add your own path to VMD here, and/or directly in the vmd() function.
% Note that you might need a '\' to skip spaces

function PATH2VMD = PATH2VMD()

PATH2VMD = '/Applications/VMD\ 1.9.2.app/Contents/MacOS/startup.command';
PATH=getenv('PATH');
PATH=strrep(PATH,'PATH2VMD','');
setenv('PATH', [getenv('PATH'),':',PATH2VMD]);

end

