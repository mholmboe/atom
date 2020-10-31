%% The GROMACS path on your computer
% Add your own path to GROMACS here. Note that you might need a '\' to skip spaces

function PATH2GMX = PATH2GMX() 

PATH2GMX ='/usr/local/gromacs-2018/bin/'; 
% PATH2GMX ='/usr/local/gromacs-2016.2/bin/'; 
% PATH2GMX ='/usr/local/gromacs-5.1.1/bin/';

% Does this work for the 'base' ?
PATH=getenv('PATH');
PATH=strrep(PATH,'PATH2GMX','');
setenv('PATH', [getenv('PATH'),':',PATH2GMX]);
% Now you could try
% system('gmx editconf -h')
end