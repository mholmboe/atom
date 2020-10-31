function PATH2GMX = PATH2GMX() 

% Add your own path to GROMACS here...
% PATH2GMX ='/usr/local/gromacs-2018/bin/'; 
PATH2GMX ='/usr/local/gromacs-2016.2/bin/'; 

% % KEBNEKAISE
% PATH2GMX ='/hpc2n/eb/software/MPI/GCC/6.4.0-2.28/impi/2018.1.163/GROMACS/2018/bin/'; 

% Does this work for the 'base' ?
PATH=getenv('PATH');

PATH=strrep(PATH,'PATH2GMX','');

setenv('PATH', [getenv('PATH'),':',PATH2GMX]);

% assignin('caller','PATH2GMX',PATH2GMX);

% Now you could try
% system('gmx editconf -h')

end