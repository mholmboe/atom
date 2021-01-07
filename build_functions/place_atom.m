%% place_atom.m
% * This function places the atom struct according to the position vector called position, trying to use the COM of the molecule
%
%% Version
% 2.09
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Function arguments
% * position is a [1x3] vector
%
%% Examples
% * atom = place_atom(atom,position)
%
function atom=place_atom(atom,position)

if size(atom,2) < 500
    atom = COM_atom(atom); % This generates the COM position through assignin. You should use unwrapped molecules
    
    atom=rmfield(atom,'COM_x');
    atom=rmfield(atom,'COM_y');
    atom=rmfield(atom,'COM_z');
    atom=rmfield(atom,'element');
    atom=rmfield(atom,'Mw');
    
else
    COM=[mean([atom.x]) mean([atom.y]) mean([atom.z])]; % Since COM_atom is a bit slow for large molecules, we do this for big molecules
end

atom = translate_atom(atom,-COM+position,'all');
