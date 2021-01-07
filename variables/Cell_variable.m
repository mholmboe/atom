%% Cell2Box_dim.m
% * This function transforms the 1x6 Cell variable containing the a, b, c 
% cell values and  the alfa, beta, gamma angle values as used in a typical 
% .pdb file, into a 1x3 or the 1x9  Box_dim variable  
%
%% Version
% 2.082
%

%% Example
% # Cell = [10 20 30 90 90 90]; % orthogonal cell
% # Cell = [10 20 30 90 100.2 90]; % triclinic cell

%% Try this: 
% Set a Cell variable, then run:
Cell = [10 20 30 90 100.2 90];
Simbox = draw_box_atom(Cell)
