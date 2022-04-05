%% Box_dim
% * This variable is a 1x1 or 1x3 or 1x9 vector, and holds the size
% parameters of the system box/simulation cell. If <Box_dim.html Box_dim>
% is a 1x1 vector, the box is a cube expanding from (0,0,0) to
% (Box_dim(1),Box_dim(1),Box_dim(1)). If Box_dim is a 1x3 vector, the box 
% expands from (0,0,0) to (Box_dim(1),Box_dim(2),Box_dim(3)). If Box_dim is 
% a 1x9 vector, the box is a triclinic box defined by the 1 to 9 box size 
% matrix parameters [lx ly lz 0 0 xy 0 xz yz], where the latter three 
% non-zero parameters are the tilt factors of the triclinic box.
%
%% Version
% 2.11
%
%% Example
% # Box_dim = [10];
% # Box_dim = [10 20 30];
% # Box_dim = [42.435 24.245 28.452 0 0 1.142 0 2.345 1.97312];

%% Try this: 
% Set a Box_dim variable, then run:
Box_dim = [10 20 30];
Simbox = draw_box_atom(Box_dim)
