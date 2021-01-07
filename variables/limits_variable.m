%% limits
% * This variable is a 1x1 or 1x3 or 1x6 vector, representing an orthogonal 
% volume. If a 1x1 vector, the volume is a cube expanding from (0,0,0) to
% (limits(1),limits(1),limits(1)). If a 1x3 vector, the volume expands from
% (0,0,0) to (limits(1),limits(2),limits(3)). If a 1x6 vector, the volume
% is a region defined by <limits.html limits>, as in 
% [xlo ylo zlo xhi yhi zhi]. Note that a 1x6 vector can alo be used to
% repsresent any arbitrary plane, as in the last example below where
% zlo equals zhi.
%
%% Version
% 2.082
%
%% Example
% # limits = [10]; % a cube from (0,0,0) to (10,10,10)
% # limits = [10 20 30];
% # limits = [10 10 0 30 40 20];
% # limits = [10 10 0 30 40 0]; % a xy plane from (10,10) to (30,40) at z=0
limits = [10 10 0 30 40 20];