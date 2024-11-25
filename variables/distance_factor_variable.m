%% distance_factor
% * This variable is used to find bonds or nearest neighbours, and is used 
% in combination with atomtypes crystal or vdw radii, where typical values 
% are 1.2 and 0.6, respectively.
% 
% * Explanation: If the sum of two atomtypes radii, multiplied by the 
% distance_factor is larger than the distance between two atoms, the atoms 
% are considered to be bonded/neighbours.
%
%% Version
% 3.00
%

%% Example
% distance_factor = 1.2 % Typical value when using crystal radii
distance_factor = 1.2