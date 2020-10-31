%% XYZ_data
% * This variable is a nx3 array holding the xyz coordinates in Ångström
% corresponding to the variable [[atom.x]' [atom.y]' [atom.z]']. The
% motivation for this variable is simply to have all coordinates in the
% same variable. The coordinates can be 'put back' to the <atom.html atom>
% struct using the function <xyz2atom.html xyz2atom>


%% Example
% # XYZ_data = [1 2 3; 4 5 6 ; 7 8 9]; % X data is then 1 4 7... and
% likewise for Y and Z