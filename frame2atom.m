%% frame2atom.m
% * This function extracts a frame to the trajectory matrix
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * atom = frame2atom(atom,traj,frame,Box_dim)
% * atom = frame2atom(atom,traj,frame,Box_dim,'smooth')

function new_atom = frame2atom(atom,traj,frame,Box_dim,varargin)
%%
new_atom=atom;
XYZ=traj(frame,:);
new_XYZ_data=reshape(XYZ,3,[])';

% To smooth out .xtc data due to low precision
if nargin>4
    new_XYZ_data=new_XYZ_data+(rand(size(new_XYZ_data))-.5)/100;
end

for i=1:size(new_XYZ_data,1);
    new_atom(i).x=new_XYZ_data(i,1);
    new_atom(i).y=new_XYZ_data(i,2);
    new_atom(i).z=new_XYZ_data(i,3);
end

% new_atom=update_atom(new_atom);

% function atom = atom_frame(atom,frame)
% %% This function is used to extract a frame to the atom struct
%
% X=num2cell(arrayfun(@(x) x.x(frame,:), atom));
% Y=num2cell(arrayfun(@(x) x.y(frame,:), atom));
% Z=num2cell(arrayfun(@(x) x.z(frame,:), atom));
%
% [atom.x]=deal(X{:});
% [atom.y]=deal(Y{:});
% [atom.z]=deal(Z{:});

%assignin('caller','frame',frame);