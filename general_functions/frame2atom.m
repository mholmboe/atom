%% frame2atom.m
% * This function extracts a frame to the trajectory matrix
%
%% Version
% 2.03
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% atom = frame2atom(atom,traj,frame,Box_dim)
% atom = frame2atom(atom,traj,frame,Box_dim,'smooth')
%
function new_atom = frame2atom(atom,traj,frame,varargin)

new_atom=atom;
XYZ=traj(frame,:);
% new_XYZ_data=reshape(XYZ,3,[])';

% To smooth out .xtc data due to low precision
if nargin>3
    XYZ=XYZ+(rand(size(XYZ))-.5)/100;
end

% This seems also to work but is 10% slower...
% X=num2cell(new_XYZ_data(:,1));
% Y=num2cell(new_XYZ_data(:,2));
% Z=num2cell(new_XYZ_data(:,3));
X=num2cell(XYZ(1,1:3:end));
Y=num2cell(XYZ(1,2:3:end));
Z=num2cell(XYZ(1,3:3:end));
[new_atom.x]=X{:};
[new_atom.y]=Y{:};
[new_atom.z]=Z{:};

% for i=1:size(new_XYZ_data,1)
%     new_atom(i).x=new_XYZ_data(i,1);
%     new_atom(i).y=new_XYZ_data(i,2);
%     new_atom(i).z=new_XYZ_data(i,3);
% end

% new_atom=update_atom(new_atom);

% function atom = atom_frame(atom,frame)
% %% This function is used to extract a frame to the atom struct
% tic;
% X=num2cell(arrayfun(@(x) x.x(frame,:), atom));
% Y=num2cell(arrayfun(@(x) x.y(frame,:), atom));
% Z=num2cell(arrayfun(@(x) x.z(frame,:), atom));
% 
% [new_atom.x]=deal(X{:});
% [new_atom.y]=deal(Y{:});
% [new_atom.z]=deal(Z{:});
% toc;

%assignin('caller','frame',frame);