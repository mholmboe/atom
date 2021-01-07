%% show_axis.m
% * This function draws the axis in a plot
%
%% Version
% 2.09
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # show_axis(varargin)
% # show_axis(axis_origin,len,text,transparancy)

function show_axis(varargin)


if nargin==0
    axis_origin=[-3 -3 -3];
else
    axis_origin=varargin{1};
end

if numel(axis_origin)==1
    axis_origin=[axis_origin axis_origin axis_origin];
end

if nargin<2
    len=1;
else
    len=varargin{2};
end

xyz_labels=0;
if nargin>2 && varargin{3}>0
    xyz_labels=1;
else
    xyz_labels=0;
end

transparancy=0.1;
if nargin>3 && varargin{4}>0
    transparancy=varargin{4};
end

hold on;

x = show_arrow(axis_origin,[axis_origin(1)     axis_origin(2)     axis_origin(3)+len],'color',[0 0 1],'facealpha',transparancy);
y = show_arrow(axis_origin,[axis_origin(1)     axis_origin(2)+len axis_origin(3)    ],'color',[0 1 0],'facealpha',transparancy);
z = show_arrow(axis_origin,[axis_origin(1)+len axis_origin(2)     axis_origin(3)    ],'color',[1 0 0],'facealpha',transparancy);

hold on;

x = show_arrow(axis_origin,[axis_origin(1)     axis_origin(2)     axis_origin(3)+len],'color',[0 0 1],'facealpha',transparancy);
y = show_arrow(axis_origin,[axis_origin(1)     axis_origin(2)+len axis_origin(3)    ],'color',[0 1 0],'facealpha',transparancy);
z = show_arrow(axis_origin,[axis_origin(1)+len axis_origin(2)     axis_origin(3)    ],'color',[1 0 0],'facealpha',transparancy);

xlim([-5+axis_origin(1) 15+axis_origin(1)+len]);
ylim([-5+axis_origin(2) 15+axis_origin(2)+len]);
zlim([-5+axis_origin(3) 15+axis_origin(3)+len]);


if xyz_labels>0
    text(axis_origin(1)+1.05*len,axis_origin(2),axis_origin(3),'X','FontSize',16);
    text(axis_origin(1),axis_origin(2)+1.05*len,axis_origin(3),'Y','FontSize',16);
    text(axis_origin(1),axis_origin(2),axis_origin(3)+1.05*len,'Z','FontSize',16);
end

rotate3d on;

