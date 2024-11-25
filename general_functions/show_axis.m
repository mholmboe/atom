%% show_axis.m
% * This function draws the axis in a plot
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # show_axis(varargin)
% # show_axis(axis_origin,len)
% # show_axis(axis_origin,len,text)
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
    len=5;
else
    len=varargin{2};
end

if length(len)==1
    len=[len len len];
end

xyz_labels=0;
if nargin>2 && varargin{3}>0
    xyz_labels=1;
else
    xyz_labels=0;
end

transparancy=0.2;
if nargin>3 && varargin{4}>0
    transparancy=varargin{4};
end


if size(len,2)==3 || size(len,2)==9

    Lx = len(1);
    Ly = len(2);
    Lz = len(3);
    if numel(len)==3
        xy = 0;
        xz = 0;
        yz = 0;
    else
        xy = len(6);
        xz = len(8);
        yz = len(9);
    end

    a = Lx;
    b = (Ly^2 + xy^2)^0.5;
    c = (Lz^2 + xz^2 + yz^2)^0.5;
    alfa = acosd((xy*xz+Ly*yz)/(b*c));
    beta = acosd(xz/c);
    gamma = acosd(xy/b);

elseif size(len,2)==6

    disp('Using angles')
    a=len(1);
    b=len(2);
    c=len(3);
    alfa=len(4);
    beta=len(5);
    gamma=len(6);
    Lx = a;
    xy = b * cos(deg2rad(gamma));
    Ly = (b^2-xy^2)^.5;
    xz = c*cos(deg2rad(beta));
    yz = (b*c*cos(deg2rad(alfa))-xy*xz)/Ly;
    Lz = (c^2 - xz^2 - yz^2)^0.5;

end

hold on;

x = show_arrow(axis_origin,[axis_origin(1)     axis_origin(2)     axis_origin(3)+len(1)],'color',[0 0 1],'facealpha',transparancy);
y = show_arrow(axis_origin,[axis_origin(1)     axis_origin(2)+len(2) axis_origin(3)    ],'color',[0 1 0],'facealpha',transparancy);
z = show_arrow(axis_origin,[axis_origin(1)+len(3) axis_origin(2)     axis_origin(3)    ],'color',[1 0 0],'facealpha',transparancy);

hold on;

x = show_arrow(axis_origin,[axis_origin(1)     axis_origin(2)     axis_origin(3)+len(1)],'color',[0 0 1],'facealpha',transparancy);
y = show_arrow(axis_origin,[axis_origin(1)     axis_origin(2)+len(2) axis_origin(3)    ],'color',[0 1 0],'facealpha',transparancy);
z = show_arrow(axis_origin,[axis_origin(1)+len(3) axis_origin(2)     axis_origin(3)    ],'color',[1 0 0],'facealpha',transparancy);

xlim([-5+axis_origin(1) 15+axis_origin(1)+len(1)]);
ylim([-5+axis_origin(2) 15+axis_origin(2)+len(2)]);
zlim([-5+axis_origin(3) 15+axis_origin(3)+len(3)]);


if xyz_labels>0
    text(axis_origin(1)+1.05*len(1),axis_origin(2),axis_origin(3),'X','FontSize',16);
    text(axis_origin(1),axis_origin(2)+1.05*len(2),axis_origin(3),'Y','FontSize',16);
    text(axis_origin(1),axis_origin(2),axis_origin(3)+1.05*len(3),'Z','FontSize',16);
end

rotate3d on;

