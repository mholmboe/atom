%% draw_box_atom.m
% * This function draws the simulation box, have not used it in a while,
% does it work with triclinic sim box?
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # draw_box_atom(Box_dim,[0.5 0.5 0.5],2)
%
function Simbox = draw_box_atom(Box_dim,varargin)

if nargin>1
    LineColor=varargin{1};
else
    LineColor=[0 0 1];
end

if nargin>2
    LineThickness=varargin{2};
else
    LineThickness=1;
end

if size(Box_dim,2)==3 || size(Box_dim,2)==9
    
    Lx = Box_dim(1);
    Ly = Box_dim(2);
    Lz = Box_dim(3);
    if numel(Box_dim)==3
        xy = 0;
        xz = 0;
        yz = 0;
    else
        xy = Box_dim(6);
        xz = Box_dim(8);
        yz = Box_dim(9);
    end
    
    a = Lx;
    b = (Ly^2 + xy^2)^0.5;
    c = (Lz^2 + xz^2 + yz^2)^0.5;
    alfa = acosd((xy*xz+Ly*yz)/(b*c));
    beta = acosd(xz/c);
    gamma = acosd(xy/b);
    Volume=a*b*c* ((1-cos(alfa*pi()/180)^2-cos(beta*pi()/180)^2-cos(gamma*pi()/180)^2)+2*(cos(alfa*pi()/180)*cos(beta*pi()/180)*cos(gamma*pi()/180)))^0.5;
    
    Simbox = [a b c alfa beta gamma Volume];
    assignin('base','Simbox',Simbox);
    %copy(Simbox);
    
    x = [0 a a+xy xy 0 nan xz a+xz a+xy+xz xy+xz   xz nan 0 xz nan xy  xy+xz nan a+xy a+xy+xz nan a a+xz];
    y = [0 0 b    b  0 nan yz yz   b+yz    b+yz yz nan 0 yz nan b   b+yz  nan b    b+yz    nan 0 yz  ];
    z = [0 0 0    0  0 nan c  c    c       c    c  nan 0 c  nan 0   c     nan 0    c       nan 0 c   ];
    hold on;
    plot3(x,y,z,'color',LineColor,'linewidth', LineThickness);
    % xlabel('X'); ylabel('Y'); zlabel('Z')
    
elseif size(Box_dim,2)==6
    
    Lx = Box_dim(4)-Box_dim(1);
    Ly = Box_dim(5)-Box_dim(2);
    Lz = Box_dim(6)-Box_dim(3);
    
    xy = 0;
    xz = 0;
    yz = 0;
    
    a = Lx;
    b = (Ly^2 + xy^2)^0.5;
    c = (Lz^2 + xz^2 + yz^2)^0.5;
    alfa = acosd((xy*xz+Ly*yz)/(b*c));
    beta = acosd(xz/c);
    gamma = acosd(xy/b);
    Volume=a*b*c* ((1-cos(alfa*pi()/180)^2-cos(beta*pi()/180)^2-cos(gamma*pi()/180)^2)+2*(cos(alfa*pi()/180)*cos(beta*pi()/180)*cos(gamma*pi()/180)))^0.5;
    
    Simbox = [a b c alfa beta gamma Volume];
    assignin('base','Simbox',Simbox);
    %copy(Simbox);
    
    x = Box_dim(1)+[0 a a+xy xy 0 nan xz a+xz a+xy+xz xy+xz   xz nan 0 xz nan xy  xy+xz nan a+xy a+xy+xz nan a a+xz];
    y = Box_dim(2)+[0 0 b    b  0 nan yz yz   b+yz    b+yz yz nan 0 yz nan b   b+yz  nan b    b+yz    nan 0 yz  ];
    z = Box_dim(3)+[0 0 0    0  0 nan c  c    c       c    c  nan 0 c  nan 0   c     nan 0    c       nan 0 c   ];
    hold on;
    plot3(x,y,z,'color',LineColor,'linewidth', LineThickness);
    % xlabel('X'); ylabel('Y'); zlabel('Z')
end






