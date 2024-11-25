%% show_box.m
% * This function draws the simulation box, have not used it in a while,
% does it work with triclinic sim box?
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # show_box(Box_dim)
% # show_box(Box_dim,[0.5 0.5 0.5])
% # show_box(Box_dim,[0.5 0.5 0.5],2)
% # show_box(Box_dim,[0.5 0.5 0.5],2,[1 2 3])
%
function Simbox = show_box(Box_dim,varargin)

if nargin>1
    LineColor=varargin{1};
    if ~size(LineColor,2)==3
        LineColor='k';
    end
else
    LineColor=[0 0 1];
end

if nargin>2
    LineThickness=varargin{2};
else
    LineThickness=1;
end

replicate=[1 1 1];
if nargin>3
    replicate=varargin{3};
    if numel(replicate)==1
        replicate(1)=replicate(1);
        replicate(2)=replicate(1);
        replicate(3)=replicate(1);
    end
    replicate(replicate==0)=1;
end

highlight=[0 0 0];
if nargin>4
    highlight=varargin{4};
    if numel(highlight)==1
        highlight(1)=highlight(1);
        highlight(2)=highlight(1);
        highlight(3)=highlight(1);
    end
    highlight(highlight==0)=1;
end

if size(Box_dim,2) == 1
    Box_dim(2)=Box_dim(1);
    Box_dim(3)=Box_dim(2);
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
    
elseif size(Box_dim,2)==6
    
    disp('Using angles')
    a=Box_dim(1);
    b=Box_dim(2);
    c=Box_dim(3);
    alfa=Box_dim(4);
    beta=Box_dim(5);
    gamma=Box_dim(6);
    Lx = a;
    xy = b * cos(deg2rad(gamma));
    Ly = (b^2-xy^2)^.5;
    xz = c*cos(deg2rad(beta));
    yz = (b*c*cos(deg2rad(alfa))-xy*xz)/Ly;
    Lz = (c^2 - xz^2 - yz^2)^0.5;
    
    Volume=a*b*c* ((1-cos(alfa*pi()/180)^2-cos(beta*pi()/180)^2-cos(gamma*pi()/180)^2)+2*(cos(alfa*pi()/180)*cos(beta*pi()/180)*cos(gamma*pi()/180)))^0.5;
    
    Simbox = [a b c alfa beta gamma Volume];
    assignin('base','Simbox',Simbox);
    
end

xmin=1000;
ymin=1000;
zmin=1000;
xmax=-1000;
ymax=-1000;
zmax=-1000;

for nz=1:replicate(3)
    for ny=1:replicate(2)
        for nx=1:replicate(1)
            
            x = [0 Lx Lx+xy xy 0 nan xz  Lx+xz Lx+xy+xz xy+xz   xz nan 0 xz nan xy  xy+xz nan Lx+xy Lx+xy+xz nan Lx Lx+xz]+(nx-1)*Lx+(ny-1)*xy+(nz-1)*xz;
            y = [0 0  Ly    Ly 0 nan yz  yz    Ly+yz    Ly+yz   yz nan 0 yz nan Ly  Ly+yz nan Ly    Ly+yz    nan 0  yz    ]+(ny-1)*Ly+(nz-1)*yz;
            z = [0 0  0     0  0 nan Lz  Lz    Lz       Lz      Lz nan 0 Lz nan 0   Lz    nan 0     Lz       nan 0  Lz    ]+(nz-1)*Lz;
            
            xmin=min([x xmin]);
            ymin=min([y ymin]);
            zmin=min([z zmin]);
            xmax=max([x xmax]);
            ymax=max([y ymax]);
            zmax=max([z zmax]);
            
            hold on;
            
            if [nx ny nz]==highlight
                plot3(x,y,z,'b-','linewidth', 2*LineThickness);
            else
                if size(LineColor,2)<3
                    plot3(x,y,z,'k--');
                else
                    plot3(x,y,z,'color',LineColor,'linewidth', LineThickness);
                end
            end
        end
    end
end

assignin('caller','xmin',xmin);
assignin('caller','ymin',ymin);
assignin('caller','zmin',zmin);
assignin('caller','xmax',xmax);
assignin('caller','ymax',ymax);
assignin('caller','zmax',zmax);

cameratoolbar
rotate3d on;
set(gcf,'Visible','on','Color',[1 1 1]);
set(gca,'Color',[1 1 1],'FontSize',24); %
xlabel('X [Å]'); ylabel('Y [Å]'); zlabel('Z [Å]');
view([0,0]);
axis normal tight equal
rotate3d on





