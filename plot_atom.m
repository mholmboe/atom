%%  plot_atom
% * This function draws the atom struct in 3D. Its very simplistic with no cool features
% * whatsoever, so for fancier plots, use the vmd(atom,Box_dim) function
% * if you have VMD installed and the PATH2VMD() function setup

%% Examples
% * plot_atom(atom,Box_dim)
% * plot_atom(atom,Box_dim,'axis')

function plot_atom(atom,Box_dim,varargin)

Atom_label=unique([atom.type]);

if numel(Box_dim)==1
    Box_dim(1)=Box_dim(1);
    Box_dim(2)=Box_dim(1);
    Box_dim(3)=Box_dim(1);
end

% Sets plot limits for the data
xlo = floor(min([-2 -2+min([atom.x])])); xhi = 2 + ceil(max([max([atom.x]) Box_dim(1)])/10)*10;
ylo = floor(min([-2 -2+min([atom.y])])); yhi = 2 + ceil(max([max([atom.y]) Box_dim(2)])/10)*10;
zlo = floor(min([-2 -2+min([atom.z])])); zhi = 2 + ceil(max([max([atom.z]) Box_dim(3)])/10)*10;

if xhi < 20
    xhi=20;
end

if yhi < 20
    yhi=20;
end

if zhi < 20
    zhi=20;
end

% figure('units','normalized','position',[0 .2 .6 .8])
rotate3d on;
camlight(220,210,'infinite');
set(gca,'PlotBoxAspectRatio',[(xhi-xlo)/(zhi-zlo) (yhi-ylo)/(zhi-zlo) (zhi-zlo)/(zhi-zlo)],'FontSize',21);
set(gcf,'Color',[1,1,1]);

hold on

% Also draw the actual simulation box...
Simbox = draw_box_atom(Box_dim,[0 0 1],2);

for i = 1:length(Atom_label)
%     disp('%%%%')
%     Atom_label(i)
    radii = 2*radius_ion(Atom_label(i));
    color =  element_color(Atom_label(i));
%     disp('%%%%')
    ind=strncmpi([atom.type],Atom_label(i),1);
    
    plot3([atom(ind).x],[atom(ind).y],[atom(ind).z],...
        'o',...
        'LineWidth',0.5,...
        'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceColor',color,...  
        'MarkerSize',6*radii);
    

    %             'MarkerFaceColor',[1 0 0],...
end

% Draw the axis in the lower left corner
if nargin > 2
   plot3([xlo+1 xlo+5],[ylo+1 ylo+1],[zlo+1 zlo+1],...
        'r-',...
        'LineWidth',2); 
    plot3([xlo+1 xlo+1],[ylo+1 ylo+5],[zlo+1 zlo+1],...
        'g-',...
        'LineWidth',2); 
    plot3([xlo+1 xlo+1],[ylo+1 ylo+1],[zlo+1 zlo+5],...
        'b-',...
        'LineWidth',2); 
end

xlabel('X'); ylabel('Y'); zlabel('Z');
axis([xlo xhi ylo yhi zlo zhi]);
view([0,0]);



