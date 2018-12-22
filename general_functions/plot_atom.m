%%  plot_atom
% * This function draws the atom struct in 3D. Its very basic. 
% * It can display bonds over the pbc if supplied with a Bond_index from 
% bond_angle_atom(). It can also indicate the axis directions.
% * For fancier plots, use the vmd(atom,Box_dim) function if you have VMD
% installed and the PATH2VMD() function set up accordingly
%
%% Version
% 2.0
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # plot_atom(atom,Box_dim)
% # plot_atom(atom,Box_dim,2) % 2 is a scale factor
% # plot_atom(atom,Box_dim,1,Bond_index) % Bond_index from bond_angle_atom()
% # plot_atom(atom,Box_dim,1,[],'axis') % Will indicate the axis ind=find(ismember(Bond_index(:,1),6)|ismember(Bond_index(:,2),6)),:)
% # plot_atom(atom,Box_dim,Bond_index(ind)) To see all bonds to index=6
%
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

if xhi < 10
    xhi=15;
end

if yhi < 10
    yhi=15;
end

if zhi < 10
    zhi=15;
end

% figure('units','normalized','position',[0 .2 .6 .8])
hold on
rotate3d on;
camlight(220,210,'infinite');
set(gca,'PlotBoxAspectRatio',[(xhi-xlo)/(zhi-zlo) (yhi-ylo)/(zhi-zlo) (zhi-zlo)/(zhi-zlo)],'FontSize',21);
set(gcf,'Color',[1,1,1]);


if nargin > 2
   scalefactor=20*varargin{1};
else
   scalefactor=20; 
end

% Draw bonds
if nargin > 3
    if numel(varargin{2})>0
    Bond_index=varargin{2};
    for i=1:size(Bond_index,1)
        Bond_index(i,1);
        Bond_index(i,2);
        plot3([atom(Bond_index(i,1)).x atom(Bond_index(i,2)).x],...
              [atom(Bond_index(i,1)).y atom(Bond_index(i,2)).y],...
              [atom(Bond_index(i,1)).z atom(Bond_index(i,2)).z],...
              'Color',[.5 .5 .5],...
              'LineWidth',3);
    end
    end
end

% Also draw the actual simulation box...
Simbox = draw_box_atom(Box_dim,[0 0 0],2);
maxY=max([atom.y])+1;
for i = 1:length(Atom_label)
%     disp('%%%%')
%     Atom_label(i)
    radii = scalefactor*2*abs(radius_ion(Atom_label(i)));
    color =  element_color(Atom_label(i));
%     disp('%%%%')
    ind=strncmpi([atom.type],Atom_label(i),1);
    
%     plot3([atom(ind).x],[atom(ind).y],[atom(ind).z],...
%         'o',...
%         'LineWidth',0.5,...
%         'MarkerEdgeColor',[0 0 0],...
%         'MarkerFaceColor',color,...  
%         'MarkerSize',6*radii);
    %             'MarkerFaceColor',[1 0 0],...
        scatter3([atom(ind).x],[atom(ind).y],[atom(ind).z],...
        radii,...   
        'MarkerEdgeColor',[.1 .1 .1],...
        'MarkerFaceColor',color,... 
        'MarkerFaceAlpha',0.7... 
        );
end

% Draw the axis in the lower left corner
if nargin > 4
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

fig = gcf;fig.Color = [0.85 0.85 0.85];
set(gca,'Color',[0.85 0.85 0.85]);
xlabel('X'); ylabel('Y'); zlabel('Z');
axis([xlo xhi ylo yhi zlo zhi]);
view([0,0]);

end



