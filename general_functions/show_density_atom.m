%% show_density_atom.m
% * This function is used to plot density profiles in the X|Y|Z-direction
% * It's fancier than plot_density_atom, but slower..
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # show_density_atom(atom,Box_dim)
% # show_density_atom(atom,Box_dim,0.1)
% # show_density_atom(atom,Box_dim,0.1,2)
% # show_density_atom(atom,Box_dim,0.1,2,[1:10:100])

function show_density_atom(atom,Box_dim,varargin)

hold off

if nargin > 2
    if varargin{1}==0
        disp('Will not draw any box!')
    else
        try
            Simbox = draw_box_atom(Box_dim,[0 0 0.8],2);
        catch
            disp('Could draw the box!')
        end
    end
end


show_atom(atom,Box_dim,'ballstick');

if nargin>3
    scalehist=varargin{2};
else
    scalehist=1;
end

scalehist=1/(0.01*scalehist);

if nargin>4
    ind=varargin{3};
else
    ind=1:size(atom,2);
end

[dx,x,dy,y,dz,z] = hist_atom(atom,Box_dim,.2,[0 0 0],0,1); % binsize,[center_x center_y center_z],symmetric,gaussian deconv)

hold on

if numel(unique([atom.type]))==1
    color = element_color([atom(1).type]);
    plot3(dx,Box_dim(1,2)*ones(numel(x),1),x/scalehist+Box_dim(1,3)+1,'Color',color,'LineWidth',2);
    plot3(zeros(numel(y),1),dy,y/scalehist+Box_dim(1,3)+1,'Color',color,'LineWidth',2);
    plot3(z/scalehist+Box_dim(1,1)+1,Box_dim(1,2)*ones(numel(z),1),dz,'Color',color,'LineWidth',1);
else
    plot3(dx,Box_dim(1,2)*ones(numel(x),1),x/scalehist+Box_dim(1,3)+1,'r','LineWidth',2);
    plot3(zeros(numel(y),1),dy,y/scalehist+Box_dim(1,3)+1,'b','LineWidth',2);
    plot3(flipud(z)/scalehist+Box_dim(1,1)+1,Box_dim(1,2)+ones(numel(z),1),dz,'k','LineWidth',1);
end
axis([-5 ceil(max(z/scalehist+Box_dim(1,1)+1)/10)*10 -5 ceil(max(z/scalehist+Box_dim(1,2)+1)/10)*10 -5 ceil(max(x/scalehist+Box_dim(1,3)+1)/10)*10],'equal');

fig = gcf;fig.Color = [1 1 1];
set(gca,'Color',[1 1 1]);
xlabel('X [�]'); ylabel('Y [�]'); zlabel('Z [�]');
