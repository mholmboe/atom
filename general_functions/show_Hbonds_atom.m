%%  show_hbonds_atom.m
% * This function draws the atom struct in 3D. It neglects bonds over the
% * pbc for clarity (this can be changed on line 54).
% * For less fancier plots, use the plot_atom(atom,Box_dim) function
% * This function is inspired by molecule3D.m, written by Andr? Ludwig (aludwig@phys.ethz.ch)
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # show_hbonds_atom(atom)
% # show_hbonds_atom(atom,Box_dim)
% # show_hbonds_atom(atom,Box_dim,) % representation style, should be either 'ballstick' (default),'licorice','halfvdw','vdw', 'crystal', 'ionic', 'lines', 'labels' or 'index'
% # show_hbonds_atom(atom,Box_dim,1) % 1 | 0 Will/will not show the unit cell/box
% # show_hbonds_atom(atom,Box_dim,ick',0,0.3) % Will use 30% transparency
% # show_hbonds_atom(atom,Box_dim,'ballstick',0,0,[2.25 0.6]) % Will set the rmaxlong cutoff and alternativel a distance_Factor
% # show_hbonds_atom(atom,Box_dim,'ballstick',0,0,[],[0 0 -50]) % Will translate the XYZ coordinates
% # show_hbonds_atom(atom,Box_dim,'ballstick',0,0,[],[],[0.5 0.5 0.5]) % Single color as given by the 1x3 RGB vector


function show_hbonds_atom(varargin)
tic

%% Fetch either a .pdb|.gro file or use an atom struct with its Box_dim
if ischar(varargin{1})
    filename=varargin{1};
    if regexp(filename,'.gro') > 1
        disp('Found .gro file');
        atom = import_atom_gro(filename);
    elseif regexp(filename,'.pdb') > 1
        disp('Found .pdb file');
        atom = import_atom_pdb(filename); % Does the pdb come with occupancy and B-factor info?
    end
    assignin('caller','atom_xrd',atom);
    assignin('caller','Box_dim_xrd',Box_dim)
else
    atom=varargin{1};
    if nargin>1
        Box_dim=varargin{2};
    end
end

% disp('Choose between these representations:')
% disp('ballstick licorice smallvdw halfvdw vdw contour crystal ionic polyhedra lines labels charge index')

if nargin>2
    style = char(varargin{3}); %'ballstick','licorice','halfvdw','vdw'
else
    style = 'licorice';
end

if ~ismember(style,{'ballstick' 'small' 'smallvdw' 'licorice' 'halfvdw' 'vdw' 'contour' 'crystal' 'ionic' 'lines' 'labels' 'charge' 'index' 'poly' 'polyhedra' 'filled'})
    style = 'licorice';
end

bond_radii = 0.12; % bond radii
resolution = 30;  % higher looks better, takes more time
element=element_atom(atom);
XYZ_labels=[element.type]';
nAtoms = size(XYZ_labels,1);


radii = 1/5*abs(radius_vdw(XYZ_labels));
color =  1*element_color(XYZ_labels);

if nargin>4
    alpha=1-varargin{5};
else
    alpha=1; % Transperacy
end

rmaxlong=2.25
distance_factor=0.6;
if nargin>5
    if numel(varargin{6})>0
        rmaxlong=varargin{6}; % Dummy value
        if numel(rmaxlong)>1
            distance_factor=rmaxlong(2);
            rmaxlong=rmaxlong(1);
        end
    end
end

if nargin>6
    trans_vec=varargin{7};
    if numel(trans_vec)==3
        atom=translate_atom(atom,trans_vec(1:3));
        if numel(trans_vec)==4
            atom=wrap_atom(atom,Box_dim);
        end
    end
end

if nargin>7
    color=varargin{8};
    color=repmat(color,nAtoms,1);
end

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];

water_ind=find(ismember(XYZ_labels,{'Ow' 'OW' 'Hw' 'HW' 'HW1' 'HW2'}));
radii(water_ind)=bond_radii;%0.5*radii(water_ind);

assignin('caller','radii',radii)
assignin('caller','color',color)

if nargin>1

    Box_dim=varargin{2};

    if numel(Box_dim)>0

        if numel(Box_dim)==1
            Box_dim(1)=Box_dim(1);
            Box_dim(2)=Box_dim(1);
            Box_dim(3)=Box_dim(1);
        end

        if size(atom,2)>39 && size(atom,2) < 10000 && strcmp(style,'ballstick') || strcmp(style,'licorice')
            disp('Scanning intramolecular bonds, neglecting the PBC')
            atom = bond_atom(atom,1.1*Box_dim,rmaxlong,distance_factor); % the factor 10% makes sure there are no bonds over the pbc!
            % disp('Scanning intramolecular bonds, across the PBC')
            % atom = bond_atom(atom,1.0*Box_dim,rmaxlong,distance_factor); % the factor 10% makes sure there are no bonds over the pbc!
        else
            atom = bond_atom(atom,Box_dim,rmaxlong,distance_factor);

            if max([atom.molid])<max([atom.index])
                i=1;
                while numel(Bond_index)==0 && i < 10
                    rmaxlong=2.2+i/4;
                    distance_factor=6+i/10;
                    atom = bond_atom(atom,Box_dim,rmaxlong,distance_factor);
                    i=i+1;
                end
                disp('Used this rmaxlong and distance_factor:')
                rmaxlong
                distance_factor
            end
        end


    end
    % Sets plot limits for the data
    xlo = floor(min([-5 min([atom.x])-max(radii)])); xhi = ceil(max([max([atom.x])+max(radii) Box_dim(1)])/5)*5;
    ylo = floor(min([-5 min([atom.y])-max(radii)])); yhi = ceil(max([max([atom.y])+max(radii) Box_dim(2)])/5)*5;
    zlo = floor(min([-5 min([atom.z])-max(radii)])); zhi = ceil(max([max([atom.z])+max(radii) Box_dim(3)])/5)*5;
else
    xlo = floor(min([-5 min([atom.x])-max(radii)])); xhi = ceil(max(max([atom.x]))/5)*5;
    ylo = floor(min([-5 min([atom.y])-max(radii)])); yhi = ceil(max(max([atom.y]))/5)*5;
    zlo = floor(min([-5 min([atom.z])-max(radii)])); zhi = ceil(max(max([atom.z]))/5)*5;
end

if xhi < 5
    xhi=5;
end

if yhi < 5
    yhi=5;
end

if zhi < 5
    zhi=5;
end

% xhi=35
% yhi=40
% zhi=45
hold on;
cameratoolbar
rotate3d on;
camlight(220,210,'infinite');
set(gcf,'Visible','on','Color',[1 1 1]);
set(gca,'Color',[1 1 1],'PlotBoxAspectRatio',[(xhi-xlo)/(zhi-zlo) (yhi-ylo)/(zhi-zlo) (zhi-zlo)/(zhi-zlo)],'FontSize',24); %

fig=gcf;
fig.Color = [1 1 1];
ax = fig.CurrentAxes;
ax.XLim = [xlo xhi];
ax.YLim = [ylo yhi];
ax.ZLim = [zlo zhi];

xlabel('X [�]'); ylabel('Y [�]'); zlabel('Z [�]');
view([0,0]);
toc
tic
if strncmpi(style,'licorice',4) || strncmpi(style,'ballstick',4)
    disp('Drawing the atoms')

    [rx,ry,rz] = sphere(resolution);
    r_all=[];
    for i = 1:size(XYZ_data,1)

        switch style

            case 'licorice'
                r_temp = bond_radii;
            case 'ballstick'
                r_temp = radii(i);

        end
        r_all=[r_all; r_temp];
        color_temp = color(i,:);

        surface(XYZ_data(i,1) + r_temp*rx,XYZ_data(i,2) + r_temp*ry, ...
            XYZ_data(i,3) + r_temp*rz,'FaceColor',color_temp, ...
            'EdgeColor','none','FaceLighting','gouraud','FaceAlpha',alpha,...
            'AmbientStrength',.6,'DiffuseStrength',.3,'SpecularStrength',0);


        if mod(i,1000)==1
            if i > 1
                i-1
                drawnow limitrate
            end
        end
    end
end
toc
tic
if exist('Bond_index','var') && numel(Bond_index)>0 && ismember(style,{'ballstick' 'licorice'})
    rdist = Bond_index(:,3);
    % draw cylinders for each bond
    disp('Drawing the bonds')
    for i = 1:size(Bond_index,1) % draw sticks for all bounds
        r1 = XYZ_data(Bond_index(i,1),:); % coordinates atom 1
        r2 = XYZ_data(Bond_index(i,2),:); % coordinates atom 2

        % bond angles in spherical coordinates
        v = (r2-r1)/norm(r2-r1);
        phi = atan2d(v(2),v(1));
        theta = -asind(v(3));

        % bond distance minus sphere radii
        bd = rdist(i) - radii(Bond_index(i,1)) - radii(Bond_index(i,2));
        cyl2 = radii(Bond_index(i,1)) + bd/2; % length half bond cylinder
        cyl1 = rdist(i); % length full bond cylinder

        % get colors of both atoms
        color_temp1 = color(Bond_index(i,2),:);
        color_temp2 = color(Bond_index(i,1),:);

        % prototype cylinders for bond
        [z,y,x] = cylinder(bond_radii,resolution/2); % full bond cylinder
        x(2,:) = x(2,:) * cyl1; % adjust length
        [z2,y2,x2] = cylinder(bond_radii*1.01,resolution/2); % half bond cylinder, thicker
        x2(2,:) = x2(2,:) * cyl2; % adjust length

        % rotate cylinders to match bond vector v
        for kk = 1:numel(x)
            vr = [x(kk); y(kk); z(kk);];
            vr = rotz(phi)*roty(theta)*vr;
            x(kk) = vr(1);
            y(kk) = vr(2);
            z(kk) = vr(3);

            vr = [x2(kk); y2(kk); z2(kk);];
            vr = rotz(phi)*roty(theta)*vr;
            x2(kk) = vr(1);
            y2(kk) = vr(2);
            z2(kk) = vr(3);
        end

        % full bond color 1
        surface(r1(1) + x,r1(2) + y,r1(3) + z,...
            'FaceColor',color_temp1,...
            'EdgeColor','none','FaceLighting','gouraud','FaceAlpha',alpha,...
            'AmbientStrength',.6,'DiffuseStrength',.1,'SpecularStrength',0);
        %'EdgeColor','none',...
        %'FaceLighting','gouraud','FaceAlpha',alpha)

        % half bond color 2
        surface(r1(1) + x2,r1(2) + y2,r1(3) + z2,...
            'FaceColor',color_temp2,...
            'EdgeColor','none','FaceLighting','gouraud','FaceAlpha',alpha,...
            'AmbientStrength',.6,'DiffuseStrength',.3,'SpecularStrength',0);

        if mod(i,100)==1
            if i > 1
                i-1
                drawnow limitrate
            end
        end

    end

end
toc

if nargin>3
    if varargin{4}==0
        disp('Will not draw any box!')
    else
        try
            Simbox = draw_box_atom(Box_dim,[0 0 0.8],2);
        catch
            disp('Could draw the box!')
        end
    end
end

hold off;

end

function rotmat = roty(beta)
% rotate in the direction of z->x, counter-clockwise
rotmat = [cosd(beta) 0 sind(beta); 0 1 0; -sind(beta) 0 cosd(beta)];
end

function rotmat = rotz(gamma)
% rotate in the direction of x->y, counter-clockwise
rotmat = [cosd(gamma) -sind(gamma) 0; sind(gamma) cosd(gamma) 0; 0 0 1];
end
