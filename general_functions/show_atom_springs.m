%% show_atom_springs.m
% * This function draws the atom struct in 3D like show_atom(), but draws helical
% * (coil/shock) springs instead of the usual bond cylinders, to illustrate the
% * harmonic interactions of a force field. See show_spring() for the spring itself.
% *
% * Two kinds of springs can be drawn, set with the 'spring' option:
% *   'bonds'   one spring per pair in Bond_index, i.e. the harmonic bonds
% *   'angles'  one spring per triplet in Angle_index, drawn between atom 1 and
% *             atom 3 of the triplet, i.e. across the angle, with the central
% *             atom 2 left untouched
% *   'both'    both of the above
% *   'none'    no springs, just the atoms
% *
% * Instead of letting bond_atom() find them, an explicit list can be handed over
% * with 'spring_index' as an nx2 matrix of pairs or an nx3 matrix of triplets.
% *
% * The spring color, the thickness of the entire spring ('springWidth', the coil
% * diameter) and the thickness of its wire ('wireWidth' for the 3D tube style, or
% * 'linewidth' for the fast line style) are all set by the user.
% *
% * For a triplet the spring spans atom 1 to atom 3 by default, but it can be pulled
% * in towards the central atom so that it sits along and between the two bonds:
% *   'angleSpan'  a fraction of each bond, 1 = out at the atoms, 0.4 = nearer the vertex
% *   'angleDist'  a fixed distance from the vertex instead, in Angstrom
% *
% * The angle itself can be marked with a surface between atoms 1, 2 and 3, set with
% * 'angleFace':
% *   'triangle'  the flat triangle through the three atoms
% *   'sector'    a circular sector seated at the vertex, swept from one bond to the
% *               other. Its size is set either with 'sectorRadius', an absolute
% *               radius in Angstrom, or with 'sectorSpan', a fraction of the shorter
% *               of the two bonds (default 0.5). sectorRadius wins when both are set,
% *               and a radius longer than a bond simply reaches past that atom.
% * Its look is set with 'faceColor' (default: the spring color), 'faceAlpha' and
% * 'faceEdge'.
% *
% * The atom sizes follow the representation style, but can be overridden with
% * 'radiiScale', a plain multiplier, or with 'radii', an explicit scalar or one
% * value per atom. Note that the radii variable that this function writes back to
% * the calling workspace is an OUTPUT only, so setting it there has no effect, use
% * these options instead.
% *
% * Set 'atoms' to false to draw the springs on their own, without any atoms. The
% * springs then span the full distance, since there are no spheres to clear.
% *
% * Springs stretching over the periodic boundary are dropped by default, since they
% * otherwise shoot straight across the whole box. This is set with 'pbc':
% *   'skip'  leave them out, the default
% *   'wrap'  draw them towards the nearest periodic image instead
% *   'draw'  draw them as they are, from atom to atom
% * The test compares the direct separation with the minimum image, so it is exact
% * and also correct for triclinic cells.
% *
% * Note that a full clay or water system has far more angles than bonds, so use
% * 'spring_index', 'spring_types' or 'maxsprings' to keep the number sensible.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # show_atom_springs(atom,Box_dim)
% # show_atom_springs(atom,Box_dim,'licorice') % Any show_atom representation style
% # show_atom_springs(atom,Box_dim,'spring','angles')
% # show_atom_springs(atom,Box_dim,'spring','both','springColor',[1 0 0])
% # show_atom_springs(atom,Box_dim,'springWidth',1.0,'wireWidth',0.15,'coils',8)
% # show_atom_springs(atom,Box_dim,'springStyle','line','linewidth',2)
% # show_atom_springs(atom,Box_dim,'spring_index',[1 2; 3 4]) % Explicit pairs
% # show_atom_springs(atom,Box_dim,'spring_index',[1 2 3]) % Explicit triplet
% # show_atom_springs(atom,Box_dim,'spring','angles','spring_types',{'Sit'})
% # show_atom_springs(atom,Box_dim,'radiiScale',0.25) % Quarter sized atoms
% # show_atom_springs(atom,Box_dim,'radii',0.3) % Every atom drawn with radius 0.3 A
% # show_atom_springs(atom,Box_dim,'atoms',false) % Springs only, no atoms drawn
% # show_atom_springs(atom,Box_dim,'pbc','wrap') % Draw PBC springs to the nearest image
% # show_atom_springs(atom,Box_dim,'pbc','draw') % Keep the springs crossing the box
% # show_atom_springs(atom,Box_dim,'spring','angles','angleSpan',0.45) % Nearer the vertex
% # show_atom_springs(atom,Box_dim,'spring','angles','angleDist',1.0) % 1 A from the vertex
% # show_atom_springs(atom,Box_dim,'spring','angles','angleFace','triangle')
% # show_atom_springs(atom,Box_dim,'spring','angles','angleFace','sector','faceAlpha',0.4)
% # show_atom_springs(atom,Box_dim,'spring','angles','angleFace','sector','sectorRadius',1.2)
% # show_atom_springs(atom,Box_dim,'spring','angles','angleFace','sector','sectorSpan',0.8)
%

function show_atom_springs(atom,Box_dim,varargin)
%%

if nargin < 2, Box_dim = []; end

% --- an optional positional representation style, as in show_atom ---
known_styles = {'ballstick','small','smallvdw','licorice','halfvdw','vdw','contour', ...
                'crystal','ionic','lines','labels','charge','index','filled'};
style = 'ballstick';
if numel(varargin) > 0 && (ischar(varargin{1}) || isstring(varargin{1})) ...
        && any(strcmpi(varargin{1}, known_styles))
    style = char(varargin{1});
    varargin(1) = [];
end

% --- options (name/value) ---
p = struct('style',style, ...
           'atoms',true, ...                % false = draw the springs only, no atoms
           'pbc','skip', ...                % skip | draw | wrap, for springs over the PBC
           'spring','bonds', ...            % bonds | angles | both | none
           'spring_index',[], ...           % explicit nx2 pairs or nx3 triplets
           'spring_types',{{}}, ...         % restrict to springs touching these types
           'springColor',[0.25 0.25 0.25], ...
           'springWidth',0.7, ...           % thickness of the entire spring
           'wireWidth',0.10, ...            % thickness of the spring wire
           'linewidth',1.5, ...             % wire thickness for springStyle 'line'
           'springStyle','tube', ...        % tube | line
           'coils',6, ...
           'springGap',[], ...              % [] = shrink by the atom radii
           'springAlpha',1, ...
           'angleSpan',1, ...               % triplets: fraction along the two bonds
           'angleDist',[], ...              % triplets: fixed distance from the vertex
           'angleFace','none', ...          % none | triangle | sector
           'faceColor',[], ...              % [] = the spring color
           'faceAlpha',0.25, ...
           'faceEdge','none', ...           % edge color of the face
           'sectorRadius',[], ...           % absolute sector radius in Angstrom
           'sectorSpan',0.5, ...            % or a fraction of the shorter bond
           'sectorPoints',24, ...
           'maxsprings',5000, ...           % safety cap
           'drawbonds',false, ...           % also draw the usual bond cylinders
           'box',1, ...
           'transparency',0, ...
           'radii',[], ...                  % explicit atom radii, scalar or per atom
           'radiiScale',1, ...              % multiplier on whatever radii the style gives
           'rmaxlong',2.45, ...
           'distance_factor',0.6, ...
           'trans_vec',[], ...
           'color',[]);
for k=1:2:numel(varargin)
    fn = fieldnames(p);
    hit = find(strcmpi(fn,varargin{k}),1);
    if isempty(hit)
        error('show_atom_springs:UnknownOption','Unknown option "%s".', num2str(varargin{k}));
    end
    p.(fn{hit}) = varargin{k+1};
end
style = char(p.style);
if ~ismember(style,{'ballstick' 'small' 'smallvdw' 'licorice' 'halfvdw' 'vdw' ...
                    'contour' 'crystal' 'ionic' 'lines' 'labels' 'charge' 'index' 'filled'})
    style = 'ballstick';
end

if ~isempty(p.trans_vec) && numel(p.trans_vec)>=3
    atom = translate_atom(atom, p.trans_vec(1:3));
end

bond_radii = 0.12;
resolution = 30;
XYZ_labels = [atom.type]';
nAtoms     = size(XYZ_labels,1);

if strncmpi(style,'crystal',3) || strncmpi(style,'filled',4)
    radii = 1/5*abs(radius_crystal(XYZ_labels));
elseif strncmpi(style,'ionic',3)
    radii = 1/5*abs(radius_ion(XYZ_labels));
else
    radii = 1/5*abs(radius_vdw(XYZ_labels));
end
radii = radii(:);

color = 1*element_color(XYZ_labels);
dark  = ismember(XYZ_labels,{'Alt' 'Fet'});     % darken the tetrahedral substitutions
color(dark,:) = 0.3*color(dark,:);              % note: row-wise, unlike show_atom
if ~isempty(p.color)
    color = repmat(p.color,nAtoms,1);
end
alpha = 1 - p.transparency;

XYZ_data = [[atom.x]' [atom.y]' [atom.z]'];
water_ind = find(ismember(XYZ_labels,{'Ow' 'OW' 'Hw' 'HW' 'HW1' 'HW2'}));
radii(water_ind) = bond_radii;

% --- user control over the atom radii ---
% Applied last, so an explicit radii also wins over the water override above. The
% scale is carried by the vector itself, so the axis limits and the default
% springGap follow along with it.
if ~isempty(p.radii)
    if isscalar(p.radii)
        radii = p.radii*ones(nAtoms,1);
    elseif numel(p.radii) == nAtoms
        radii = p.radii(:);
    else
        error('show_atom_springs:BadRadii', ...
              'radii must be a scalar or have one value per atom (%d).', nAtoms);
    end
end
radii = p.radiiScale * radii;

%% --- find the bonds and angles, unless an explicit list was given ---
Bond_index = []; Angle_index = [];
need_scan = isempty(p.spring_index) && ~strcmpi(p.spring,'none');
if (need_scan || p.drawbonds) && ~isempty(Box_dim)
    if numel(Box_dim)==1, Box_dim = [Box_dim Box_dim Box_dim]; end
    disp('Scanning intramolecular bonds and angles, neglecting the PBC')
    evalc('atom = bond_atom(atom,1.1*Box_dim,p.rmaxlong,p.distance_factor);');
end

%% --- assemble the spring list ---
springs = [];   % nx2 (pairs) or nx3 (triplets, spring drawn between col 1 and col 3)
if ~isempty(p.spring_index)
    springs = p.spring_index;
elseif ~strcmpi(p.spring,'none')
    if any(strcmpi(p.spring,{'bonds','both'})) && ~isempty(Bond_index)
        springs = [springs; Bond_index(:,1:2)];
    end
    if any(strcmpi(p.spring,{'angles','both'})) && ~isempty(Angle_index)
        if isempty(springs)
            springs = Angle_index(:,1:3);
        else
            springs = [springs, nan(size(springs,1),1); Angle_index(:,1:3)];
        end
    end
end

if ~isempty(springs)
    % keep only springs touching one of spring_types, if asked for
    if ~isempty(p.spring_types)
        keep = false(size(springs,1),1);
        for r = 1:size(springs,1)
            idx = springs(r,~isnan(springs(r,:)));
            keep(r) = any(ismember(XYZ_labels(idx), p.spring_types));
        end
        springs = springs(keep,:);
    end
    if size(springs,1) > p.maxsprings
        fprintf('Capping the springs at maxsprings = %d of %d\n', p.maxsprings, size(springs,1));
        springs = springs(1:p.maxsprings,:);
    end
end

%% --- set up the figure, as in show_atom ---
if ~isempty(Box_dim)
    xlo = floor(min([-5 min([atom.x])-max(radii)])); xhi = ceil(max([max([atom.x])+max(radii) Box_dim(1)])/5)*5;
    ylo = floor(min([-5 min([atom.y])-max(radii)])); yhi = ceil(max([max([atom.y])+max(radii) Box_dim(2)])/5)*5;
    zlo = floor(min([-5 min([atom.z])-max(radii)])); zhi = ceil(max([max([atom.z])+max(radii) Box_dim(3)])/5)*5;
else
    xlo = floor(min([-5 min([atom.x])-max(radii)])); xhi = ceil(max(max([atom.x]))/5)*5;
    ylo = floor(min([-5 min([atom.y])-max(radii)])); yhi = ceil(max(max([atom.y]))/5)*5;
    zlo = floor(min([-5 min([atom.z])-max(radii)])); zhi = ceil(max(max([atom.z]))/5)*5;
end
xhi = max(xhi,5); yhi = max(yhi,5); zhi = max(zhi,5);

hold on; rotate3d on;
% camlight(220,210,'infinite');
set(gcf,'Visible','on','Color',[1 1 1]);
set(gca,'Color',[1 1 1], ...
    'PlotBoxAspectRatio',[(xhi-xlo)/(zhi-zlo) (yhi-ylo)/(zhi-zlo) 1],'FontSize',24);
ax = gca; ax.XLim=[xlo xhi]; ax.YLim=[ylo yhi]; ax.ZLim=[zlo zhi];
xlabel('X [Å]'); ylabel('Y [Å]'); zlabel('Z [Å]');
view([0,0]);

%% --- draw the atoms, unless only the springs were asked for ---
if ~p.atoms
    % nothing drawn here
elseif strncmpi(style,'labels',5) || strncmpi(style,'charge',5) || strncmpi(style,'index',5)
    if strncmpi(style,'labels',5)
        text(XYZ_data(:,1)-.2,XYZ_data(:,2)-.2,XYZ_data(:,3)+.1,[atom.type],'FontSize',12);
    elseif strncmpi(style,'charge',5)
        try
            text(XYZ_data(:,1)-.4,XYZ_data(:,2)-.4,XYZ_data(:,3)+.1, ...
                 strsplit(num2str([atom.charge])),'FontSize',16);
        catch
            disp('Did not find any charge...')
        end
    else
        text(XYZ_data(:,1)-.2,XYZ_data(:,2)-.2,XYZ_data(:,3)+.1, ...
             strsplit(num2str([atom.index])),'FontSize',16);
    end
elseif ~strncmpi(style,'lines',4)
    disp('Drawing the atoms')
    [rx,ry,rz] = sphere(resolution);
    for i = 1:nAtoms
        switch style
            case {'ballstick','small','smallvdw'}, r_temp = p.radiiScale*radii(i);
            case 'licorice',                       r_temp = p.radiiScale*bond_radii;
            case {'halfvdw','contour'},            r_temp = 5/2*p.radiiScale*radii(i);
            otherwise,                             r_temp = 5*p.radiiScale*radii(i);
        end
        surface(XYZ_data(i,1) + r_temp*rx, XYZ_data(i,2) + r_temp*ry, ...
                XYZ_data(i,3) + r_temp*rz,'FaceColor',color(i,:), ...
                'EdgeColor','none','FaceLighting','gouraud','FaceAlpha',alpha, ...
                'AmbientStrength',.6,'DiffuseStrength',.3,'SpecularStrength',0);
        if mod(i,1000)==1 && i>1, drawnow limitrate; end
    end
end

%% --- draw the springs ---
if ~isempty(springs)
    if ~isempty(p.spring_index), what = 'explicit spring_index'; else, what = p.spring; end
    fprintf('Drawing %d springs (%s)\n', size(springs,1), what);

    % --- box vectors, for the minimum-image test below ---
    havebox = ~isempty(Box_dim) && numel(Box_dim) >= 3;
    if havebox
        lx=Box_dim(1); ly=Box_dim(2); lz=Box_dim(3);
        if numel(Box_dim) >= 9
            xy=Box_dim(6); xz=Box_dim(8); yz=Box_dim(9);
        else
            xy=0; xz=0; yz=0;
        end
    end

    nskipped = 0;
    for i = 1:size(springs,1)
        row = springs(i,:);
        row = row(~isnan(row));
        istriplet = numel(row) >= 3;
        if istriplet
            i1 = row(1); iv = row(2); i3 = row(3);
            if any([i1 iv i3]<1) || any([i1 iv i3]>nAtoms), continue; end
        elseif numel(row) == 2
            i1 = row(1); i3 = row(2); iv = [];
            if any([i1 i3]<1) || any([i1 i3]>nAtoms), continue; end
        else
            continue
        end

        % --- geometry, taken relative to the vertex for a triplet ---
        % Anchoring on the vertex means the minimum image is applied to both bonds,
        % which is what keeps a wrapped triplet together and its angle correct.
        if istriplet
            rv = XYZ_data(iv,:);
            d1 = XYZ_data(i1,:) - rv;
            d3 = XYZ_data(i3,:) - rv;
            crossed = false;
            if havebox
                m1 = mic_(d1, lx,ly,lz,xy,xz,yz);
                m3 = mic_(d3, lx,ly,lz,xy,xz,yz);
                crossed = norm(d1-m1) > 1e-6 || norm(d3-m3) > 1e-6;
            end
        else
            rv = XYZ_data(i1,:);
            d1 = [0 0 0];
            d3 = XYZ_data(i3,:) - rv;
            crossed = false;
            if havebox
                m1 = d1;
                m3 = mic_(d3, lx,ly,lz,xy,xz,yz);
                crossed = norm(d3-m3) > 1e-6;
            end
        end

        % --- interactions stretching over the periodic boundary ---
        if crossed && ~strncmpi(p.pbc,'draw',4)
            if strncmpi(p.pbc,'wrap',4)
                d1 = m1; d3 = m3;          % use the nearest images instead
            else
                nskipped = nskipped + 1;
                continue                   % 'skip', the default
            end
        end
        P1 = rv + d1;  P3 = rv + d3;       % the atom positions actually drawn

        % --- where the spring ends sit ---
        % For a triplet the ends may be pulled in along the two bonds towards the
        % vertex, either by a fraction of each bond or at a fixed distance from it.
        if istriplet
            if ~isempty(p.angleDist)
                f1 = min(1, p.angleDist/max(eps,norm(d1)));
                f3 = min(1, p.angleDist/max(eps,norm(d3)));
            else
                f1 = p.angleSpan; f3 = p.angleSpan;
            end
            r1 = rv + f1*d1;  r2 = rv + f3*d3;
            pulled_in = (f1 < 1) || (f3 < 1);
        else
            r1 = P1; r2 = P3; pulled_in = false;
        end

        % --- the angle face, drawn between atoms 1, 2 and 3 ---
        if istriplet && ~strncmpi(p.angleFace,'none',4)
            draw_face_(rv, d1, d3, p);
        end

        if isempty(p.springGap)
            if p.atoms && ~pulled_in
                gap = min([radii(i1) radii(i3) 0.45*norm(r2-r1)]);  % stop at the atom surfaces
            else
                gap = 0;     % no atoms to clear, or the spring already starts inside them
            end
        else
            gap = p.springGap;
        end

        show_spring(r1,r2,'color',p.springColor,'springWidth',p.springWidth, ...
                    'wireWidth',p.wireWidth,'linewidth',p.linewidth, ...
                    'coils',p.coils,'style',p.springStyle,'gap',gap, ...
                    'facealpha',p.springAlpha);
        if mod(i,200)==1 && i>1, drawnow limitrate; end
    end
    if nskipped > 0
        fprintf('Skipped %d spring(s) stretching over the PBC\n', nskipped);
    end
end

%% --- optionally also the usual bond cylinders ---
if p.drawbonds && ~isempty(Bond_index)
    disp('Drawing the bonds')
    for i = 1:size(Bond_index,1)
        r1 = XYZ_data(Bond_index(i,1),:);
        r2 = XYZ_data(Bond_index(i,2),:);
        v  = (r2-r1)/norm(r2-r1);
        phi = atan2d(v(2),v(1)); theta = -asind(v(3));
        [z,y,x] = cylinder(bond_radii,resolution/2);
        x(2,:) = x(2,:)*Bond_index(i,3);
        for kk = 1:numel(x)
            vr = rotz_(phi)*roty_(theta)*[x(kk); y(kk); z(kk)];
            x(kk)=vr(1); y(kk)=vr(2); z(kk)=vr(3);
        end
        surface(r1(1)+x, r1(2)+y, r1(3)+z,'FaceColor',color(Bond_index(i,2),:), ...
                'EdgeColor','none','FaceLighting','gouraud','FaceAlpha',alpha, ...
                'AmbientStrength',.6,'DiffuseStrength',.1,'SpecularStrength',0);
    end
end

%% --- the simulation box ---
if p.box ~= 0 && ~isempty(Box_dim)
    try
        Simbox = draw_box_atom(Box_dim,[0 0 0.8],2); %#ok<NASGU>
    catch
        disp('Could not draw the box!')
    end
end

assignin('caller','springs',springs);
assignin('caller','radii',radii);
assignin('caller','color',color);
hold off;

end  % main function


% ===== local helpers =====
function draw_face_(rv, d1, d3, p)
% Mark the angle at the vertex rv spanned by the two bond vectors d1 and d3,
% either as the flat triangle through atoms 1, 2 and 3, or as a circular sector
% seated at the vertex and swept from one bond to the other.
if isempty(p.faceColor), fc = p.springColor; else, fc = p.faceColor; end
n1 = norm(d1); n3 = norm(d3);
if n1 < eps || n3 < eps, return; end

if strncmpi(p.angleFace,'triangle',3)
    V = [rv; rv+d1; rv+d3];
    patch('Faces',[1 2 3],'Vertices',V,'FaceColor',fc,'FaceAlpha',p.faceAlpha, ...
          'EdgeColor',p.faceEdge,'FaceLighting','none');
    return
end

% --- sector: fan of triangles from the vertex, along the arc between the bonds ---
u1 = d1/n1; u3 = d3/n3;
nrm = cross(u1,u3);
if norm(nrm) < 1e-10, return; end          % collinear, no angle to show
nrm = nrm/norm(nrm);
th  = atan2(norm(cross(u1,u3)), dot(u1,u3));
% sectorRadius is an absolute radius, sectorSpan a fraction of the shorter bond
if ~isempty(p.sectorRadius)
    R = p.sectorRadius;
else
    R = p.sectorSpan * min(n1,n3);
end
m   = max(3,p.sectorPoints);
tt  = linspace(0,th,m);
arc = zeros(m,3);
for k = 1:m
    arc(k,:) = rv + R*(u1*cos(tt(k)) + cross(nrm,u1)*sin(tt(k)) + ...
                       nrm*dot(nrm,u1)*(1-cos(tt(k))));   % Rodrigues about nrm
end
V = [rv; arc];
F = [ones(m-1,1), (2:m)', (3:m+1)'];
patch('Faces',F,'Vertices',V,'FaceColor',fc,'FaceAlpha',p.faceAlpha, ...
      'EdgeColor',p.faceEdge,'FaceLighting','none');
end

function d = mic_(d, lx,ly,lz,xy,xz,yz)
% Triclinic minimum-image (GROMACS Box_dim convention), matching bond_angle_type.
rx=d(:,1); ry=d(:,2); rz=d(:,3);
gt=rz>lz/2;  lt=rz<-lz/2;
rz(gt)=rz(gt)-lz; rz(lt)=rz(lt)+lz;
rx(gt)=rx(gt)-xz; rx(lt)=rx(lt)+xz;
ry(gt)=ry(gt)-yz; ry(lt)=ry(lt)+yz;
gt=ry>ly/2;  lt=ry<-ly/2;
ry(gt)=ry(gt)-ly; ry(lt)=ry(lt)+ly;
rx(gt)=rx(gt)-xy; rx(lt)=rx(lt)+xy;
gt=rx>lx/2;  lt=rx<-lx/2;
rx(gt)=rx(gt)-lx; rx(lt)=rx(lt)+lx;
d=[rx ry rz];
end

function rotmat = roty_(beta)
rotmat = [cosd(beta) 0 sind(beta); 0 1 0; -sind(beta) 0 cosd(beta)];
end

function rotmat = rotz_(gamma)
rotmat = [cosd(gamma) -sind(gamma) 0; sind(gamma) cosd(gamma) 0; 0 0 1];
end
