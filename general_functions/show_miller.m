%%  show_miller.m
% * This function draws up the h,k,l Miller plane within the Box_dim/Cell

%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% #  show_miller(h,k,l)
% #  show_miller(h,k,l,Box_dim)             % 3, 6, or 9 elements
% #  show_miller(h,k,l,Box_dim,Color)       % Color: 'b' or [r g b]
% #  show_miller(h,k,l,Box_dim,Color,FaceAlpha)   % FaceAlpha in [0..1]
%
% Name–Value options (pass after the positionals above):
%   'SinglePlane'     true/false (default false)
%   'PlaneLevel'      integer n or 'auto' (default 'auto') for hx+ky+lz = n
%   'ExcludeBoundary' true/false (default false): suppress edge/vertex-only contacts
%   'Box_dim'         override positional Box_dim if desired (3/6/9)
%   'Color'           override positional Color
%   'FaceAlpha'       override positional FaceAlpha
%
% Notes
% • Handles ANY sign combination (e.g., -1 -1 1, -1 -1 -1).
% • Determines all integer n that intersect the unit fractional cell by scanning
%   s=hx+ky+lz over the 8 corners; optionally draws a single plane.
% • Box_dim supports [a b c], [a b c α β γ], or [lx ly lz 0 0 xy 0 xz yz].
%
% Returns: hPatches (column vector of patch handles).

function hPatches = show_miller(h,k,l,Box_dim,color,faceAlpha, varargin)

    % -------- Defaults for positionals --------
    if nargin<4 || isempty(Box_dim),   Box_dim   = [1 1 1]; end
    if nargin<5 || isempty(color),     color     = [0 0 1]; end
    if nargin<6 || isempty(faceAlpha), faceAlpha = 0.5;     end

    % -------- Manual NV parsing (robust; no inputParser) --------
    opts.SinglePlane     = false;
    opts.PlaneLevel      = 'auto';
    opts.ExcludeBoundary = false;

    i = 1;
    while i <= numel(varargin)
        name = varargin{i};
        if ~(ischar(name) || (isstring(name)&&isscalar(name)))
            error('show_miller:InvalidNameValue', ...
                  'Expected name–value pairs after positional arguments.');
        end
        if i==numel(varargin)
            error('show_miller:MissingValue','Missing value for "%s".', name);
        end
        val = varargin{i+1};
        switch lower(char(name))
            case 'singleplane'
                opts.SinglePlane = logical(val);
            case 'planelevel'
                if ischar(val) || isstring(val)
                    opts.PlaneLevel = char(val);
                else
                    opts.PlaneLevel = round(val);
                end
            case 'excludeboundary'
                opts.ExcludeBoundary = logical(val);
            case 'box_dim'
                Box_dim = val;
            case 'color'
                color = val;
            case 'facealpha'
                faceAlpha = val;
            otherwise
                error('show_miller:UnknownOption','Unknown option "%s".', name);
        end
        i = i + 2;
    end

    % -------- Validate & normalize inputs --------
    h = round(h); k = round(k); l = round(l);
    assert(any([h k l]~=0), 'At least one of h,k,l must be non-zero.');

    if ischar(color) || (isstring(color)&&isscalar(color))
        color = local_colorSpecToRGB(char(color));
    end

    [FromFrac,~,~,~,~,~,~] = local_from_Boxdim(Box_dim);


    % -------- Fractional cube & edges --------
    Vfrac = [0 0 0; 1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1];
    E = [1 2;3 4;5 6;7 8; 1 3;2 4;5 7;6 8; 1 5;2 6;3 7;4 8]; % 12 edges

    % Integer n range that can intersect the cell
    svals = Vfrac*[h; k; l];
    nmin  = ceil(min(svals));
    nmax  = floor(max(svals));

    hold_state = ishold; hold on;

    % Draw wireframe of the (possibly skew) cell
    Vcart = (FromFrac * Vfrac')';
    local_draw_cell_wireframe(Vcart, E);

    hPatches = gobjects(0,1);

    % -------- Single plane mode --------
    if opts.SinglePlane
        if ischar(opts.PlaneLevel) || isstring(opts.PlaneLevel)
            % 'auto': pick usable n nearest to 1
            nlist = nmin:nmax;
            if isempty(nlist), if ~hold_state, hold off; end, return; end
            [n_sel,Pc_sel] = local_pick_single_level(nlist, E, Vfrac, FromFrac, [h k l], opts.ExcludeBoundary);
            if isempty(n_sel), if ~hold_state, hold off; end, return; end
            hPatches(1,1) = patch('XData',Pc_sel(:,1),'YData',Pc_sel(:,2),'ZData',Pc_sel(:,3), ...
                                  'FaceColor',color,'FaceAlpha',faceAlpha,'EdgeColor','none');
        else
            n = round(opts.PlaneLevel);
            [Pfrac,Pc,ok] = local_intersect_polygon(E,Vfrac,FromFrac,[h k l],n);
            if ok && local_polygon_ok(Pfrac,opts.ExcludeBoundary)
                hPatches(1,1) = patch('XData',Pc(:,1),'YData',Pc(:,2),'ZData',Pc(:,3), ...
                                      'FaceColor',color,'FaceAlpha',faceAlpha,'EdgeColor','none');
            end
        end
    else
        % -------- All planes mode --------
        for n = nmin:nmax
            [Pfrac,Pc,ok] = local_intersect_polygon(E,Vfrac,FromFrac,[h k l],n);
            if ~ok || ~local_polygon_ok(Pfrac,opts.ExcludeBoundary), continue; end
            hPatches(end+1,1) = patch('XData',Pc(:,1),'YData',Pc(:,2),'ZData',Pc(:,3), ...
                                      'FaceColor',color,'FaceAlpha',faceAlpha,'EdgeColor','none'); %#ok<AGROW>
        end
    end

    % -------- Cosmetics --------
    axis equal; view(3);
    camlight(220,210,'infinite'); rotate3d on;
    set(gcf,'Color',[1 1 1]); set(gca,'Color',[1 1 1],'FontSize',12);
    xlabel('a [Å]'); ylabel('b [Å]'); zlabel('c [Å]');
    axis normal tight equal

    if ~hold_state, hold off; end
end

% ========================= Helpers =========================

function tf = local_polygon_ok(Pfrac,excludeBoundary)
    if size(Pfrac,1) < 3, tf=false; return; end
    if ~excludeBoundary, tf=true; return; end
    area = local_polyarea3d(Pfrac);
    tf = area > 1e-12;
end

function [Pfrac,Pc,ok] = local_intersect_polygon(E,Vfrac,FromFrac,hkl,n)
    h=hkl(1); k=hkl(2); l=hkl(3);
    P = zeros(0,3);
    for e = 1:size(E,1)
        v0 = Vfrac(E(e,1),:); v1 = Vfrac(E(e,2),:);
        dv = v1 - v0;
        denom = h*dv(1) + k*dv(2) + l*dv(3);
        num   = n - (h*v0(1) + k*v0(2) + l*v0(3));
        if abs(denom) < 1e-14, continue; end
        t = num/denom;
        if t >= -1e-12 && t <= 1+1e-12
            p = v0 + t*dv;
            if all(p >= -1e-12 & p <= 1+1e-12)
                P(end+1,:) = min(max(p,0),1); %#ok<AGROW>
            end
        end
    end
    if isempty(P), Pfrac=[]; Pc=[]; ok=false; return; end
    P = uniquetol(P,1e-10,'ByRows',true);
    if size(P,1) < 3, Pfrac=[]; Pc=[]; ok=false; return; end

    % order vertices in-plane
    c0   = mean(P,1);
    nvec = [h k l];
    ref  = [1 0 0];
    if norm(cross(nvec,ref)) < 1e-12, ref = [0 1 0]; end
    if norm(cross(nvec,ref)) < 1e-12, ref = [0 0 1]; end
    u = cross(nvec,ref); nu = norm(u);
    if nu < 1e-14
        [~,~,Vsvd] = svd(bsxfun(@minus,P,c0),'econ');
        u = Vsvd(:,1).'; v = Vsvd(:,2).';
    else
        u = u/nu; v = cross(nvec,u); v = v/norm(v);
    end
    A = [(P - c0)*u' (P - c0)*v'];
    [~,idx] = sort(atan2(A(:,2),A(:,1)));
    Pfrac = P(idx,:);
    Pc = (FromFrac * Pfrac')';  % to Cartesian
    ok = true;
end

function [n_sel,Pc_sel] = local_pick_single_level(nlist, E, Vfrac, FromFrac, hkl, excludeBoundary)
    [~,order] = sort(abs(nlist - 1),'ascend');  % nearest to n=1
    n_sel = []; Pc_sel = [];
    for ii = 1:numel(order)
        n = nlist(order(ii));
        [Pfrac,Pc,ok] = local_intersect_polygon(E,Vfrac,FromFrac,hkl,n);
        if ~ok || ~local_polygon_ok(Pfrac,excludeBoundary), continue; end
        n_sel = n; Pc_sel = Pc; return
    end
end

function area = local_polyarea3d(Pfrac)
    c0 = mean(Pfrac,1); m = size(Pfrac,1); acc = [0 0 0];
    for i=1:m
        a = Pfrac(i,:) - c0;
        b = Pfrac(mod(i,m)+1,:) - c0;
        acc = acc + cross(a,b);
    end
    area = 0.5 * norm(acc);
end

function [FromFrac,a,b,c,alpha_deg,beta_deg,gamma_deg] = local_from_Boxdim(Box_dim)
    if numel(Box_dim)==3
        a=Box_dim(1); b=Box_dim(2); c=Box_dim(3);
        alpha_deg=90; beta_deg=90; gamma_deg=90;
    elseif numel(Box_dim)==6
        a=Box_dim(1); b=Box_dim(2); c=Box_dim(3);
        alpha_deg=Box_dim(4); beta_deg=Box_dim(5); gamma_deg=Box_dim(6);
    elseif numel(Box_dim)==9
        lx=Box_dim(1); ly=Box_dim(2); lz=Box_dim(3);
        xy=Box_dim(6); xz=Box_dim(8); yz=Box_dim(9);
        a=lx; b = hypot(ly,xy); c = sqrt(lz^2 + xz^2 + yz^2);
        alpha_deg = rad2deg(acos((ly*yz+xy*xz)/(b*c)));
        beta_deg  = rad2deg(acos(xz/c));
        gamma_deg = rad2deg(acos(xy/b));
    else
        error('Box_dim must have 3, 6, or 9 elements.');
    end
    ga = deg2rad(gamma_deg); be = deg2rad(beta_deg); al = deg2rad(alpha_deg);
    v = sqrt(1 - cos(al)^2 - cos(be)^2 - cos(ga)^2 + 2*cos(al)*cos(be)*cos(ga));
    FromFrac = [ a, b*cos(ga),  c*cos(be); ...
                 0, b*sin(ga),  c*(cos(al)-cos(be)*cos(ga))/sin(ga); ...
                 0,        0,   c*v/sin(ga) ];
end

function local_draw_cell_wireframe(Vcart, E)
    lw = 2.0; ec = [0 0 1];
    for i = 1:size(E,1)
        p = Vcart(E(i,1),:); q = Vcart(E(i,2),:);
        line([p(1),q(1)],[p(2),q(2)],[p(3),q(3)],'Color',ec,'LineWidth',lw);
    end
end

function rgb = local_colorSpecToRGB(c)
    switch char(c)
        case 'b', rgb=[0 0 1];
        case 'g', rgb=[0 1 0];
        case 'r', rgb=[1 0 0];
        case 'c', rgb=[0 1 1];
        case 'm', rgb=[1 0 1];
        case 'y', rgb=[1 1 0];
        case 'k', rgb=[0 0 0];
        case 'w', rgb=[1 1 1];
        otherwise, error('Unknown color spec: %s', c);
    end
end


