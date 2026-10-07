%% show_spring.m
% * This function draws a 3D helical (coil/shock) spring between two points p1 and p2.
% * It is meant to illustrate a harmonic interaction, either a bond between a pair of
% * atoms or an angle between the two outer atoms of a triplet, see show_atom_springs().
% *
% * The spring is drawn either as a swept tube ('tube', looks best with lighting) or as
% * a simple 3D line ('line', much faster for many springs). The coil amplitude is
% * tapered to zero over the first and last 'lead' fraction, so the spring starts and
% * ends exactly on the axis p1->p2.
% *
% * springWidth sets the thickness of the entire spring (the coil diameter), while
% * wireWidth sets the thickness of the wire itself for the 'tube' style, and linewidth
% * does the same for the 'line' style.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # show_spring([0 0 0],[3 0 0])
% # show_spring([0 0 0],[3 0 0],'color',[1 0 0])
% # show_spring([0 0 0],[3 0 0],'springWidth',0.8,'wireWidth',0.12,'coils',8)
% # show_spring([0 0 0],[3 0 0],'style','line','linewidth',2,'color','r')
% # show_spring([0 0 0],[3 0 0],'gap',0.4) % Shrink both ends, e.g. to clear the atoms
%

function h = show_spring(p1,p2,varargin)
%%

% --- options (name/value) ---
p = struct('color',[0.2 0.2 0.2], ...  % spring color
           'springWidth',0.7, ...      % thickness of the entire spring (coil diameter)
           'wireWidth',0.10, ...       % thickness of the wire itself, 'tube' style
           'linewidth',1.5, ...        % wire thickness in points, 'line' style
           'coils',6, ...              % number of turns
           'style','tube', ...         % 'tube' or 'line'
           'lead',0.12, ...            % straight/tapered fraction at each end
           'resolution',10, ...        % cross-section points, 'tube' style
           'pointsPerCoil',24, ...     % centerline points per turn
           'gap',0, ...                % shorten both ends by this length
           'facealpha',1);
for k=1:2:numel(varargin)
    nm = varargin{k};
    fn = fieldnames(p);
    hit = find(strcmpi(fn,nm),1);
    if isempty(hit)
        error('show_spring:UnknownOption','Unknown option "%s".', num2str(nm));
    end
    p.(fn{hit}) = varargin{k+1};
end

p1 = p1(:)'; p2 = p2(:)';
ax = p2 - p1;
L  = norm(ax);
if L <= 0, h = []; return; end
w  = ax / L;

% --- optionally clear the ends, e.g. so the spring starts at the atom surfaces ---
if p.gap > 0
    if 2*p.gap >= L, h = []; return; end
    p1 = p1 + p.gap*w;
    p2 = p2 - p.gap*w;
    L  = L - 2*p.gap;
end

% --- an orthonormal frame around the axis ---
tmp = [0 0 1];
if abs(dot(tmp,w)) > 0.9, tmp = [1 0 0]; end
u1 = cross(w,tmp); u1 = u1/norm(u1);
u2 = cross(w,u1);

% --- helical centerline, with the amplitude tapered to zero at both ends ---
M   = max(16, round(p.pointsPerCoil * p.coils));
t   = linspace(0,1,M)';
R   = p.springWidth/2;
amp = R*ones(M,1);
ld  = max(eps, p.lead);
il  = t < ld;        amp(il) = R * (t(il)/ld);
ir  = t > 1-ld;      amp(ir) = R * ((1-t(ir))/ld);
ph  = 2*pi*p.coils*t;
C   = p1 + t*L.*w + (amp.*cos(ph))*u1 + (amp.*sin(ph))*u2;

if strncmpi(p.style,'line',4)
    h = plot3(C(:,1),C(:,2),C(:,3),'-','Color',p.color,'LineWidth',p.linewidth);
    return
end

% --- sweep a circular cross-section along the centerline (parallel transport) ---
T = gradient_(C);
N = zeros(M,3);
n0 = u1 - dot(u1,T(1,:))*T(1,:);
if norm(n0) < 1e-12, n0 = u2 - dot(u2,T(1,:))*T(1,:); end
N(1,:) = n0/norm(n0);
for k = 2:M
    v = cross(T(k-1,:),T(k,:));
    nv = norm(v);
    if nv < 1e-12
        N(k,:) = N(k-1,:);
    else
        v  = v/nv;
        th = atan2(nv, dot(T(k-1,:),T(k,:)));
        N(k,:) = rodrigues_(N(k-1,:), v, th);      % rotation-minimizing frame
    end
    N(k,:) = N(k,:) - dot(N(k,:),T(k,:))*T(k,:);
    nn = norm(N(k,:));
    if nn < 1e-12, N(k,:) = N(k-1,:); else, N(k,:) = N(k,:)/nn; end
end
B = cross(T,N,2);

rw  = p.wireWidth/2;
phi = linspace(0,2*pi,max(4,p.resolution));
X = C(:,1) + rw*(cos(phi).*N(:,1) + sin(phi).*B(:,1));
Y = C(:,2) + rw*(cos(phi).*N(:,2) + sin(phi).*B(:,2));
Z = C(:,3) + rw*(cos(phi).*N(:,3) + sin(phi).*B(:,3));

h = surface(X,Y,Z,'FaceColor',p.color,'EdgeColor','none', ...
            'FaceLighting','gouraud','FaceAlpha',p.facealpha, ...
            'AmbientStrength',.6,'DiffuseStrength',.3,'SpecularStrength',.1);

end  % main function


% ===== local helpers =====
function T = gradient_(C)
% Unit tangents along a 3D curve, by central differences.
T = zeros(size(C));
T(2:end-1,:) = C(3:end,:) - C(1:end-2,:);
T(1,:)       = C(2,:) - C(1,:);
T(end,:)     = C(end,:) - C(end-1,:);
nn = sqrt(sum(T.^2,2)); nn(nn==0) = 1;
T = T ./ nn;
end

function vr = rodrigues_(v, k, th)
% Rotate the row vector v about the unit axis k by th radians.
vr = v*cos(th) + cross(k,v)*sin(th) + k*dot(k,v)*(1-cos(th));
end
