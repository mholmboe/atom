%% hydroxyl_orientation_atom.m
% * Orientation of the structural O-H vectors of a clay/mica layer relative to the
% * layer normal (by default z, i.e. the c* direction of a 00l-oriented cell).
% * For each H the bonded hydroxyl O is found within a short cutoff, and the angle
% * between the O->H vector and the normal is returned.
% *
% * Because an O-H pointing "up" and one pointing "down" are equivalent by symmetry,
% * the angle is folded into [0 90] deg, i.e. theta = acosd(|cos|). theta = 0 means
% * the O-H points along the normal, theta = 90 means it lies in the ab plane. The
% * complement, info.theta_plane = 90 - theta, is the angle to the (001) plane that
% * clay studies usually quote (the OH dipole tilt, rho).
% *
% * In a dioctahedral smectite the O-H points towards the vacant octahedral site and
% * thus lies close to the ab plane (theta near 90, theta_plane small), while in a
% * trioctahedral one it points nearly along the normal (theta near 0).
% *
% * The angles are also reported per hydroxyl-O type in info.by_type, so that for
% * instance Oh (Al-OH-Al) and Ohmg (Al-OH-Mg) are separated.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # theta = hydroxyl_orientation_atom(atom,Box_dim)
% # [theta,angles,info] = hydroxyl_orientation_atom(atom,Box_dim)
% # theta = hydroxyl_orientation_atom(atom,Box_dim,'normal',[0 0 1])
% # theta = hydroxyl_orientation_atom(atom,Box_dim,'oh_cutoff',1.3)
%

function [theta,angles,info] = hydroxyl_orientation_atom(atom,Box_dim,varargin)
%%

% --- options (name/value) ---
p = struct('oh_cutoff',1.3, ...      % max O-H distance in Angstrom
           'h_types',{{}}, ...       % default: any type/element starting with H
           'normal',[0 0 1], ...     % layer normal
           'fold',true);             % fold the angle into [0 90]
for k=1:2:numel(varargin)
    p.(varargin{k}) = varargin{k+1};
end

X=[atom.x]'; Y=[atom.y]'; Z=[atom.z]';
T=[atom.type];
N=numel(atom);
E=T;
if isfield(atom,'element')
    try
        Etmp=[atom.element];
        if numel(Etmp)==N && all(~cellfun(@isempty,Etmp)), E=Etmp; end
    catch
    end
end

lx=Box_dim(1); ly=Box_dim(2); lz=Box_dim(3);
if numel(Box_dim)>=9
    xy=Box_dim(6); xz=Box_dim(8); yz=Box_dim(9);
else
    xy=0; xz=0; yz=0;
end

nrm = p.normal(:)'; nrm = nrm/norm(nrm);

theta=NaN; angles=[];
info=struct('theta_std',NaN,'theta_median',NaN,'theta_plane',NaN, ...
            'theta_signed',NaN,'n_OH',0,'n_H',0,'by_type',struct([]),'note','');

% --- hydrogens and oxygens ---
if ~isempty(p.h_types)
    is_H = ismember(T,p.h_types);
else
    is_H = strcmp(E,'H') | strncmpi(T,'H',1);
end
is_O = strcmp(E,'O') | strncmpi(T,'O',1);
h_idx=find(is_H); o_idx=find(is_O);
info.n_H=numel(h_idx);
if isempty(h_idx) || isempty(o_idx)
    info.note='no hydrogens and/or oxygens found';
    return;
end

Xo=X(o_idx); Yo=Y(o_idx); Zo=Z(o_idx);
ang=[]; signed=[]; otype={};
for ii=1:numel(h_idx)
    h=h_idx(ii);
    d = mic_([Xo-X(h), Yo-Y(h), Zo-Z(h)], lx,ly,lz,xy,xz,yz);   % H->O vectors
    r = sqrt(sum(d.^2,2));
    [rmin,kk] = min(r);
    if rmin >= p.oh_cutoff, continue; end                        % unbound H, e.g. water
    vOH = -d(kk,:);                                              % O->H
    nv  = norm(vOH);
    if nv==0, continue; end
    c = dot(vOH/nv, nrm);
    signed(end+1) = acosd(max(-1,min(1,c)));                     %#ok<AGROW> 0..180
    if p.fold
        ang(end+1) = acosd(min(1,abs(c)));                       %#ok<AGROW> 0..90
    else
        ang(end+1) = signed(end);                                %#ok<AGROW>
    end
    otype{end+1} = T{o_idx(kk)};                                 %#ok<AGROW>
end

if isempty(ang)
    info.note='no O-H pairs found within oh_cutoff';
    return;
end

angles = ang(:);
theta  = mean(angles);
info.n_OH        = numel(angles);
info.theta_std   = std(angles,1);
info.theta_median= median(angles);
info.theta_plane = 90 - theta;          % angle to the (001) plane
info.theta_signed= mean(signed);

% --- breakdown per hydroxyl-O type ---
ut = unique(otype);
bt = struct('type',{},'theta',{},'theta_std',{},'theta_plane',{},'n',{});
for kk=1:numel(ut)
    sel = strcmp(otype, ut{kk});
    bt(kk).type        = ut{kk};
    bt(kk).theta       = mean(angles(sel));
    bt(kk).theta_std   = std(angles(sel),1);
    bt(kk).theta_plane = 90 - mean(angles(sel));
    bt(kk).n           = sum(sel);
end
info.by_type = bt;

end  % main function


% ===== local helpers =====
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
