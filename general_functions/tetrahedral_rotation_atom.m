%% tetrahedral_rotation_atom.m
% * Tetrahedral rotation angle alpha (deg) — the ditrigonal distortion of a
% * phyllosilicate/clay/mica tetrahedral sheet (Bailey 1984, Rev. Mineral. 19;
% * Radoslovich 1961). alpha = 0 is the ideal hexagonal siloxane ring; the
% * geometric maximum is 30 deg; smectites/micas fall in ~1-15 deg
% * (dioctahedral > trioctahedral). MATLAB/Octave port of
% * atomipy.distortion.tetrahedral_rotation.
% *
% * Method (per-tetrahedron, ring-perception-free). For each tetrahedral cation
% * (Si/Sit/Alt) the basal (bridging) and apical oxygens are found by atom-type
% * (Ob/Op) with a connectivity fallback (an O bonded to 2 tetrahedral cations is
% * basal, to 1 is apical). The local sheet normal is the metal->apical direction
% * (so the result is orientation/tilt independent, and the basal-triplet plane
% * normal is checked against it). alpha is then the in-plane angle between the
% * metal->basalO bond and the metal->neighbour-metal line — this equals the
% * canonical (1/2)*mean|120-phi| ring definition.
%
%% Version
% 1.00
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # alpha = tetrahedral_rotation_atom(atom,Box_dim)
% # [alpha,angles,info] = tetrahedral_rotation_atom(atom,Box_dim,'bond_cutoff',1.9)
% # alpha = tetrahedral_rotation_atom(atom,Box_dim,'tet_types',{'Sit','Alt'})

function [alpha,angles,info] = tetrahedral_rotation_atom(atom,Box_dim,varargin)
%%

% --- options (name/value) ---
p = struct('bond_cutoff',1.9, ...
           'tet_types',{{'Si','Sit','Alt','Tit','Fee3'}}, ...
           'tet_elements',{{'Si'}}, ...
           'basal_types',{{}}, 'apical_types',{{}}, ...
           'align_tol',0.5, 'tet_tet_cutoff',3.6);
for k=1:2:numel(varargin)
    p.(varargin{k}) = varargin{k+1};
end

X=[atom.x]'; Y=[atom.y]'; Z=[atom.z]';
T=[atom.type];
N=numel(atom);
E=T;                          % default: use the atom type as an element proxy
if isfield(atom,'element')    % .pdb import provides elements; .gro import usually does not
    try
        Etmp=[atom.element];
        if numel(Etmp)==N && all(~cellfun(@isempty,Etmp)), E=Etmp; end
    catch
    end
end

lx=Box_dim(1); ly=Box_dim(2); lz=Box_dim(3);
if numel(Box_dim)>=9
    % Box_dim=[lx ly lz 0 0 xy 0 xz yz]
    xy=Box_dim(6); xz=Box_dim(8); yz=Box_dim(9);
else
    xy=0; xz=0; yz=0;
end

% --- select tetrahedral cations and oxygens ---
is_tet = ismember(T,p.tet_types) | ismember(E,p.tet_elements);
is_O   = strcmp(E,'O') | strncmp(T,'O',1);
tet_idx = find(is_tet);
o_idx   = find(is_O);

alpha=NaN; angles=[]; info=struct('n_tet',numel(tet_idx),'n_tet_used',0, ...
    'n_bonds',0,'mean_triplet_alignment',NaN,'alpha_std',NaN,'alpha_median',NaN, ...
    'apical_tilt',NaN,'apical_tilt_std',NaN,'tau',NaN,'tau_std',NaN, ...
    'dz_corrugation',NaN,'dz_corrugation_std',NaN,'n_sheets',0, ...
    'psi',NaN,'psi_std',NaN,'n_oct',0,'note','');
if isempty(tet_idx) || isempty(o_idx)
    info.note='no tetrahedral cations and/or oxygens found';
    return;
end

% --- bond graph (minimum image): o_tet{o}=[tet vx vy vz]; tet_o{t}=[o vx vy vz] ---
o_tet = cell(N,1);
tet_o = cell(N,1);
Xo=X(o_idx); Yo=Y(o_idx); Zo=Z(o_idx);
for ii=1:numel(tet_idx)
    t=tet_idx(ii);
    d = mic_([Xo-X(t), Yo-Y(t), Zo-Z(t)], lx,ly,lz,xy,xz,yz);
    r = sqrt(sum(d.^2,2));
    near = find(r < p.bond_cutoff);
    for kk=1:numel(near)
        o = o_idx(near(kk));
        o_tet{o}(end+1,:) = [t, d(near(kk),:)];
        tet_o{t}(end+1,:) = [o, d(near(kk),:)];
    end
end

% --- per-tetrahedron rotation ---
aligns=[]; n_used=0; global_normal=[0 0 1];
apical_tilts=[]; tau_angles=[]; basal_z=zeros(N,1); basal_seen=false(N,1);
for ii=1:numel(tet_idx)
    M=tet_idx(ii);
    nb=tet_o{M};
    if isempty(nb), continue; end
    isb=false(size(nb,1),1); isa=false(size(nb,1),1);
    for r=1:size(nb,1)
        rl = role_(char(atom(nb(r,1)).type), size(o_tet{nb(r,1)},1), p);
        isb(r)=strcmp(rl,'basal'); isa(r)=strcmp(rl,'apical');
    end
    basal=nb(isb,:); apic=nb(isa,:);
    if size(basal,1)<3 || isempty(apic), continue; end
    % local sheet normal = metal->apical
    apex=mean(apic(:,2:4),1); na=norm(apex);
    if na>0, nM=apex/na; else nM=global_normal; end
    % quality: basal-triplet plane normal vs metal->apical
    b3=basal(1:3,2:4);
    tri=cross(b3(2,:)-b3(1,:), b3(3,:)-b3(1,:)); ntri=norm(tri);
    if ntri>0, align=abs(dot(tri/ntri,nM)); else align=0; end
    aligns(end+1)=align;
    if align < p.align_tol, continue; end
    n_used=n_used+1;
    % --- companion parameters (same tetrahedron) ---
    apical_tilts(end+1) = acosd(min(1, abs(dot(nM, global_normal))));  % apical bond vs normal
    ov = [basal(:,2:4); apic(:,2:4)];                                  % M->O vectors
    for ia=1:size(ov,1)
        for ib=ia+1:size(ov,1)
            nn=norm(ov(ia,:))*norm(ov(ib,:));
            if nn>0, tau_angles(end+1)=acosd(max(-1,min(1,dot(ov(ia,:),ov(ib,:))/nn))); end
        end
    end
    for r=1:size(basal,1)
        oi=basal(r,1);
        basal_z(oi)=dot([X(oi) Y(oi) Z(oi)], global_normal); basal_seen(oi)=true;
    end
    for r=1:size(basal,1)
        o=basal(r,1); vMO=basal(r,2:4);
        links=o_tet{o};
        oth=links(links(:,1)~=M,:);
        if isempty(oth), continue; end
        vM2O=oth(1,2:4);
        vMM=vMO-vM2O;            % M->M' anchored on the shared oxygen
        if norm(vMM) > p.tet_tet_cutoff, continue; end
        ang=angle_in_plane_(vMO,vMM,nM);
        if ~isnan(ang), angles(end+1)=ang; end
    end
end

if isempty(angles)
    info.note='no valid tetrahedra (need 3 basal + 1 apical O per cation)';
    return;
end
angles=angles(:);
alpha=mean(angles);
info.n_tet_used=n_used;
info.n_bonds=numel(angles);
info.alpha_std=std(angles,1);
info.alpha_median=median(angles);
if ~isempty(aligns), info.mean_triplet_alignment=mean(aligns); end

% --- companion distortion parameters ---
if ~isempty(apical_tilts)
    info.apical_tilt=mean(apical_tilts); info.apical_tilt_std=std(apical_tilts,1);
end
if ~isempty(tau_angles)
    info.tau=mean(tau_angles); info.tau_std=std(tau_angles,1);          % ideal 109.47 deg
end
[info.dz_corrugation, info.dz_corrugation_std, info.n_sheets] = corrugation_(basal_z(basal_seen));
[info.psi, info.psi_std, info.n_oct] = octahedral_flattening_(X,Y,Z,T,E, ...
    tet_idx, global_normal, lx,ly,lz,xy,xz,yz);   % ideal 54.74 deg

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

function rl = role_(t, ntet, p)
% basal / apical from atom-type prefix (Ob*/Op*/Omg) with a bridging-count fallback.
if ~isempty(p.apical_types) && any(strcmp(t,p.apical_types)), rl='apical'; return; end
if ~isempty(p.basal_types)  && any(strcmp(t,p.basal_types)),  rl='basal';  return; end
if numel(t)>=2 && strcmpi(t(1:2),'Ob'), rl='basal';  return; end
if (numel(t)>=2 && strcmpi(t(1:2),'Op')) || strcmpi(t,'Omg'), rl='apical'; return; end
if ntet>=2, rl='basal'; elseif ntet==1, rl='apical'; else rl='other'; end
end

function a = angle_in_plane_(v1, v2, n)
% Angle (deg) between v1 and v2 after projecting both onto the plane normal n.
p1=v1-dot(v1,n)*n; p2=v2-dot(v2,n)*n;
n1=norm(p1); n2=norm(p2);
if n1==0 || n2==0, a=NaN; return; end
c=dot(p1,p2)/(n1*n2);
a=acosd(max(-1,min(1,c)));
end

function [dz, dz_std, nsheets] = corrugation_(coords, gap)
% Basal-oxygen corrugation Delta_z (Angstrom): cluster the sheet-normal coordinates into
% sheets (split where the sorted gap exceeds 'gap'), take max-min within each sheet, and
% return (mean, std across sheets, n_sheets).
if nargin<2, gap=2.0; end
c=sort(coords(:));
if numel(c)<3, dz=NaN; dz_std=NaN; nsheets=0; return; end
splits=find(diff(c)>gap);
starts=[1; splits+1]; ends=[splits; numel(c)];
dzs=[];
for g=1:numel(starts)
    seg=c(starts(g):ends(g));
    if numel(seg)>=3, dzs(end+1)=max(seg)-min(seg); end
end
if isempty(dzs), dz=NaN; dz_std=NaN; nsheets=0;
else dz=mean(dzs); dz_std=std(dzs,1); nsheets=numel(dzs); end
end

function [psi, psi_std, noct] = octahedral_flattening_(X,Y,Z,T,E, tet_idx, normal, lx,ly,lz,xy,xz,yz)
% Mean octahedral flattening angle psi (deg): cos(psi)=t_oct/(2*<M-O>), t_oct = octahedral
% sheet thickness (separation of the two O triangles along the normal). Ideal psi=54.74.
OCT={'Al','Mg','Fe','Li','Ti','Mn','Ni','Co','Cr','Zn','Ca'};
N=numel(X);
isO = strcmp(E,'O') | strncmp(T,'O',1);
o_idx=find(isO); Xo=X(o_idx); Yo=Y(o_idx); Zo=Z(o_idx);
tetset=false(N,1); tetset(tet_idx)=true;
normal=normal(:);
psis=[];
for M=1:N
    if tetset(M), continue; end
    if ~any(strncmpi(E{M}, OCT, 2)), continue; end   % octahedral element (Al/Mg/Fe/... prefix)
    d=mic_([Xo-X(M), Yo-Y(M), Zo-Z(M)], lx,ly,lz,xy,xz,yz);
    r=sqrt(sum(d.^2,2));
    sel=find(r<2.5);
    if numel(sel)<5 || numel(sel)>7, continue; end   % ~octahedral (6, allow 5-7)
    v=d(sel,:);
    mo=mean(sqrt(sum(v.^2,2)));
    proj=sort(v*normal);
    k=floor(numel(sel)/2);
    t_oct=mean(proj(end-k+1:end))-mean(proj(1:k));
    if mo>0, psis(end+1)=acosd(min(1,max(0, t_oct/(2*mo)))); end
end
if isempty(psis), psi=NaN; psi_std=NaN; noct=0;
else psi=mean(psis); psi_std=std(psis,1); noct=numel(psis); end
end
