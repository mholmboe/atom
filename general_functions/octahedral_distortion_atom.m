%% octahedral_distortion_atom.m
% * Distortion parameters of the octahedrally coordinated cations of a clay/mica
% * layer, i.e. the Alo-Op/Oh and Mgo-Omg/Ohmg coordination of a smectite.
% *
% * Per octahedron the following standard measures are computed:
% *   MO      mean M-O bond length (Angstrom)
% *   Delta   Baur's bond length distortion index, mean(|d_i - <d>|)/<d>, 0 = regular
% *   lambda  Robinson's quadratic elongation, mean((d_i/d0)^2), 1 = regular, where d0
% *           is the centre-to-vertex distance of a regular octahedron of equal volume
% *   sigma2  Robinson's bond angle variance, sum((theta_i-90)^2)/11 over the 12 cis
% *           O-M-O angles (deg^2), 0 = regular
% *   psi     octahedral flattening angle, cos(psi) = t_oct/(2*<M-O>), ideal 54.74 deg,
% *           with t_oct the separation of the two O triangles along the layer normal
% *
% * Robinson, Gibbs & Ribbe, Science 172 (1971) 567; Baur, Acta Cryst B30 (1974) 1195.
% *
% * Results are returned for all octahedra together and, in info.by_cation, split per
% * cation type (Al, Mgo, Feo, ...). info.by_ligand gives the mean M-O distance per
% * coordinating oxygen type, so Alo-Op and Alo-Oh are separated.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # Delta = octahedral_distortion_atom(atom,Box_dim)
% # [Delta,per_oct,info] = octahedral_distortion_atom(atom,Box_dim)
% # Delta = octahedral_distortion_atom(atom,Box_dim,'oct_cutoff',2.5)
% # Delta = octahedral_distortion_atom(atom,Box_dim,'oct_types',{'Al','Mgo'})
%

function [Delta,per_oct,info] = octahedral_distortion_atom(atom,Box_dim,varargin)
%%

% --- options (name/value) ---
p = struct('oct_cutoff',2.5, ...     % max M-O distance in Angstrom
           'oct_types',{{}}, ...     % default: the oct_elements prefixes below
           'oct_elements',{{'Al','Mg','Fe','Li','Ti','Mn','Ni','Co','Cr','Zn'}}, ...
           'tet_types',{{'Si','Sit','Alt','Tit','Fet','Fee3'}}, ...
           'normal',[0 0 1], ...     % layer normal, for psi
           'min_coord',5, 'max_coord',7);
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
nrm = p.normal(:); nrm = nrm/norm(nrm);

Delta=NaN; per_oct=[];
info=struct('MO',NaN,'MO_std',NaN,'Delta',NaN,'Delta_std',NaN, ...
            'lambda',NaN,'lambda_std',NaN,'sigma2',NaN,'sigma2_std',NaN, ...
            'psi',NaN,'psi_std',NaN,'n_oct',0,'n_oct_skipped',0, ...
            'by_cation',struct([]),'by_ligand',struct([]),'note','');

% --- pick the octahedral cations ---
is_tet = ismember(T,p.tet_types);
if ~isempty(p.oct_types)
    is_oct = ismember(T,p.oct_types);
else
    is_oct = false(1,N);
    for kk=1:numel(p.oct_elements)
        e=p.oct_elements{kk};
        is_oct = is_oct | strncmpi(T,e,numel(e)) | strncmpi(E,e,numel(e));
    end
end
is_oct = is_oct & ~is_tet;
is_O   = strcmp(E,'O') | strncmpi(T,'O',1);
m_idx=find(is_oct); o_idx=find(is_O);
if isempty(m_idx) || isempty(o_idx)
    info.note='no octahedral cations and/or oxygens found';
    return;
end

Xo=X(o_idx); Yo=Y(o_idx); Zo=Z(o_idx);
MO=[]; DD=[]; LA=[]; SG=[]; PS=[]; ctype={}; nskip=0;
lig_type={}; lig_d=[]; lig_cat={};
for ii=1:numel(m_idx)
    M=m_idx(ii);
    d = mic_([Xo-X(M), Yo-Y(M), Zo-Z(M)], lx,ly,lz,xy,xz,yz);
    r = sqrt(sum(d.^2,2));
    sel = find(r < p.oct_cutoff);
    nc  = numel(sel);
    if nc < p.min_coord || nc > p.max_coord, nskip=nskip+1; continue; end

    v  = d(sel,:);            % M->O vectors
    dd = r(sel);              % M-O distances
    dm = mean(dd);
    if dm<=0, nskip=nskip+1; continue; end

    MO(end+1)    = dm;                                   %#ok<AGROW>
    DD(end+1)    = mean(abs(dd-dm))/dm;                  %#ok<AGROW> Baur
    ctype{end+1} = T{M};                                 %#ok<AGROW>

    % --- ligand bookkeeping, so Alo-Op and Alo-Oh can be separated ---
    for kk=1:nc
        lig_type{end+1} = T{o_idx(sel(kk))};             %#ok<AGROW>
        lig_d(end+1)    = dd(kk);                        %#ok<AGROW>
        lig_cat{end+1}  = T{M};                          %#ok<AGROW>
    end

    % --- Robinson quadratic elongation, only for a proper 6-coordination ---
    if nc==6
        Vol = octvol_(v);
        if ~isnan(Vol) && Vol>0
            % A regular octahedron with centre-to-vertex distance a has V = 4/3*a^3,
            % so the equal-volume reference bond length is d0 = (3V/4)^(1/3).
            d0 = (3*Vol/4)^(1/3);
            LA(end+1) = mean((dd/d0).^2);                %#ok<AGROW>
        else
            LA(end+1) = NaN;                             %#ok<AGROW>
        end
        % --- Robinson bond angle variance over the 12 cis angles ---
        na=[];
        for a1=1:nc-1
            for a2=a1+1:nc
                na(end+1) = acosd(max(-1,min(1, ...
                    dot(v(a1,:),v(a2,:))/(dd(a1)*dd(a2))))); %#ok<AGROW>
            end
        end
        na = sort(na);                                   % 15 angles: 12 cis + 3 trans
        cis = na(1:12);
        SG(end+1) = sum((cis-90).^2)/11;                 %#ok<AGROW>
    else
        LA(end+1) = NaN;  SG(end+1) = NaN;               %#ok<AGROW>
    end

    % --- flattening angle psi ---
    proj = sort(v*nrm);
    k    = floor(nc/2);
    t_oct= mean(proj(end-k+1:end)) - mean(proj(1:k));
    PS(end+1) = acosd(min(1,max(0, t_oct/(2*dm))));      %#ok<AGROW>
end

info.n_oct_skipped = nskip;
if isempty(MO)
    info.note='no octahedra with an acceptable coordination number';
    return;
end

Delta = mean(DD);
per_oct = [MO(:), DD(:), LA(:), SG(:), PS(:)];     % one row per octahedron
info.n_oct      = numel(MO);
info.MO         = mean(MO);        info.MO_std     = std(MO,1);
info.Delta      = mean(DD);        info.Delta_std  = std(DD,1);
info.lambda     = mean(LA,'omitnan');  info.lambda_std = std(LA(~isnan(LA)),1);
info.sigma2     = mean(SG,'omitnan');  info.sigma2_std = std(SG(~isnan(SG)),1);
info.psi        = mean(PS);        info.psi_std    = std(PS,1);

% --- per cation type ---
uc = unique(ctype);
bc = struct('type',{},'MO',{},'Delta',{},'lambda',{},'sigma2',{},'psi',{},'n',{});
for kk=1:numel(uc)
    sel = strcmp(ctype, uc{kk});
    bc(kk).type   = uc{kk};
    bc(kk).MO     = mean(MO(sel));
    bc(kk).Delta  = mean(DD(sel));
    bc(kk).lambda = mean(LA(sel),'omitnan');
    bc(kk).sigma2 = mean(SG(sel),'omitnan');
    bc(kk).psi    = mean(PS(sel));
    bc(kk).n      = sum(sel);
end
info.by_cation = bc;

% --- per cation-ligand pair, e.g. Al-Op, Al-Oh, Mgo-Omg, Mgo-Ohmg ---
pair = strcat(lig_cat, '-', lig_type);
up   = unique(pair);
bl   = struct('pair',{},'MO',{},'MO_std',{},'n',{});
for kk=1:numel(up)
    sel = strcmp(pair, up{kk});
    bl(kk).pair   = up{kk};
    bl(kk).MO     = mean(lig_d(sel));
    bl(kk).MO_std = std(lig_d(sel),1);
    bl(kk).n      = sum(sel);
end
info.by_ligand = bl;

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

function V = octvol_(v)
% Volume of the polyhedron spanned by the 6 M->O vectors, via its convex hull.
V=NaN;
try
    [~,V] = convhull(v(:,1),v(:,2),v(:,3));
catch
end
end
