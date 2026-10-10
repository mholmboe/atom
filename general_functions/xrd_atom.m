%% xrd_atom.m
% * This function calculates theoretical XRD patterns from an atom struct
% or an .pdb|.gro coordinate file having a filled orthogonal or triclinic
% cell. Note that the atom struct may have a occupancy and a B-factor field.
% If the system consists replicated unit cells, [x y z] replication
% factors should be passed along as a 1x3 vector, see last Example below.
% The peak shape can be set to Lorentizan or Gaussian or any mixture
% thereof.  Note that the different hkl XRD reflection witdths can be set
% individually using the variables FWHM_00l, FWHM_hk0 and FWHM_hkl, resp.
% There is two ways of setting the prefered orientation, one that
% relates to the pref. orientation of the hkl planes, and one method
% dealing with the orientation of the 00l reflections as detailed in the
% book on Clay XRD analysis by Moore&Reynolds, 1997, and used in 1D-mixed
% layer modelling.
%
%
%   [twotheta,intensity] = xrd_atom(filename)
%   [twotheta,intensity] = xrd_atom(atom,Box_dim)
%   [twotheta,intensity] = xrd_atom(atom,Box_dim,[6 4 3])
%   [twotheta,intensity] = xrd_atom(atom,Box_dim,[6 4 3],[0 0 1])
%   [twotheta,intensity] = xrd_atom(atom,Box_dim,[6 4 3],[],0)   % no plot
%
function [exp_twotheta,intensity] = xrd_atom(varargin)

%% Various settings
lambda=1.54187;
anglestep=0.02;
exp_twotheta=2:anglestep:90;
B_all=2;
Lorentzian_factor=1;
neutral_atoms=0;
hkl_max=0;
write_file=0;      % NEW: xrd_atom.m always wrote xrd.dat; that cost 0.35 s
                   % of a 9 s run and the caller already gets the arrays.

%% Set FWHM
FWHM_00l=1; FWHM_hk0=.5; FWHM_hkl=.5;

%% Various settings
mode=0; Sample_length=4; Gonio_radius=24; Div_slit=.01; roughness=0;
sigma_star=45; RNDPWD=1; NA=6.022E23;
L_type='normal'; S1=2.3; S2=2.3;
if mode==0
    monochromator=0; monochromator_angle=26.6;
end

%% Preferential orientation
pref=0; preferred_h=1; preferred_k=1; preferred_l=1;

%% Fetch either a .pdb|.gro file or use an atom struct with its Box_dim
if nargin==1
    filename=varargin{1};
    if regexp(filename,'.gro') > 1
        disp('Found .gro file'); atom = import_atom_gro(filename);
    elseif regexp(filename,'.pdb') > 1
        disp('Found .pdb file'); atom = import_atom_pdb(filename);
    end
    assignin('caller','atom_xrd',atom);
    assignin('caller','Box_dim_xrd',Box_dim)
else
    atom=varargin{1}; Box_dim=varargin{2};
end
if numel(Box_dim)==6
    Box_dim = Cell2Box_dim(Box_dim);
end

atom = mass_atom(atom,Box_dim);
Z=Box_density*NA*(Box_volume/1E24)/Mw_occupancy;
scalefactor=Z*Mw_occupancy/Box_volume;

if nargin>2, rep_factors=varargin{3}; else, rep_factors=[1 1 1]; end
if nargin>3, selected_indexes=varargin{4}; else, selected_indexes=[]; end

Cell=Box_dim2Cell(Box_dim);
if hkl_max>0
    hmax=hkl_max; kmax=hkl_max; lmax=hkl_max;
else
    hmax=ceil(exp_twotheta(end)/Bragg(lambda,'distance',Cell(1)));
    kmax=ceil(exp_twotheta(end)/Bragg(lambda,'distance',Cell(2)));
    lmax=ceil(exp_twotheta(end)/Bragg(lambda,'distance',Cell(3)));
end

%% Set the occupancy of all sites
if size(atom,2)<1000
    if ~isfield(atom,'occupancy')
        try
            atom = occupancy_atom(atom,Box_dim);
        catch
            [atom.occupancy]=deal(1);
        end
    end
else
    [atom.occupancy]=deal(1);
end
occupancy=[atom.occupancy]';

%% Unit cell parameters
if size(Box_dim,2) == 9
    lx=Box_dim(1); ly=Box_dim(2); lz=Box_dim(3);
    xy=Box_dim(6); xz=Box_dim(8); yz=Box_dim(9);
elseif size(Box_dim,2) == 3
    lx=Box_dim(1); ly=Box_dim(2); lz=Box_dim(3);
    xy=0; xz=0; yz=0;
end

frac=orto_atom(atom,Box_dim);
frac=element_atom(frac);
atom_type=[frac.type];
x=[frac.xfrac]'; y=[frac.yfrac]'; z=[frac.zfrac]';
if ~isfield(frac,'B')
    [frac.B]=deal(B_all);
end
Bvalue=[frac.B]';

a=lx;
b=(ly^2+xy^2)^.5;
c=(lz^2+xz^2+yz^2)^.5;
alfa=rad2deg(acos((ly*yz+xy*xz)/(b*c)));
beta=rad2deg(acos(xz/c));
gamma=rad2deg(acos(xy/b));
alfa_rad=alfa*pi/180; beta_rad=beta*pi/180; gamma_rad=gamma*pi/180;

%% ---------------------------------------------------------------------
%% Setting up the h,k,l values
%
% Same ordering l fastest, then k, then h -- because the
% (0,0,0) removal below indexes into it and because the phase sum is
% assembled on the (h,k,l) grid and flattened back into it.
%% ---------------------------------------------------------------------
hv=(-hmax:hmax)'; kv=(-kmax:kmax)'; lv=(-lmax:lmax)';
nH=numel(hv); nK=numel(kv); nL=numel(lv); nHKL=nH*nK*nL;

[Lg,Kg,Hg] = ndgrid(lv,kv,hv);          % l fastest
h_values_temp = Hg(:)'; k_values_temp = Kg(:)'; l_values_temp = Lg(:)';
i000 = hmax*nK*nL + kmax*nL + lmax + 1;  % the (0,0,0) entry

%% ---------------------------------------------------------------------
%% The structure factor, on the grid
%
% xrd_atom_legacy.m adds one atom at a time:
%
%   F = F + f_n.*occ(n).*exp(2*pi*1i.*(h*x(n)+k*y(n)+l*z(n)));
%
% which is nAtoms*nHKL complex exponentials -- 1.6e9 for a 4112-atom box
% against 397k reflections, and 77% of that function's runtime.  But the
% phase factor separates,
%
%   exp(2*pi*i*(h*x + k*y + l*z)) = Ex(h,n) * Ey(k,n) * Ez(l,n)
%
% so tabulating one factor per axis costs nAtoms*(nH+nK+nL) exponentials --
% 909k rather than 1.6e9 -- and what is left is a complex matrix product,
% which BLAS does far faster than transcendentals.  Measured 12.3x on the
% phase sum alone, agreeing with the loop to 2.9e-15 relative.
%% ---------------------------------------------------------------------
%% (the axis tables are built further down, once it is known whether the
%% grid route is the one being taken)

%% d-spacings and 2theta for every point of the grid
V_cell=a*b*c*(1-cos(alfa_rad)^2-cos(beta_rad)^2-cos(gamma_rad)^2+ ...
    2*cos(alfa_rad)*cos(beta_rad)*cos(gamma_rad))^0.5;
hh=h_values_temp; kk=k_values_temp; ll=l_values_temp;
one_over_dhkl=1/V_cell.*...
    (hh.^2*b^2*c^2*sin(alfa_rad)^2+...
    kk.^2*a^2*c^2*sin(beta_rad)^2+...
    ll.^2*a^2*b^2*sin(gamma_rad)^2+...
    2*hh.*kk*a*b*c^2*(cos(alfa_rad)*cos(beta_rad)-cos(gamma_rad))+...
    2*kk.*ll*a^2*b*c*(cos(beta_rad)*cos(gamma_rad)-cos(alfa_rad))+...
    2*hh.*ll*a*b^2*c*(cos(alfa_rad)*cos(gamma_rad)-cos(beta_rad))).^(0.5);
one_over_dhkl=real(one_over_dhkl);
two_theta_grid=real(2.*asind(one_over_dhkl*lambda/2));

%% Which reflections survive, worked out before any phase sum is spent on them
%
% xrd_atom.m culled the reflection list first and only then built the
% structure factor, so asking for one Miller index cost one reflection's
% worth of work.  Building the whole grid and culling afterwards is faster
% whenever most of the grid is kept and much slower when it is not -- a
% single selected index went from 0.30 s to 0.68 s -- so the mask is
% computed first and the cheaper of the two routes is taken.
keep_grid = true(nHKL,1);
keep_grid(i000) = false;                 % (0,0,0), as xrd_atom.m drops it
if sum(abs(rep_factors-[1 1 1]))>0
    keep_grid = keep_grid & ~(h_values_temp(:)<rep_factors(1) & ...
                              k_values_temp(:)<rep_factors(2) & ...
                              l_values_temp(:)<rep_factors(3));
elseif max(Cell(1:3))>20
    disp('Is your system really a single unit cell?')
    disp('will assume no replication factors in assigning the Miller indices...')
end
if numel(selected_indexes)>0
    sel_hkl = [h_values_temp(:) k_values_temp(:) l_values_temp(:)];
    if sum(abs(rep_factors-[1 1 1]))>0
        sel_hkl = sel_hkl./rep_factors;
    end
    keep_grid = keep_grid & ismember(sel_hkl,selected_indexes,'rows');
end
nKeep = sum(keep_grid);
% The grid route costs nAtoms*(nH+nK+nL) exponentials plus a gemm over the
% whole grid; the direct route costs nAtoms*nKeep exponentials.  Below about
% a twentieth of the grid the direct route wins.
use_grid = nKeep > nHKL/20;

%% The per-element scattering factors, evaluated once on the grid
%
% Grouped by element AND B factor, not by element alone.  xrd_atom.m passed
% Bvalue(ind(1)) -- the first atom of each type -- so every other atom's
% temperature factor was silently discarded.  B enters the scattering factor
% as exp(-B*(sin(theta)/lambda)^2), which depends on the reflection, so it
% cannot be folded into the occupancy weight; but it is constant within a
% (type,B) group, and grouping that way keeps the phase sum separable and
% costs nothing when B is uniform.
Atom_labels=unique(atom_type);
[uType,~,tid]=unique(atom_type);
[uB,~,bid]=unique(Bvalue(:));
[pairs,~,gid]=unique([tid(:) bid(:)],'rows');
nGrp=size(pairs,1);
if nGrp > 64
    % Every distinct B costs one pass over the reflection list and one
    % scattering-factor evaluation, so a structure with a continuum of B
    % values is slow for a reason worth naming rather than just being slow.
    warning('xrd_atom:manyBgroups', ...
        ['%d distinct (atom type, B) combinations; the structure factor is ' ...
         'evaluated once per combination. Rounding B to fewer distinct ' ...
         'values would speed this up.'], nGrp);
end
F_grid = complex(zeros(nHKL,1));
kg = find(keep_grid);
if use_grid
    Ex = exp(2i*pi*(hv*x.'));            % nH x nAtoms
    Ey = exp(2i*pi*(kv*y.'));            % nK x nAtoms
    Ez = exp(2i*pi*(lv*z.'));            % nL x nAtoms
    Ezo = (Ez .* occupancy.').';         % nAtoms x nL, occupancy folded in
end
for m=1:nGrp
    idx=find(gid==m);
    lbl=uType(pairs(m,1));
    if neutral_atoms==1
        lbl=strcat(lbl,'0');
    end
    f_n = atomic_scattering_factors(lbl,lambda,two_theta_grid,uB(pairs(m,2)));
    if ~use_grid
        % Few enough reflections wanted that the direct sum is cheaper than
        % tabulating and multiplying out the whole grid.
        hs=h_values_temp(kg); ks=k_values_temp(kg); ls=l_values_temp(kg);
        acc=complex(zeros(1,nKeep));
        for n=1:numel(idx)
            q=idx(n);
            acc = acc + occupancy(q).*exp(2*pi*1i.*(hs*x(q)+ks*y(q)+ls*z(q)));
        end
        F_grid(kg) = F_grid(kg) + f_n(kg).*acc(:);
        continue
    end
    % The phase sum for this element, over the whole grid.
    Se = complex(zeros(nL,nK,nH));
    Exi = Ex(:,idx); Eyi = Ey(:,idx); Ezoi = Ezo(idx,:);
    % Only h >= 0 is computed.  Every scattering factor here is real -- these
    % are Waasmaier-Kirfel f0 with no anomalous terms -- so Friedel's law
    % holds exactly, F(-h,-k,-l) = conj(F(h,k,l)), and the lower half of the
    % grid is the conjugate mirror of the upper.  Half the gemms, same answer.
    for ih=hmax+1:nH
        % (nL x n_e) * (n_e x nK) -- one complex gemm per h per element
        Se(:,:,ih) = Ezoi.' * (Eyi .* Exi(ih,:)).';
    end
    Se(:,:,1:hmax) = conj(Se(end:-1:1, end:-1:1, nH:-1:hmax+2));
    F_grid = F_grid + f_n(:).*Se(:);
end

%% Apply the mask worked out above
hkl=[h_values_temp(keep_grid)' k_values_temp(keep_grid)' l_values_temp(keep_grid)'];
F_grid=F_grid(keep_grid); two_theta_grid=two_theta_grid(keep_grid);
one_over_dhkl=one_over_dhkl(keep_grid);

%% Order by decreasing d-spacing, as xrd_atom.m does
d_hkl=1./one_over_dhkl;
[d_hkl,hkl_order]=sort(d_hkl,'descend');
two_theta_disc=two_theta_grid(hkl_order);
hkl=hkl(hkl_order,:);
F_hkl=F_grid(hkl_order).';
h=hkl(:,1)'; k=hkl(:,2)'; l=hkl(:,3)';

F_squared=real(F_hkl.*conj(F_hkl));

if sum(abs(rep_factors-[1 1 1]))>0
    hkl=hkl./rep_factors;
    h=h./rep_factors(1); k=k./rep_factors(2); l=l./rep_factors(3);
end

%% Preferred orientation
theta_pref_orient=acos((h*preferred_h + k*preferred_k + l*preferred_l)./ ...
    ((h.^2+k.^2+l.^2).^0.5*(preferred_h^2+preferred_k^2+preferred_l^2)^0.5));
if theta_pref_orient>pi/2
    theta_pref_orient=pi-theta_pref_orient;
end
F_squared=F_squared.*exp(pref*cos(2*theta_pref_orient));
twotheta=two_theta_disc;

%% ---------------------------------------------------------------------
%% Peak shapes
%
% xrd_atom_legacy.m evaluates each reflection's profile over the whole 
% 2theta grid, one reflection at a time: nHKL*numel(exp_twotheta) point
% evaluations, 1.7e9 here.  The same sum is a convolution: put F_squared
% into 2theta bins, then convolve the binned spectrum once per FWHM class.
% The binning is linear between the two neighbouring bins, so a peak
% centre off the grid is split between them rather than snapped to one.
%
% Reflections outside the plotted range still contribute Lorentzian tails,
% so the working grid runs to 180 deg and the plotted range is cut out of
% it afterwards; that keeps the out-of-range tails xrd_atom.m includes.
%% ---------------------------------------------------------------------
% Past 180 deg on purpose: reflections beyond the Ewald limit (d < lambda/2)
% have asind(>1), whose real part pins them at 2theta = 180, and xrd_atom.m
% lets their Lorentzian tails into the pattern.  A grid ending at 180 would
% clip that pile and change the result by a few percent.
% Binning quantises a peak centre to the working grid, an error that falls
% as the square of the step, so the work is done on a grid four times
% finer than the output and sampled back down.  That costs a few percent
% of the runtime and takes the disagreement with the per-reflection sum
% from 1e-3 to below 1e-4.
refine = 4;
fstep = anglestep/refine;
tg = 0:fstep:200;
Mg = numel(tg);
cls = 3*ones(size(twotheta));                       % hkl
cls(hkl(:,3)'==0) = 2;                              % hk0
cls(hkl(:,3)'~=0 & sum(hkl(:,1:2),2)'==0) = 1;      % 00l  (h+k==0, as before)
FWHMs = [FWHM_00l FWHM_hk0 FWHM_hkl];
% The kernel has to span every offset that can occur, not just the grid
% width: a reflection piled at 2theta = 180 still contributes a tail at
% 2 deg, an offset of -178.  A kernel of the grid's own length would
% truncate that and lose a percent or two of the pattern.
koff = ((0:2*Mg-2)-(Mg-1))*fstep;                   % +-(Mg-1)*fstep

gauss_component=0; lorentz_component=0;
for cc = 1:3
    sel = (cls==cc);
    if ~any(sel), continue; end
    binned = bin_to_grid(twotheta(sel), F_squared(sel), tg, fstep);
    if Lorentzian_factor<1
        c_g = FWHMs(cc)/(2*(2*log(2))^0.5);
        gk  = exp(-koff.^2/(2*c_g^2));
        gauss_component = gauss_component + conv(binned, gk, 'same');
    end
    if Lorentzian_factor>0
        lk = 1./(koff.^2+(0.5*FWHMs(cc))^2);
        lorentz_component = lorentz_component + conv(binned, lk, 'same');
    end
end
% Cut the plotted range out of the working grid
i0 = round((exp_twotheta(1)-tg(1))/fstep)+1;
sl = i0:refine:(i0+refine*(numel(exp_twotheta)-1));
if Lorentzian_factor<1
    gauss_component = gauss_component(sl);
    gauss_part = gauss_component/max(gauss_component);
else
    gauss_part = 0;
end
if Lorentzian_factor>0
    lorentz_component = lorentz_component(sl);
    lorentz_part = lorentz_component/max(lorentz_component);
else
    lorentz_part = 0;
end

intensity=scalefactor*(Lorentzian_factor*lorentz_part+(1-Lorentzian_factor)*gauss_part);

%% Divergence slit
if Div_slit>0
    DIV=Sample_length/(Gonio_radius*Div_slit*pi()/180).*sin(exp_twotheta./2*pi()/180);
    DIV(DIV>1)=1;
else
    DIV=[sin(exp_twotheta/2*pi()/180)];
end
%% Surface roughness
SR=0.5*(1+(sin((exp_twotheta/2-roughness)*pi()/180)./sin((exp_twotheta/2+roughness)*pi()/180)));

if mode==1
    S_bar=((S1/2)^2+(S2/2)^2)^.5;
    Q=S_bar./(2*2^0.5*sin(exp_twotheta/2*pi()/180)*sigma_star);
    PSI=erf(Q)*(2*pi())^.5/(2*sigma_star*S_bar)-2*sin(exp_twotheta/2*pi()/180)/S_bar^2.*(1-exp(-Q.^2));
    Lorentz=(1+cos(exp_twotheta*pi()/180).^2);
    SingXtalLorentz=sin(exp_twotheta/2*pi()/180);
    if strcmp(L_type,'Reynolds')
        LP=Lorentz./ ( sin(exp_twotheta*pi()/180).*(sin(exp_twotheta/2*pi()/180)).^0.8);
    else
        LP=Lorentz./SingXtalLorentz.*PSI;
    end
    if RNDPWD == 1
        LP_random = Lorentz./(sin(exp_twotheta/2*pi()/180)) * 1./sin(exp_twotheta*pi()/180);
        LP=LP_random;
    end
    intensity=SR.*DIV.*LP.*intensity;
    intensity=real(intensity/max(intensity));
else
    if monochromator==0
        intensity=SR.*DIV.*intensity.*(1+cos(exp_twotheta*pi/180).^2)./(2*sin(exp_twotheta/2*pi/180).*sin(exp_twotheta*pi/180));
    else
        intensity=SR.*DIV.*intensity.*(1+cos(exp_twotheta*pi/180).^2)*(cos(monochromator_angle*pi/180).^2)./(cos(exp_twotheta/2*pi/180).*sin(exp_twotheta/2*pi/180).^2);
    end
    intensity=real(intensity/max(intensity));
end

assignin('caller','atom_xrd',atom)
assignin('caller','F_squared',F_squared)
assignin('caller','twotheta_disc',two_theta_disc)
assignin('caller','intensity',intensity)
assignin('caller','twotheta',exp_twotheta)
assignin('caller','h',h); assignin('caller','k',k); assignin('caller','l',l);
assignin('caller','hkl',hkl); assignin('caller','d_hkl',d_hkl);

if write_file
    writematrix([exp_twotheta' 100*intensity'],'xrd.dat','Delimiter','tab');
end

if nargin<5
    hold on;
    plot(exp_twotheta,intensity,'LineWidth',1);
    [peaks_int,locs_twotheta]=findpeaks(intensity,exp_twotheta,'MinPeakProminence',.05*max(intensity));
    if numel(peaks_int)<10
        [peaks_int,locs_twotheta]=findpeaks(intensity,exp_twotheta,'MinPeakProminence',.01*max(intensity));
    end
    if numel(peaks_int)<10
        [peaks_int,locs_twotheta]=findpeaks(intensity,exp_twotheta,'MinPeakProminence',.001*max(intensity));
    end
    assignin('caller','peaks_int',peaks_int)
    assignin('caller','locs_twotheta',locs_twotheta)
    intensity_disc=interp1(exp_twotheta,intensity,two_theta_disc);
    [~,ind_Intensity]=maxk(intensity_disc./max(intensity_disc),20*numel(peaks_int));
    two_theta_disc_Intensity_max=two_theta_disc(ind_Intensity);
    hkl_max_Intensity=hkl(ind_Intensity,:);
    hkl_abs_sorted=sort(abs(hkl),2,'descend');
    assignin('caller','hkl_abs_sorted',hkl_abs_sorted);
    hkl_ind=[];
    for i=1:numel(locs_twotheta)
        [dd, ind] = min(abs(two_theta_disc_Intensity_max-locs_twotheta(i)));
        if dd<1
            Miller_index=num2str(abs(hkl_max_Intensity(ind,:)));
            seq=sort(abs(hkl_max_Intensity(ind,:)),2,'descend');
            multiplicity=numel(find(ismember(hkl_abs_sorted,seq,'rows')));
            text(two_theta_disc_Intensity_max(ind)-0.32,peaks_int(i)+0.06,strcat('(',Miller_index(~isspace(Miller_index)),')'),'FontSize',14);
            if size(atom,2)<100
                text(two_theta_disc_Intensity_max(ind)-0.32,peaks_int(i)+0.12,num2str(multiplicity),'FontSize',14);
            end
            hkl_ind=[hkl_ind i];
        end
    end
    stem(locs_twotheta(hkl_ind),peaks_int(hkl_ind),'Color','black','MarkerEdgeColor','none');
    stem(locs_twotheta(hkl_ind),-0.03*ones(numel(locs_twotheta(hkl_ind))),'Color','black','MarkerEdgeColor','none');
    stem(two_theta_disc_Intensity_max,-0.03*ones(numel(two_theta_disc_Intensity_max)),'Color','black','MarkerEdgeColor','none');
    xlim([0 max(exp_twotheta)]);
    try, ylim([-.1 max(intensity)*1.15]); catch, end
    set(gca,'LineWidth',2,'FontName','Arial','FontSize',22);
    xlabel('Two-theta','FontSize',24);
    ylabel('Norm. intensity','FontSize',24);
end
end

function s = bin_to_grid(pos, wt, tg, step)
%% Accumulate weights onto a grid, splitting each between its two
%% neighbouring bins so a peak centre off the grid is not snapped onto it.
s = zeros(1,numel(tg));
u = (pos - tg(1))/step + 1;
u = u(:); wt = wt(:);
ok = isfinite(u) & isfinite(wt) & u>=1 & u<=numel(tg)-1;
u = u(ok); wt = wt(ok);
if isempty(u), return; end
i0 = floor(u); fr = u - i0;
s = s + accumarray(i0,  wt.*(1-fr), [numel(tg) 1]).';
s = s + accumarray(i0+1,wt.*fr,     [numel(tg) 1]).';
end
