%% xrd_atom.m
% * This function calculates theoretical XRD patterns from a .pdb|.gro etc.
% * The script was heavily inpired by MOF-FIT:
% * http://www.rsc.org/suppdata/ee/c3/c3ee40876k/c3ee40876k.pdf
%
%% Version
% 2.0
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # [twotheta,intensity] = xrd_atom(filename)
% # [twotheta,intensity] = xrd_atom(atom,Box_dim))
%
function [twotheta,intensity] = xrd_atom(varargin)

%% Various settings
mode=1; % Activate sigma_star, DIV, surface roughness
FWHM=.8; % Specify the full width at half maximum of your choice
num_hkl=10;
lambda=1.54187; % Ångstrom
anglestep=0.02; % The incremental twotheta angle step
exp_twotheta=2:anglestep:70; % The twotheta range of interest
Bfactor=2; % Debye-Waller factor Ångstrom
Lorentzian_factor=1; % Enter the fraction of the calculated pattern you would like to have described by a lorentzian function vs. a gaussian function
Sample_length = 4; % cm
Gonio_radius = 24; % cm
Div_slit = .1; % Divergence slit setting, 0 for automatic
roughness = 0; % Surface roughness
% misalignment=0.0; % Not yet implemented
% alfa_strain = 0; % Not yet implemented
sigma_star=45; % Reynolds 00l mean preferred orientation, 45 [deg] is random, 1 [deg] is the opposite
RNDPWD = 0; % Random powder
% mu_star=45; % Not yet implemented
L_type='normal'; % Lorent polarization type, else 'Reynolds';
S1=2.3;S2=2.3; % Primary and secondary Soller slit , in deg

%% Fetch either a .pdb|.gro file or use an atom struct with its Box_dim
if nargin==1
    filename=varargin{1};
    atom=import_atom(filename);
else
    atom=varargin{1};
    Box_dim=varargin{2};
end
%% Import the structure into an atom struct
atom=wrap_atom(atom,Box_dim);
% atom=replicate_atom(atom,Box_dim,[4 4 4]); % If we want to replicate the struct
%%
frac=orto_atom(atom,Box_dim);
% frac = round_atom(frac,Box_dim,3,'orto');
frac=element_atom(frac);
atom_type=[frac.type];% atom_type(2:length(atom_type));
x=[frac.xfrac]';y=[frac.yfrac]';z=[frac.zfrac]';
if ~isfield(frac,'occupancy')
    [frac.occupancy]=deal(1);
end

occupancy=[frac.occupancy]';

%% Enter the maximum h, k, and l values you would like to calculate. The calculated peaks will be for -hmax <= h <= hmax; -kmax<= k <= kmax; -lmax <= l <= lmax
hmax=num_hkl;
kmax=num_hkl;
lmax=num_hkl;
%% Enter the degree (number greater than or equal to 0) and direction of preferential orientation
pref=1;
preferred_h=0;
preferred_k=0;
preferred_l=1;

%% Set the unit cell parameters
if size(Box_dim,2) == 9
    lx=Box_dim(1);    ly=Box_dim(2);    lz=Box_dim(3);
    xy=Box_dim(6);    xz=Box_dim(8);    yz=Box_dim(9);
elseif size(Box_dim,2) == 3
    lx=Box_dim(1);    ly=Box_dim(2);    lz=Box_dim(3);
    xy=0;    xz=0;    yz=0;
end

a=lx;
b=(ly^2+xy^2)^.5;
c=(lz^2+xz^2+yz^2)^.5;
alfa=rad2deg(acos((ly*yz+xy*xz)/(b*c)));
beta=rad2deg(acos(xz/c));
gamma=rad2deg(acos(xy/b));

%% Converting to radians
%%%%% From MOF-FIT %%%%%%
alfa_rad=alfa*pi/180;
beta_rad=beta*pi/180;
gamma_rad=gamma*pi/180;
%% Setting up the different h,k,l values
l_values_temp=repmat(-1*lmax:lmax,1,(2*kmax+1)*(2*hmax+1));
k_repeat_unit=repmat(-1*kmax:kmax,(2*lmax+1),1);
k_values_temp=repmat(reshape(k_repeat_unit,1,(2*lmax+1)*(2*kmax+1)),1,(2*hmax+1));
h_values_temp=reshape(repmat(-1*hmax:hmax,(2*kmax+1)*(2*lmax+1),1),1,(2*hmax+1)*(2*kmax+1)*(2*lmax+1));

h_values_temp(hmax*(2*kmax+1)*(2*lmax+1)+kmax*(2*lmax+1)+lmax+1)=[];
k_values_temp(hmax*(2*kmax+1)*(2*lmax+1)+kmax*(2*lmax+1)+lmax+1)=[];
l_values_temp(hmax*(2*kmax+1)*(2*lmax+1)+kmax*(2*lmax+1)+lmax+1)=[];

h=h_values_temp;
k=k_values_temp;
l=l_values_temp;
hkl=[h' k' l'];

%% Now we have all the h,k,l values with 0,0,0 taken out
%% Determine the twotheta values for each h,k,l
V_cell=a*b*c*(1-cos(alfa_rad)^2-cos(beta_rad)^2-cos(gamma_rad)^2+2*cos(alfa_rad)*cos(beta_rad)*cos(gamma_rad))^0.5;
one_over_dhkl=1/V_cell.*...
    (h.^2*b^2*c^2*sin(alfa_rad)^2+...
    k.^2*a^2*c^2*sin(beta_rad)^2+...
    l.^2*a^2*b^2*sin(gamma_rad)^2+...
    2*h.*k*a*b*c^2*(cos(alfa_rad)*cos(beta_rad)-cos(gamma_rad))+...
    2*k.*l*a^2*b*c*(cos(beta_rad)*cos(gamma_rad)-cos(alfa_rad))+...
    2*h.*l*a*b^2*c*(cos(alfa_rad)*cos(gamma_rad)-cos(beta_rad))).^(0.5);
one_over_dhkl=real(one_over_dhkl);
%%%%% End MOF-FIT %%%%%%
%      assignin('caller','one_over_dhkl1',one_over_dhkl);
%% Order the hkl after decreasing d-spacings and increasing two_theta
d_hkl=1./one_over_dhkl;
twotheta_rad=2.*asin(one_over_dhkl*lambda/2);
[d_hkl,hkl_order]=sort(d_hkl,'ascend');
twotheta_rad=twotheta_rad(hkl_order);
hkl=hkl(hkl_order,:);
assignin('caller','hkl',hkl);
assignin('caller','d_hkl',d_hkl);
h=h(hkl_order);
k=k(hkl_order);
l=l(hkl_order);

%% Calculate the corresponding discrete twotheta in degrees
two_theta_disc=real(180*twotheta_rad/pi);

%% Calculate the structure factor for each reflection
% Can we speed this up by fetching the scattering factor only once for each atom type?
n=1;
structure_factor=0;
while (n<=length(atom_type))
    scattering_factor = atomic_scattering_factors(atom_type(n),lambda,two_theta_disc,Bfactor);
    %         structure_factor=structure_factor+scattering_factor.*occupancy(n).*exp(2*pi*1i.*(h*x(n)+k*y(n)+l*z(n)));
    structure_factor=structure_factor+scattering_factor.*occupancy(n).*exp(2*pi*1i.*(l*z(n)));
    
    n=n+1;
end

%% Square each term in the structure factor vector
F_squared=structure_factor.*conj(structure_factor);

%% Correction for preferred orientation angle between h,k,l and preferrential orientation direction
theta_pref_orient=acos((h*preferred_h + k*preferred_k + l*preferred_l)./((h.^2+k.^2+l.^2).^0.5*(preferred_h^2+preferred_k^2+preferred_l^2)^0.5));
if theta_pref_orient>pi/2
    theta_pref_orient=pi-theta_pref_orient;
end
F_squared=F_squared.*exp(pref*cos(2*theta_pref_orient));

twotheta=real(180*twotheta_rad/pi);
%% Construct the calculated pxrd pattern by adding a lorentzian fraction to a gaussian fraction
% Gaussian part
n=1;gauss_component=0;c_=FWHM/(2*(2*log(2))^0.5);
while(n<=length(twotheta))
    temp_gauss=F_squared(n).*exp(-(exp_twotheta-twotheta(n)).^2/(2*c_^2));
    gauss_component=gauss_component+temp_gauss;
    n=n+1;
end
gauss_part=gauss_component/max(gauss_component);
% Lorentzian part
n=1;lorentz_component=0;
while(n<length(twotheta))
    temp_lorentz=F_squared(n)./((exp_twotheta-twotheta(n)).^2+(0.5*FWHM)^2);
    lorentz_component=lorentz_component+temp_lorentz;
    n=n+1;
end
lorentz_part=lorentz_component/max(lorentz_component);
% Mix together the Lorentzian and Gaussian parts in the ratio specified by eta to generate the pseudo-Voigt function
intensity=Lorentzian_factor*lorentz_part+(1-Lorentzian_factor)*gauss_part;

%% Lorentz factor
S_bar=((S1/2)^2+(S2/2)^2)^.5;
Q=S_bar./(2*2^0.5*sin(exp_twotheta/2*pi()/180)*sigma_star);
PSI=erf(Q)*(2*pi())^.5/(2*sigma_star*S_bar)-2*sin(exp_twotheta/2*pi()/180)/S_bar^2.*(1-exp(-Q.^2));

%% Lorentz * Polarization factors
Lorentz=(1+cos(exp_twotheta*pi()/180).^2);
SingXtalLorentz=sin(exp_twotheta/2*pi()/180);
if strcmp(L_type,'Reynolds')
    LP=Lorentz./ ( sin(exp_twotheta*pi()/180).*(sin(exp_twotheta/2*pi()/180)).^0.8);
else
    LP=Lorentz./SingXtalLorentz.*PSI;
end
LP_random = Lorentz./(sin(exp_twotheta/2*pi()/180)) * 1./sin(exp_twotheta*pi()/180);

if RNDPWD == 1
    LP=LP_random;
end

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
    intensity=SR.*DIV.*LP.*intensity;
    intensity=intensity-min(intensity);
    intensity=real(intensity/max(intensity(1:floor(15/((exp_twotheta(end)-exp_twotheta(1))/length(exp_twotheta))))));
else
    intensity=intensity.*(1+cos(exp_twotheta*pi/180).^2)./(8*sin(exp_twotheta*pi/180/2).^2.*cos(exp_twotheta*pi/180/2));
    %     intensity=real(intensity/max(intensity));
    %    intensity=real(intensity/max(intensity(1:floor(15/((exp_twotheta(end)-exp_twotheta(1))/length(exp_twotheta))))));
end

%% Plot the results
hold on;
plot(exp_twotheta,intensity,'b');


% plot(exp_twotheta,intensity+(rand(1,length(intensity))-.5)/100,'k');
% text(locs+.02,pks,num2str((1:numel(pks))'));
[pks,locs]=findpeaks(intensity,exp_twotheta);
% stem(locs,-.05*ones(numel(locs),1),'Color','black','MarkerEdgeColor','none');
stem(two_theta_disc(two_theta_disc<exp_twotheta(end)),-.05*ones(numel(two_theta_disc(two_theta_disc<exp_twotheta(end))),1),'Color','black','MarkerEdgeColor','none');
%
ind_diff=[1 find(diff(two_theta_disc(two_theta_disc<exp_twotheta(end)+5)))];
if length(ind_diff) < 20
    for i=2:length(ind_diff)
        %     temp_ind_matrix=hkl(ind_diff(i):ind_diff(i),:);
        %     [temp_max,temp_max_ind]=max(sum(temp_ind_matrix,2));
        %     ind_diff(i)
        %     temp_max_ind
        Miller_index=num2str(hkl(ind_diff(i),:));
        text(two_theta_disc(ind_diff(i))-1.2,-.08,Miller_index(~isspace(Miller_index)));
        ind_peaks=find(ismember(exp_twotheta,locs));
    end
end
xlim([0 max(exp_twotheta)])
try
    ylim([-.2 max(intensity)])
catch
end
