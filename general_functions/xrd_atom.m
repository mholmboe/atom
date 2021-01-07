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
% * The script was inpired by MOF-FIT:
% * http://www.rsc.org/suppdata/ee/c3/c3ee40876k/c3ee40876k.pdf
%
%% Version
% 2.082
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # [twotheta,intensity] = xrd_atom(filename)
% # [twotheta,intensity] = xrd_atom(atom,Box_dim)
% # [twotheta,intensity] = xrd_atom(atom,Box_dim,[6 4 3]) % If your system has been replicated as 6x4x3
%
function [exp_twotheta,intensity] = xrd_atom(varargin)

%% Various settings
% num_hkl=72; % Maximum number of reflections, used for all h,k,l's, or edit manually later on..
lambda=1.54187; % Ångstrom
anglestep=0.02; % The incremental twotheta angle step
exp_twotheta=5:anglestep:80; % The twotheta range of interest
B_all=2; % Debye-Waller factor Ångstrom, in case no such field exist within the atom struct
Lorentzian_factor=1; % [0-1] Enter the fraction of the calculated pattern you would like to have described by a lorentzian function vs. a gaussian function
neutral_atoms=0;

%% Set FWHM
FWHM_00l=.25; % Specify the full width at half maximum of your choice
FWHM_hk0=.25; % Specify the full width at half maximum of your choice
FWHM_hkl=.25; % Specify the full width at half maximum of your choice

%% Various settings
mode=0; % Activate sigma_star, DIV, surface roughness as in Moore&Reynolds, 1997
Sample_length = 4; % cm
Gonio_radius = 24; % cm
Div_slit = .1; % Divergence slit setting, 0 for automatic
roughness = 0; % Surface roughness
sigma_star=45; % Reynolds 00l mean preferred orientation, 45 [deg] is random, 1 [deg] is the opposite
RNDPWD = 45; %  Random powder
% mu_star=45; % Not yet implemented
% alfa_strain = 0; % Not yet implemented

% misalignment=0.0; % Not yet implemented

L_type='normal'; % Lorent polarization type, else 'Reynolds';
S1=2.3;S2=2.3; % Primary and secondary Soller slit , in deg

%% Enter the degree (number greater than or equal to 0) and direction of preferential orientation
pref=0;
preferred_h=1;
preferred_k=1;
preferred_l=1;

%% Fetch either a .pdb|.gro file or use an atom struct with its Box_dim
if nargin==1
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
    Box_dim=varargin{2};
end

if nargin>2
    rep_factors=varargin{3};
else
    pause(2)
    rep_factors=[1 1 1];
end

if nargin>3
    selected_indexes=varargin{4};
else
    selected_indexes=[];
end

%% Specials section
% atom = unreplicate_atom(atom,Box_dim,rep_factors);
% for R=1:1
%     % %% Unreplicate the atom struct
%     % atom = unreplicate_atom(atom,Box_dim,rep_factors);
% 
%     %% Replicate and displace the atom struct
%     atom = replicate_atom(atom,Box_dim,[1 1 2]);
%     atom(size(atom,2)/2+1:end) = translate_atom(atom(size(atom,2)/2+1:end),[0 Box_dim(2)/3 0]);
%     atom = replicate_atom(atom,Box_dim,[2 2 1]);
%     rep_factors=rep_factors.*[2 2 2];
%     % plot_atom(atom,Box_dim);
%     % pause;
% 
%     %% Rotate some layers
%     new = replicate_atom(atom,Box_dim,[1 1 2]); % Generates a new Box_dim
%     rot = rotate_atom(new(size(new,2)/2+1:end),Box_dim,[0 0 2]);
%     rot = translate_atom(rot,[2 5 0]);
%     atom = update_atom({atom rot});
%     new = replicate_atom(atom,Box_dim,[1 1 2]); % Generates a new Box_dim
%     rot = rotate_atom(new(size(new,2)/2+1:end),Box_dim,[0 0 2]);
%     rot = translate_atom(rot,[2 3 0]);
%     atom = update_atom({atom rot});
%     rep_factors=rep_factors.*[1 1 4];
% end

% %% Slice the atom struct
% atom = slice_triclinic_atom(atom,Box_dim);
%
%% Wrap the structure into an atom struct
% atom=wrap_atom(atom,Box_dim);
%
% vmd(atom,Box_dim)
%
% pause

%% Enter the maximum h, k, and l values you would like to calculate. The calculated peaks will be for -hmax <= h <= hmax; -kmax<= k <= kmax; -lmax <= l <= lmax
% hmax=max([num_hkl rep_factors(1)*num_hkl]);
% kmax=max([num_hkl rep_factors(2)*num_hkl]);
% lmax=max([num_hkl rep_factors(3)*num_hkl]);
Cell=Box_dim2Cell(Box_dim);
hmax=ceil(exp_twotheta(end)/Bragg(lambda,'distance',Cell(1)));
kmax=ceil(exp_twotheta(end)/Bragg(lambda,'distance',Cell(2)));
lmax=ceil(exp_twotheta(end)/Bragg(lambda,'distance',Cell(3)));

%% Set the occupancy of all sites
if ~isfield(atom,'occupancy')
    try
        atom = occupancy_atom(atom,Box_dim);
    catch
        [atom.occupancy]=deal(1);
    end
end

occupancy=[atom.occupancy]';


%% /Specials section

%% Set the unit cell parameters
if size(Box_dim,2) == 9
    lx=Box_dim(1);    ly=Box_dim(2);    lz=Box_dim(3);
    xy=Box_dim(6);    xz=Box_dim(8);    yz=Box_dim(9);
elseif size(Box_dim,2) == 3
    lx=Box_dim(1);    ly=Box_dim(2);    lz=Box_dim(3);
    xy=0;    xz=0;    yz=0;
end

%%
frac=orto_atom(atom,Box_dim);
% frac = round_atom(frac,Box_dim,3,'orto');
frac=element_atom(frac);
atom_type=[frac.type];% atom_type(2:length(atom_type));
x=[frac.xfrac]';y=[frac.yfrac]';z=[frac.zfrac]';

if ~isfield(frac,'B')
    [frac.B]=deal(B_all);
end
Bvalue=[frac.B]';

%% Special section to set specific B-factors
% n=1;
% while (n<=length(atom_type))
%     if n>108
%         Bvalue(n)=11;
%         occupancy(n)=0.5;
%     end
%     n=n+1;
% end

a=lx;
b=(ly^2+xy^2)^.5;
c=(lz^2+xz^2+yz^2)^.5;
alfa=rad2deg(acos((ly*yz+xy*xz)/(b*c)));
beta=rad2deg(acos(xz/c));
gamma=rad2deg(acos(xy/b));

%% Converting to radians
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

hkl=[h_values_temp' k_values_temp' l_values_temp'];

%% Specials section

%% To remove nonsense indexes, if using a replicated system
if sum(abs(rep_factors-[1 1 1]))>0
    ind_rm=[];
    i=1;
    while i<size(hkl,1)+1
        if hkl(i,1)<rep_factors(1) && hkl(i,2)<rep_factors(2) && hkl(i,3)<rep_factors(3)
            ind_rm=[ind_rm i];
        end
        i=i+1;
    end
    hkl(unique(ind_rm),:)=[];
else
    disp('Is you system a single unit cell?')
    disp('Will assume no replication factors...')
end

%% To select only the specific indexes
hkl_selected=hkl;
if sum(abs(rep_factors-[1 1 1]))>0
    hkl_selected=hkl_selected./rep_factors;
end
if numel(selected_indexes)>0
    selected_ind=[];
    for i=1:size(selected_indexes,1)
        [row,col]=ismember(selected_indexes(i,:),hkl_selected,'rows');
        if col>0
            selected_ind=[selected_ind col];
        end
    end
    hkl=hkl(selected_ind,:);
end

% %% If calculating the XRD pattern from a supercell replicated as rep_factor=[x y z]
% %% Better to call the function with Box_dim./[rep_factors]
% if sum(abs(rep_factors-[1 1 1]))>0
%     ind_rm=[];
%     i=1;
%     while i<size(hkl,1)+1
%         if ~(mod(hkl(i,1),rep_factors(1))==0 && mod(hkl(i,2),rep_factors(2))==0 && mod(hkl(i,3),rep_factors(3))==0)
%             ind_rm=[ind_rm i];
%         end
%         i=i+1;
%     end
%     hkl(unique(ind_rm),:)=[];
% end

% ind_rm=[];
% for i=1:size(hkl,1)
%     %% Smectite special
%     if mod(hkl(i,3),10)==0 %(mod(hkl(i,1),rep_factors(1))==0 && mod(hkl(i,2),rep_factors(2))==0 &&  )
%         ind_rm=[ind_rm i];
%     end
% end
% hkl(unique(ind_rm),:)=[];

%% /Specials section
hkl = unique(hkl,'rows','first'); % extra, in case removing certain reflections.

h=hkl(:,1)';
k=hkl(:,2)';
l=hkl(:,3)';

%% Now we have all the h,k,l values with 0,0,0 taken out
%% Determine the two theta values for each h,k,l
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
[d_hkl,hkl_order]=sort(d_hkl,'descend');
twotheta_rad=twotheta_rad(hkl_order);
hkl=hkl(hkl_order,:);
h=h(hkl_order);
k=k(hkl_order);
l=l(hkl_order);

%% Calculate the corresponding discrete twotheta in degrees
two_theta_disc=real(180*twotheta_rad/pi);

%% Calculate the structure factor for each reflection
Atom_labels=unique(atom_type);
structure_factor=0;
m=1;
disp('--------------')
while m<numel(Atom_labels)+1
    ind=find(ismember(atom_type,Atom_labels(m)));
    if neutral_atoms==1
        Atom_labels(m)=strcat(Atom_labels(m),'0');
    end
    scattering_factor = atomic_scattering_factors(Atom_labels(m),lambda,two_theta_disc,Bvalue(ind(1)));
    numatoms=sum([atom(ind).occupancy])
    disp('--------------')
    assignin('caller',strcat(char(Atom_labels(m)),'_f')',scattering_factor);
    assignin('caller',strcat(char(Atom_labels(m)),'_Atomtype')',Atomtype);
    assignin('caller',strcat(char(Atom_labels(m)),'_nElectrons')',nElectrons);
    n=1;
    while n<numel(ind)+1
        structure_factor=structure_factor+scattering_factor.*occupancy(ind(n)).*exp(2*pi*1i.*(h*x(ind(n))+k*y(ind(n))+l*z(ind(n))));
        n=n+1;
    end
    m=m+1;
end

%% Run in parallell?
% for m=1:numel(Atom_labels)
%     ind=find(ismember(atom_type,Atom_labels(m)));
%     scattering_factor = atomic_scattering_factors(Atom_labels(m),lambda,two_theta_disc,Bvalue(ind(1)));
%     parfor n=1:numel(ind) % parallell for loop
%         structure_factor=structure_factor+scattering_factor.*occupancy(ind(n)).*exp(2*pi*1i.*(h*x(ind(n))+k*y(ind(n))+l*z(ind(n))));
%         %% Vectorized way of doing it. Seems slower... gives lower intensity?
%         % structure_factor=structure_factor+sum(occupancy(ind)*scattering_factor,1).*sum(exp(2*pi*1i.*(h.*x(ind)+k.*y(ind)+l.*z(ind))),1);
%     end
% end

%% Square each term in the structure factor vector by its complex conjugate
F_squared=structure_factor.*conj(structure_factor);

if sum(abs(rep_factors-[1 1 1]))>0
    hkl=hkl./rep_factors;
    h=h./rep_factors(1);
    k=k./rep_factors(2);
    l=l./rep_factors(3);
end

%% Correction for preferred orientation angle between h,k,l and preferrential orientation direction
theta_pref_orient=acos((h*preferred_h + k*preferred_k + l*preferred_l)./((h.^2+k.^2+l.^2).^0.5*(preferred_h^2+preferred_k^2+preferred_l^2)^0.5));
if theta_pref_orient>pi/2
    theta_pref_orient=pi-theta_pref_orient;
end
F_squared=F_squared.*exp(pref*cos(2*theta_pref_orient));
twotheta=real(180*twotheta_rad/pi);

%% Construct the calculated pxrd pattern by adding a lorentzian fraction to a gaussian fraction
if Lorentzian_factor<1
    % Gaussian part
    n=1;gauss_component=0;
    while(n<=length(twotheta))
        if hkl(n,3)==0
            temp_FWHM=FWHM_hk0;
        else
            if sum(hkl(n,1:2))==0
                temp_FWHM=FWHM_00l;
            else
                temp_FWHM=FWHM_hkl;
            end
        end
        %         temp_FWHM=temp_FWHM*(1+2*sin(twotheta(n)*pi/180));
        %         temp_FWHM=FWHM_hkl;
        
        c_g=temp_FWHM/(2*(2*log(2))^0.5);
        temp_gauss=F_squared(n).*exp(-(exp_twotheta-twotheta(n)).^2/(2*c_g^2));
        gauss_component=gauss_component+temp_gauss;
        n=n+1;
    end
    gauss_part=gauss_component/max(gauss_component);
else
    gauss_part=0;
end

if Lorentzian_factor>0
    % Lorentzian part
    n=1;lorentz_component=0;
    while(n<length(twotheta))
        if hkl(n,3)==0
            temp_FWHM=FWHM_hk0;
        else
            if sum(hkl(n,1:2))==0
                temp_FWHM=FWHM_00l;
            else
                temp_FWHM=FWHM_hkl;
            end
        end
        %         temp_FWHM=temp_FWHM*(1+2*sin(twotheta(n)*pi/180));
        %         temp_FWHM=FWHM_hkl;
        temp_lorentz=F_squared(n)./((exp_twotheta-twotheta(n)).^2+(0.5*temp_FWHM)^2);
        lorentz_component=lorentz_component+temp_lorentz;
        n=n+1;
    end
    lorentz_part=lorentz_component/max(lorentz_component);
else
    lorentz_part=0;
end

%% Mix together the Lorentzian and Gaussian parts in the ratio specified by eta to generate the pseudo-Voigt function
intensity=Lorentzian_factor*lorentz_part+(1-Lorentzian_factor)*gauss_part;

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
    
    if RNDPWD == 1
        LP_random = Lorentz./(sin(exp_twotheta/2*pi()/180)) * 1./sin(exp_twotheta*pi()/180);
        LP=LP_random;
    end
    intensity=SR.*DIV.*LP.*intensity;
    %     intensity=intensity-min(intensity);
    intensity=real(intensity/max(intensity));
    %     intensity=real(intensity/max(intensity(1:floor(15/((exp_twotheta(end)-exp_twotheta(1))/length(exp_twotheta))))));
else
    intensity=SR.*DIV.*intensity.*(1+cos(exp_twotheta*pi/180).^2)./(8*sin(exp_twotheta*pi/180/2).^2.*cos(exp_twotheta*pi/180/2));
    intensity=real(intensity/max(intensity));
    %     intensity=real(intensity/max(intensity(1:floor(15/((exp_twotheta(end)-exp_twotheta(1))/length(exp_twotheta))))));
end

assignin('caller','F_squared',F_squared)
assignin('caller','twotheta_rad',twotheta_rad)
assignin('caller','twotheta_disc',two_theta_disc)
assignin('caller','intensity',intensity)
assignin('caller','twotheta',exp_twotheta)
assignin('caller','h',h);
assignin('caller','k',k);
assignin('caller','l',l);
assignin('caller','hkl',hkl);
assignin('caller','d_hkl',d_hkl);

dlmwrite('xrd.dat',[exp_twotheta' 100*intensity'],'delimiter','\t','precision',5)

%% Plot the results
hold on;
%plot(exp_twotheta,intensity,'Color',[0 0 0],'LineWidth',1);
plot(exp_twotheta,intensity,'LineWidth',1);
% plot(exp_twotheta,intensity+(rand(2,length(intensity))-.5)/500,'k');

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
[peaks_Intensity,ind_Intensity]=maxk(intensity_disc./max(intensity_disc),20*numel(peaks_int));
two_theta_disc_Intensity_max=two_theta_disc(ind_Intensity);
hkl_max_Intensity=hkl(ind_Intensity,:);

[peaks_Fsq,ind_Fsq]=maxk(F_squared./max(F_squared),20*numel(peaks_int));
two_theta_disc_Fsq_max=two_theta_disc(ind_Fsq);
hkl_max_Fsq=hkl(ind_Fsq,:);

%% Miller indexes wrt the peak prominence
hkl_ind=[];
hkl_abs=abs(hkl);
hkl_abs_sorted=sort(hkl_abs,2,'descend');
assignin('caller','hkl_abs_sorted',hkl_abs_sorted);
for i=1:numel(locs_twotheta)
    [diff, ind] = min(abs(two_theta_disc_Intensity_max-locs_twotheta(i)));
    if diff<1
        
        Miller_index=num2str(abs(hkl_max_Intensity(ind,:)));
        Miller_seq=abs(hkl_max_Intensity(ind,:));
        seq=sort(abs(Miller_seq),2,'descend');
        multiplicity =numel(find(ismember(hkl_abs_sorted,seq,'rows')));
% remove        text(two_theta_disc_Intensity_max(ind)-3.2,peaks_int(i)+0.03,strcat('(',Miller_index(~isspace(Miller_index)),')'),'FontSize',14);
% remove        if size(atom,2)<100
% remove            text(two_theta_disc_Intensity_max(ind)-3.2,peaks_int(i)+0.09,num2str(multiplicity),'FontSize',14);
% remove        end
        hkl_ind=[hkl_ind i];
        
    end
end
% remove stem(locs_twotheta(hkl_ind),peaks_int(hkl_ind),'Color','black','MarkerEdgeColor','none');
% remove stem(locs_twotheta(hkl_ind),-0.03*ones(numel(locs_twotheta(hkl_ind))),'Color','black','MarkerEdgeColor','none');
% stem(two_theta_disc_Intensity_max,-0.03*ones(numel(two_theta_disc_Intensity_max)),'Color','black','MarkerEdgeColor','none');

xlim([0 max(exp_twotheta)]);
try
    ylim([-.4 max(intensity)*1.15])
catch
end

set(gca,'LineWidth',2,'FontName', 'Arial','FontSize',22);% ,'Xtick',exp_twotheta(1):10:exp_twotheta(end));%,'Xtick',...
xlabel('Two-theta','FontSize',24);
ylabel('Norm. intensity','FontSize',24);

% figure
% hold on;
% plot(exp_twotheta,intensity,'LineWidth',1);
% %% Miller indexes wrt the intensity (note LP factors included)
% for i=1:numel(ind_Intensity)
%     Miller_index=num2str(abs(hkl_max_Intensity(i,:)));
%     text(two_theta_disc_Intensity_max(i)-0.5,peaks_Intensity(i)+0.075,Miller_index(~isspace(Miller_index)));
% end
% stem(two_theta_disc_Intensity_max,peaks_Intensity,'Color','black','MarkerEdgeColor','none');
%
% xlim([0 max(exp_twotheta)]);
% try
%     ylim([-.2 max(intensity)*1.2])
% catch
% end
%
% figure
% hold on;
% plot(exp_twotheta,intensity,'LineWidth',1);
% %% Miller indexes wrt Fsq intensity (note no LP factors included)
% for i=1:numel(ind_Fsq)
%     Miller_index=num2str(abs(hkl_max_Fsq(i,:)));
%     text(two_theta_disc_Fsq_max(i)-0.5,peaks_Fsq(i)+0.075,Miller_index(~isspace(Miller_index)));
% end
% stem(two_theta_disc_Fsq_max,peaks_Fsq,'Color','black','MarkerEdgeColor','none');
% % stem(two_theta_disc(two_theta_disc<exp_twotheta(end)),-.05*ones(numel(two_theta_disc(two_theta_disc<exp_twotheta(end))),1),'Color','black','MarkerEdgeColor','none');
% % text(locs+.05,pks,num2str((1:numel(pks))'));
% xlim([0 max(exp_twotheta)]);
% try
%     %     ylim([-.2 max(intensity)*1.2])
% catch
% end

% assignin('caller','peaks_int',peaks_int)
% assignin('caller','locs_twotheta',locs_twotheta)

