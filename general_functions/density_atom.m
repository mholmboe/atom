%% density_atom.m
% * This function is used to calculate density profiles along X|Y|Z. It
% will also try to calc. the electron density profile from the an
% halfly-oxidized/ionic lattice using the X-ray scattering factors from
% atomic_scattering_factors(). This means its not exacly the same thing as
% the electron density, rather the electron density 'seen' by X-rays from
% an 'ionic lattice'.
% If the atom struct has a field/attribute atom.charge, this function also
% calculates the charge density, the electric field (symmetrized) and the
% electrostatic potential.
%
% Be(a)ware of the units...
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # [Bins,Element_density] = density_atom(atom,Box_dim)
% # [Bins,Element_density] = density_atom(atom,Box_dim,0.05)
% # [Bins,Element_density] = density_atom(atom,Box_dim,0.05,'z')
% # [Bins,Element_density] = density_atom(atom,Box_dim,0.05,'z',2.3) % sigma [>0] value for Gaussian smoothing
% # [Bins,Element_density] = density_atom(atom,Box_dim,0.05,'z',1,1) % mirror/symmetrize [0/1]
% # [Bins,Element_density] = density_atom(atom,Box_dim,0.05,'z',1,1,4.57) % Arbitrary value to shift all coords before mirroring etc

function [Bins,Total_Density] = density_atom(atom,Box_dim,varargin)

Fontsize=18;
q=1.602176634E-19; % 1 charge eq is q Coloumb (C)
epsilon_null = 8.854187817E-12; % A2·s4·kg-1·m-3
epsilon=1; % Dielectric constant
Na=6.02214086E23; % Avogadros constant

Atom_label = unique([atom.type]); % Cell containing strings with the atomtypes you want to include in the analysis
% If you want to exclude some atomtypes...
% Atoms_include     = unique([atom.type]);
% Atoms_exclude     = unique([atom(strcmp([atom.resname],'LAC')).type]);
% Atom_label        = setdiff(Atoms_include,Atoms_exclude);

if size(Box_dim,2)>6
    if find(Box_dim(1,6:end))>0
        atom = orto_atom(atom,Box_dim);
        Box_dim=orto_Box_dim;
    end
end

Limits=Box_dim(1:3);
if size(Limits,2)==3
    Limits(4)=Limits(1);
    Limits(5)=Limits(2);
    Limits(6)=Limits(3);
    Limits(1)=0;
    Limits(2)=0;
    Limits(3)=0;
end

if nargin<3
    ds=0.02;
else
    ds=varargin{1};
end

if nargin<4
    dimension='z';
else
    dimension=varargin{2};
end

Distance=0:ds:Box_dim(3);

if strcmpi(dimension,'x')
    Distance=0:ds:Box_dim(1);
    Area=(Limits(5)-Limits(2))*(Limits(6)-Limits(3));
elseif strcmpi(dimension,'y')
    Distance=0:ds:Box_dim(2);
    Area=(Limits(4)-Limits(1))*(Limits(6)-Limits(3));
elseif strcmpi(dimension,'z')
    Distance=0:ds:Box_dim(3);
    Area=(Limits(4)-Limits(1))*(Limits(5)-Limits(2));
end

DistanceMax=Distance(end);
Bins = (0:ds:Distance(end)+ds)';
nBins = numel(Bins);
assignin('base','Distance',Bins);

%% Set the occupancy of all sites
occ=1;
if ~isfield(atom,'occupancy')
    occ=1;
    try
        atom = occupancy_atom(atom,Box_dim);
    catch
        [atom.occupancy]=deal(1);
    end
end
Occupancy=[atom.occupancy];

Total_Density=zeros(ceil(Distance(end)/ds),1);
Total_Density_SI=Total_Density;
Total_Electron_density=Total_Density;
Total_Electron_density_SI=Total_Density;
Total_Charge_density=Total_Density;
Total_Charge_density_SI=Total_Density;
Total_EField=Total_Density;
Total_EField_SI=Total_Density;
Total_EPot=Total_Density;
Total_EPot_SI=Total_Density;
Element_count=0;
for h=1:length(Atom_label)
    Element_count = Element_count+1;
    % Extract the trajectory data along the selected dimension dim
    ind_atom = find(strncmpi([atom.type],Atom_label(h),3))';
    if size(ind_atom,1)<1
        ind_atom = find(strncmpi([atom.type],Atom_label(h),2))';
        if size(ind_atom,1)<1
            ind_atom = find(strcmpi([atom.type],Atom_label(h)))';
        end
    end
    
    if strcmpi(dimension,'x')
        Coords=[atom(ind_atom).x]';
    elseif strcmpi(dimension,'y')
        Coords=[atom(ind_atom).y]';
    elseif strcmpi(dimension,'z')
        Coords=[atom(ind_atom).z]';
    end
    Occupancy=[atom(ind_atom).occupancy];
    
    if nargin>6
        center_vec=varargin{5};
        Coords=Coords-center_vec;
    end
    
    Coords(Coords<0) = Coords(Coords<0) + Distance(end);
    Coords(Coords>Distance(end)) = Coords(Coords>Distance(end)) - Distance(end);
    
%     Element_density=histcounts(Coords(:,1),Bins(1:end-1))';
    Element_density=zeros(length(Bins),1);
    for i=1:numel(Coords)
        Coords(i)
        bin_ind=floor((Coords(i)/DistanceMax)*nBins+1)
        Element_density(bin_ind)=Element_density(bin_ind)+Occupancy(i);
    end
    
    if numel(Element_density)~=numel(Total_Density)
        Element_density=interp1(1:numel(Element_density),Element_density,1:numel(Total_Density),'spline')';
    end
    
    Element_density=Element_density/(ds*Area*1E-30*6.022e23*1000); % mol/L
    
    if nargin>4
        sigma = varargin{3};
        if sigma>0
            window = 100;
            x = linspace(-window / 2, window / 2, window);
            gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
            gaussFilter = gaussFilter / sum(gaussFilter); % normalize
            Element_density(:)=conv(Element_density(:)', gaussFilter, 'same');
        end
    end
    
    if nargin>5
        if varargin{4}>0
            Element_density=(Element_density+flipud(Element_density))/2;
        end
    end
    
    Element_density_SI=Element_density*1000*Na; % npart/m^3.
    
    nElectrons = atomic_scattering_factors(Atom_label(h),1.542,0,0); % From X-ray scattering tables, hence not identical to electron density
    Element_electron_density=nElectrons*Element_density; % mol eq/L
    Element_electron_density_SI=q*nElectrons*Element_density_SI; % C/m^3
    
    if isfield(atom,'charge')
        if sum(unique([atom(ind_atom).charge]))>1
            disp('Averaging the atom type charge for')
            Atom_label(h)
            unique([atom(ind_atom).charge])
        end
        Element_Charge_density = Na*1E-24*Element_density*mean([atom(ind_atom).charge]); % q_ev/nm^3
        Element_EField = cumsum((Element_Charge_density + flipud(Element_Charge_density))/2)*ds*1E-10*1E27*q/(1E9*epsilon_null*epsilon);
        Element_EPot = -cumsum(Element_EField)*ds*1e-10*1E9; %
        
        Element_Charge_density_SI = q*Element_density_SI*mean([atom(ind_atom).charge]); % npart*A*s / m3  == npart*C / m3
        Element_EField_SI = cumsum((Element_Charge_density_SI + flipud(Element_Charge_density_SI))/2)*ds*1e-10/(epsilon_null*epsilon); % 1e-10m/(epsilon_null*epsilon)
        Element_EPot_SI = -cumsum(Element_EField_SI)*ds*1e-10; % *1e-10m; % V=kg*m2/(A*s2)
    end
    
    disp('Atom_label')
    Atom_label(h)
    size(ind_atom)
    sum(Element_density)
    disp('-----------')
    
    % Here we keep track of the atoms we analyze
    UC(Element_count).type=char(Atom_label(h));
    UC(Element_count).N=sum(Element_density);
    
    % Collect all the data
    Total_Density=Total_Density+Element_density;
    Total_Electron_density=Total_Electron_density+Element_electron_density;
    Total_Density_SI=Total_Density_SI+Element_density_SI;
    if isfield(atom,'charge')
        Total_Charge_density=Total_Charge_density+Element_Charge_density;
        Total_EField=Total_EField+Element_EField;
        Total_EPot=Total_EPot+Element_EPot;
        
        Total_Charge_density_SI=Total_Charge_density_SI+Element_Charge_density_SI;
        Total_EField_SI=Total_EField_SI+Element_EField_SI;
        Total_EPot_SI=Total_EPot_SI+Element_EPot_SI;
    end
    
    % Assign the data to atomtype specific variables
    if sum(Element_density) > 0
        assignin('base',strcat(char(Atom_label(h)),'_density'),[Element_density Element_density_SI]);
        assignin('base',strcat(char(Atom_label(h)),'_density'),[Element_electron_density Element_electron_density_SI]);
        if isfield(atom,'charge')
            assignin('base',strcat(char(Atom_label(h)),'_charge_density'),[Element_Charge_density Element_Charge_density_SI]);
            assignin('base',strcat(char(Atom_label(h)),'_efield'),[Element_EField Element_EField_SI]);
            assignin('base',strcat(char(Atom_label(h)),'_epot'),[Element_EPot Element_EPot_SI]);
        end
        hold on
        plot(Distance(1:length(Element_density)),Element_density)
        xlabel('Distance [Å]','FontSize',Fontsize); ylabel('mol/L','FontSize',Fontsize);
    end
end

atom = mass_atom(atom,Box_dim);
assignin('caller','Mw',Mw);
assignin('caller','Box_volume',Box_volume);
assignin('caller','Box_density',Box_density);

hold on;
plot(Distance(1:length(Total_Density)),circshift(Total_Density,[0 0]),'k')
xlabel('Distance [Å]','FontSize',Fontsize); ylabel('Concentration profile [mol/L]','FontSize',Fontsize);
legend([Atom_label {'Total'}]);
hold off

figure
plot(Distance(1:length(Total_Density)),circshift(Total_Electron_density,[0 0]),'r')
xlabel('Distance [Å]','FontSize',Fontsize); ylabel('Concentration profile [mol eqv/L]','FontSize',Fontsize);

if isfield(atom,'charge')
    figure
    plot(Distance(1:length(Total_Density)),circshift(Total_Charge_density,[0 0]))
    xlabel('Distance [Å]','FontSize',Fontsize); ylabel('Charge density [q/nm^3]','FontSize',Fontsize);
    figure
    plot(Distance(1:length(Total_Density)),circshift(Total_EField,[0 0]))
    xlabel('Distance [Å]','FontSize',Fontsize); ylabel('Electric field [V/nm]','FontSize',Fontsize);
    figure
    plot(Distance(1:length(Total_Density)),circshift(Total_EPot,[0 0]))
    xlabel('Distance [Å]','FontSize',Fontsize); ylabel('Electrostatic potential [V]','FontSize',Fontsize);
    figure
    plot(Distance(1:length(Total_Density)),circshift(Total_Charge_density_SI,[0 0]))
    xlabel('Distance [Å]','FontSize',Fontsize); ylabel('Charge density [C/m^3]','FontSize',Fontsize);
    figure
    plot(Distance(1:length(Total_Density)),circshift(Total_EField_SI,[0 0]))
    xlabel('Distance [Å]','FontSize',Fontsize); ylabel('Electric field [V/m]','FontSize',Fontsize);
    figure
    plot(Distance(1:length(Total_Density)),circshift(Total_EPot_SI,[0 0]))
    xlabel('Distance [Å]','FontSize',Fontsize); ylabel('Electrostatic potential [V]','FontSize',Fontsize);
end
% if strcmpi(Unit,'conc')
% else
%     xlabel('Distance [Å]','FontSize',Fontsize); ylabel('Num','FontSize',Fontsize);
% end

% save(strcat('density_profs_',filename))

