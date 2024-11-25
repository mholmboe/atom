%% density_atom_fast.m
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
% # [Bins,Total_Density,Total_Density_SI] = density_atom_fast(atom,Box_dim)
% # [Bins,Total_Density,Total_Density_SI] = density_atom_fast(atom,Box_dim,0.05)
% # [Bins,Total_Density,Total_Density_SI] = density_atom_fast(atom,Box_dim,0.05,'z')
% # [Bins,Total_Density,Total_Density_SI] = density_atom_fast(atom,Box_dim,0.05,'z',2.3) % sigma [>0] value for Gaussian smoothing
% # [Bins,Total_Density,Total_Density_SI] = density_atom_fast(atom,Box_dim,0.05,'z',1,1) % mirror/symmetrize [0/1]
% # [Bins,Total_Density,Total_Density_SI] = density_atom_fast(atom,Box_dim,0.05,'z',1,1,4.57) % Arbitrary value to shift all coords before mirroring etc

function [Bins,Total_Density,Total_Density_SI] = density_atom_fast(atom,Box_dim,varargin)

% Fontsize=18;
q=1.602176634E-19; % 1 charge eq is q Coloumb (C)
% epsilon_null = 8.854187817E-12; % A2·s4·kg-1·m-3
% epsilon=1; % Dielectric constant
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
        Coords(i);
        bin_ind=floor((Coords(i)/DistanceMax)*nBins+1);
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
    
%     % Collect all the data
    Total_Density=Total_Density+Element_density;
    Total_Density_SI=Total_Density_SI+Element_density_SI;

    % Assign the data to atomtype specific variables
    if sum(Element_density) > 0
        assignin('caller',strcat(char(Atom_label(h)),'_density'),[Element_density Element_density_SI]);
    end
end


atom = mass_atom(atom,Box_dim);
assignin('caller','Mw',Mw);
assignin('caller','Distance',Distance);
assignin('caller','Box_volume',Box_volume);
assignin('caller','Box_density',Box_density);

