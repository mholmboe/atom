%% rdf_atom.m
% * This function calculates the radial distribution function (RDF)
% and the coordination number, CN, between the atom1 (solute) and atom2
% (ligand/s), see Examples below. If passing a numeric fifth argunment, 
% smoothing thorugh Gaussian convolution will be invoked.
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Dependent functions
% * dist_matrix = dist_matrix_atom(atom1,atom2,Box_dim); % Calculates the
% distance matrix between all sites in atom1 and all sites in atom2
% * dist_matrix = cell_list_dist_matrix_atom(atom12,Box_dim,0,Distance(end)); % Used for systems with nAtoms>5000. Calculates the
% distance matrix between all sites in atom1+atom2.
%
%% Usage
% * Assuming an atom struct contains Na ions in water with Ow water oxygens
% first separate the Na solute struct (atom1) and the Ow ligands (atom2):
%
% # atom1 = atom(strcmpi([atom.type],'Na'))
% # atom2 = atom(strcmpi([atom.type],'Ow'))
%
% The run the rdf_atom() function as shown below.
%
%% Examples
% # [Distance,RDF,CN] = rdf_atom(atom1,atom2,Box_dim)               % atom1 is the central atomtype, atom2 is the ligand/s
% # [Distance,RDF,CN] = rdf_atom(atom1,atom2,Box_dim,[0:0.025:12])  % The full distance vector in Ångström
% # [Distance,RDF,CN] = rdf_atom(atom1,atom2,Box_dim,0.025)         % Here 0.025 is the binsize in Ångström
% # [Distance,RDF,CN] = rdf_atom(atom1,atom2,Box_dim,0.025,2)       % Here 2 is a smoothing factor

function [Distance,RDF,CN] = rdf_atom(atom1,atom2,Box_dim,varargin)

if nargin>3
    temp_var=varargin{1};
    if numel(temp_var)>1
        Distance=varargin{1}';
        binsize=Distance(2)-Distance(1);
    else
        binsize=temp_var;
    end
end

if ~exist('binsize','var')
    binsize=0.02;
end

if ~exist('Distance','var')
    Distance=[binsize/2:binsize:12+binsize/2]';
end

if Distance(end)>min([Box_dim(1:3)./2])
    disp('Note that the integration distance is larger than half the smallest system size')
end

V=4/3*pi()*Distance(end)^3;
Nsolute=size(atom1,2);
Nlig=size(atom2,2);

natom1=size(atom1,2);
natom2=size(atom2,2);
if (natom1+natom2)>50000 && numel(Box_dim)<9 % Will use the cell list routine to calc the reduced distance matrix
    atom12=[atom1 atom2];
    dist_matrix=cell_list_dist_matrix_atom(atom12,Box_dim,0,Distance(end));
    dist_matrix(natom1+1:end,:)=[];
    dist_matrix=dist_matrix(:,natom1+1:end);
    dist_matrix(dist_matrix==0)=Distance(end)+10; % Will set all distances larger than Distance(end) to some larger dummy value
else % Will calc the full distance matrix
    dist_matrix = dist_matrix_atom(atom1,atom2,Box_dim);
end

rdist=reshape(dist_matrix',1,[]);
rdist(rdist==0)=[];
rdist(rdist>Distance(end))=[];
Nlig=numel(rdist)/Nsolute;
Bins=histcounts(rdist,Distance)';

if nargin>4
    sigma = varargin{2};
    if sigma>0
        window = 100;
        x = linspace(-window / 2, window / 2, window);
        gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
        gaussFilter = gaussFilter / sum(gaussFilter); % normalize
        Bins(:)=conv(Bins(:)', gaussFilter, 'same');
    end
end

Distance=Distance(1:end-1);
RDF=Bins./(Nsolute*4*pi*Distance.^2*binsize*Nlig/V);
RDF(isnan(RDF))=0; % First index value of RDF may become NaN...
CN=cumsum(Bins./Nsolute);

% hold on
% plot(Distance,RDF)
% plot(Distance,CN)

