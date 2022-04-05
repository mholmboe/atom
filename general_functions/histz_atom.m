%% hist_atom.m
% * This function is used to calculate density profiles in the Z-direction.
% See also the density_atom function.
%
%% Version
% 2.11
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # Hist = histz_atom(atom,Box_dim,0.02)
%
function [Bins,Counts] = histz_atom(atom,Box_dim,s,varargin)

if nargin>3
    dimension=varargin{1};
else
    dimension=3;
end

if dimension==1
    Area=Box_dim(1,2)*Box_dim(1,3);
elseif dimension == 2
    Area=Box_dim(1,1)*Box_dim(1,3);
elseif dimension == 3
    Area=Box_dim(1,1)*Box_dim(1,2);
end

Bins = (0:s:Box_dim(1,dimension)+s)';
Coords=[atom.x; atom.y; atom.z;]';

if nargin>4
    center=varargin{2};
    Coords(:,dimension)=Coords(:,dimension)-center;
    Coords(Coords(:,dimension)<0,dimension)=Coords(Coords(:,dimension)<0,dimension)+Box_dim(1,dimension);
end

Counts=histcounts(Coords(:,dimension),numel(Bins))';
Counts=Counts/(s*Area*1E-30*6.022e23*1000); % mol/L

if nargin>5
    if varargin{3}>0
        Counts=(Counts+flipud(Counts))/2;
    end
end

if nargin>6
    sigma = varargin{4};
    window = 100;
    x = linspace(-window / 2, window / 2, window);
    gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
    gaussFilter = gaussFilter / sum(gaussFilter); % normalize
    Counts(:)=conv(Counts(:)', gaussFilter, 'same');
end

assignin('caller','Density',[Bins,Counts]);


