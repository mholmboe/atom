%% hist_atom.m
% * This function is used to calculate density profiles along X|Y|Z
%
%% Version
% 2.03
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # Hist = hist_atom(atom,Box_dim,0.02)
%
function [Binsx,Countsx,Binsy,Countsy,Binsz,Countsz] = hist_atom(atom,Box_dim,s,varargin)

Binsx = (0:s:Box_dim(1,1)+s)';
Binsy = (0:s:Box_dim(1,2)+s)';
Binsz = (0:s:Box_dim(1,3)+s)';

Coords=[atom.x; atom.y; atom.z;]';

if nargin>3
    center_vec=varargin{1};
    if numel(center_vec)==1
        Coords=Coords-center_vec;
    else
        Coords(:,1)=Coords(:,1)-center_vec(1);
        Coords(:,2)=Coords(:,2)-center_vec(2);
        Coords(:,3)=Coords(:,3)-center_vec(3);
    end
    Coords(Coords(:,1)<0,1)=Coords(Coords(:,1)<0,1)+Box_dim(1,1);
    Coords(Coords(:,2)<0,2)=Coords(Coords(:,2)<0,3)+Box_dim(1,2);
    Coords(Coords(:,3)<0,3)=Coords(Coords(:,3)<0,3)+Box_dim(1,3);
end

Countsx=histcounts(Coords(:,1),Binsx');
Countsy=histcounts(Coords(:,2),Binsy');
Countsz=histcounts(Coords(:,3),Binsz');
Countsx=Countsx/(s*Box_dim(1,2)*Box_dim(1,3)*1E-30*6.022e23*1000); % mol/L
Countsy=Countsy/(s*Box_dim(1,1)*Box_dim(1,3)*1E-30*6.022e23*1000); % mol/L
Countsz=Countsz/(s*Box_dim(1,1)*Box_dim(1,2)*1E-30*6.022e23*1000); % mol/L

if nargin>4 
    if varargin{2}>0
    Countsx=(Countsx+flipud(Countsx))/2;
    Countsy=(Countsy+flipud(Countsy))/2;
    Countsz=(Countsz+flipud(Countsz))/2;
    end
end

if nargin>5
    
    sigma = varargin{3};
    if sigma>0
    window = 100;
    x = linspace(-window / 2, window / 2, window);
    gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
    gaussFilter = gaussFilter / sum(gaussFilter); % normalize
    Countsx(:)=conv(Countsx(:)', gaussFilter, 'same');
    Countsy(:)=conv(Countsy(:)', gaussFilter, 'same');
    Countsz(:)=conv(Countsz(:)', gaussFilter, 'same');
    end
end

Binsx=Binsx(1:end-1)';
Binsy=Binsy(1:end-1)';
Binsz=Binsz(1:end-1)';
assignin('caller','DensityX',[Binsx,Countsx]);
assignin('caller','DensityY',[Binsy,Countsy]);
assignin('caller','DensityZ',[Binsz,Countsz]);


