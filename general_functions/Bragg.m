%% Examples
% # d = Bragg(1.5406,'twotheta',6.2) % d-spacing [same unit as lambda] from a 2theta angle [deg]
% # d = Bragg(1.5406,'theta',3.1) % The same, from a theta angle [deg]
% # twotheta = Bragg(1.5406,'d',14.2) % 2theta [deg] from a d-spacing
% # d = Bragg(1.5406,'twotheta',6.2,2) % The same, for the second order reflection
%
%% See also
% # <xrd_atom.html xrd_atom>
% # <atomic_scattering_factors.html atomic_scattering_factors>
%
function new_var = Bragg(lambda,string,var,varargin)

% Order of reflection
if nargin>3
   n=varargin{1}; 
else
    n=1;
end

if strcmp(string,'twotheta')
    new_var=n*lambda/sin(var/2*pi()/180)/2;
elseif strcmp(string,'theta')
    new_var=n*lambda/sin(var*pi()/180)/2;
else
    new_var=2*asin((n*lambda)/(2*var))*180/pi();
end