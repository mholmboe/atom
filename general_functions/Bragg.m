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