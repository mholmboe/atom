function new_var = Bragg(lambda,string,var)

if strcmp(string,'twotheta')
    new_var=lambda/sin(var/2*pi()/180)/2;
elseif strcmp(string,'theta')
    new_var=lambda/sin(var*pi()/180)/2;
else
    new_var=2*asin(lambda/(2*var))*180/pi();
end