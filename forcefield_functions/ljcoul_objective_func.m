function res = ljcoul_objective_func(param,data,r,scalefactors,varargin)

param=param./scalefactors;

q1=param(1);
q2=param(2);

sig1=param(3);
sig2=param(4);

eps1=param(5);
eps2=param(6);

if numel(param)>6
    C41=param(7);
    C42=param(7);
end

if numel(param)>7
    C41=param(7);
    C42=param(8);
end

if numel(param)<7
    
    e_mix=(eps1*eps2)^.5;
    sig_mix=(sig1+sig2)/2;
    
    lj=4*e_mix.*((sig_mix./r).^12-(sig_mix./r).^6);
    coul=(1.60217646E-19)^2*6.022E+23*q1*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;
    
else
    
    C61=4*eps1*sig1^6;
    C121=4*eps1*sig1^12;
    C62=4*eps2*sig2^6;
    C122=4*eps2*sig2^12;
    
    C12_mix=(C121*C122)^.5;
    C6_mix=(C61*C62)^.5;
    C4_mix=(C41*C42)^.5;
    
    if numel(param)>7
        lj=C12_mix./r.^12-C6_mix./r.^6-C4_mix./r.^4;
    else
        lj=C12_mix./r.^12-C6_mix./r.^6;
    end
    coul=(1.60217646E-19)^2*6.022E+23*q1*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;
    
end

Utot=lj+coul;

res=data-Utot;
