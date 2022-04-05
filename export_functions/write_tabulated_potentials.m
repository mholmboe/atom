% function ff = write_tabulated_potentials(r,q,c6,c12)

format long

Ion1='Na+'
Ion2='Na+'
filename_out=strcat('table_',Ion1,'_',Ion2,'.xvg');
ff1='ions_Merz_12_6_4_monovalent_opc3_ff.mat';
ff2='ions_Merz_12_6_4_monovalent_opc3_ff.mat'; %'water_models.mat'
% ff2='water_models.mat'
% Water_model='opc3';
% Ion2=strcat(Ion2,'_',Water_model)

s=0.002;
r=0:s:3; % nm

if strcmp(ff1,ff2)
    ff=load(ff1);
    all_ff=ff;
else
    ff1=load(ff1);
    ff2=load(ff2);
    all_ff=[ff ff2];
end

C4=0;
try
    [rout,lj,coul,Utot,q1,q2,sig1,sig2,eps1,eps2,C4] = nonbonded_ff(all_ff,{Ion1 Ion2});
    c4_mix=(C4*C4)^.5;
catch
    [rout,lj,coul,Utot,q1,q2,sig1,sig2,eps1,eps2] = nonbonded_ff(all_ff,{Ion1 Ion2});
end

% coul=(1.60217646E-19)^2*6.022E+23*q1*q2./(r*1E-9)*1/(4*3.14159*8.85E-12)/1000;coul(1)=coul(2)-(coul(3)-coul(2));
coul=1./r;coul(1)=coul(2)-(coul(3)-coul(2));

e_mix=(eps1*eps2)^.5;
sig_mix=(sig1+sig2)/2;
lj_sigeps=4*e_mix.*((sig_mix./r).^12-(sig_mix./r).^6)-c4_mix./r.^4;

c12_mix=4*e_mix*sig_mix^12;
c6_mix=4*e_mix*sig_mix^6;

lj_c6=-c6_mix./r.^6-c4_mix./r.^4;lj_c6(1)=lj_c6(2)-(lj_c6(3)-lj_c6(2));
lj_c12=c12_mix./r.^12;lj_c12(1)=lj_c12(2)-(lj_c12(3)-lj_c12(2));

Utot=lj_c6+lj_c12+coul;


% hold on
% plot(Data(:,1),Data(:,7))
%
% plot(r+.002,dlj_c12)

dcoul=coul./r;dcoul(1)=0;
dlj_c6=6*lj_c6./r;dlj_c6(1)=-0;
dlj_c12=12*lj_c12./r;dlj_c12(1)=0;

% dcoul(1)=0;
% dlj_c6(1)=0;
% dlj_c12(1)=0;
% 
% dcoul(2)=2*dcoul(3)-dcoul(4);
% dlj_c6(2)=2*dlj_c6(3)-dlj_c6(4);
% dlj_c12(2)=2*dlj_c12(3)-dlj_c12(4);
% 
% dcoul(end)=2*dcoul(end-1)-dcoul(end-2);
% dlj_c6(end)=2*dlj_c6(end-1)-dlj_c6(end-2);
% dlj_c12(end)=2*dlj_c12(end-1)-dlj_c12(end-2);

fid = fopen(filename_out, 'wt');
fprintf(fid, '%s\r\n','# Created in Matlab by MHolmboe');

i=0;
while i<numel(r)%-1
    i=i+1;
    atomtypes = [r(i) coul(i) dcoul(i) lj_c6(i) dlj_c6(i) lj_c12(i) dlj_c12(i) ];%q dq c6 dc6 c12 dc12];
    fprintf(fid, '%15.10e   %15.10e %15.10e   %15.10e %15.10e   %15.10e %15.10e\r\n', atomtypes);
end

fclose(fid);


% % %% Start test case
% % import_xvg('table6-12.xvg') % From the Gromacs dist
% % r=Data(:,1);
% % coul=Data(:,2);
% % dcoul=coul./r;
% % 
% % lj_c6=Data(:,4);
% % dlj_c6=6*lj_c6./r;
% % 
% % lj_c12=Data(:,6);
% % dlj_c12=12*lj_c12./r;
% % 
% % fid = fopen('test.xvg', 'wt');
% % fprintf(fid, '%s\r\n','# Created in Matlab by MHolmboe');
% % 
% % i=0;
% % while i<numel(r)%-1
% %     i=i+1;
% %     atomtypes = [r(i) coul(i) dcoul(i) lj_c6(i) dlj_c6(i) lj_c12(i) dlj_c12(i) ];%q dq c6 dc6 c12 dc12];
% %     fprintf(fid, '%15.10e   %15.10e %15.10e   %15.10e %15.10e   %15.10e %15.10e\r\n', atomtypes);
% % end
% % 
% % fclose(fid);
% % %% End test case



