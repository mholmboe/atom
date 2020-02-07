format long

%% Constants
eVinJ=1.60217653E-19; % SI
Na=6.022140857E+23;
kBinJK=1.3806485279E-23;

for i=1:size(ff,2)
    ff(i).radius_nm=ff(i).sigma_nm*2^(1/6);
%    ff(i).radius_nm=ff(i).radius_A/10;
    ff(i).radius_A=ff(i).radius_nm*10;
%     ff(i).sigma_nm=ff(i).radius_nm/2^(1/6);
    ff(i).sigma_A=ff(i).radius_A/2^(1/6);
    
    ff(i).e_kcalmol=ff(i).e_kJmol/4.184;
%     ff(i).e_kJmol=ff(i).e_kcalmol*4.184;
    ff(i).e_eV=ff(i).e_kJmol*1000/Na/eVinJ;
    ff(i).e_kB=ff(i).e_kJmol*1000/Na/kBinJK;
end

%% Order the records
[atomtypes,atomtypes_order]=sort([ff.type]);
ff_temp=ff;
for i=1:size(ff,2)
   ff_temp(i)=ff(atomtypes_order(i)); 
end
ff=ff_temp;

[atnum,atnum_order]=sort([ff.atnum]);
ff_temp=ff;
for i=1:size(ff,2)
   ff_temp(i)=ff(atnum_order(i)); 
end
ff=ff_temp;

%% Order the fields
defaultAttributes={'type' 'atnum' 'mass' 'charge' 'radius_nm' 'radius_A' 'sigma_nm' 'sigma_A' 'e_kJmol' 'e_kcalmol' 'e_eV' 'e_kB'}
ffAttributes=fieldnames(ff)';
indDefault=find(ismember(defaultAttributes,ffAttributes));
defaultAttributes=defaultAttributes(indDefault);
ind_ff=find(ismember(ffAttributes,defaultAttributes));
ffAttributes=ffAttributes(ind_ff);
ff=orderfields(ff,unique({defaultAttributes{:} ffAttributes{:}},'stable'));

% save('ions_tip4p','ff')