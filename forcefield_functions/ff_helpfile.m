format long

ffname='Merz_12-6-4_divalent'

watermodels={'OPC3' 'OPC' 'tip3p_fb' 'tip4p_fb'}

for n=1:numel(watermodels)
    
    watermodel=char(watermodels(n))
    
    Ions=PolyvalentIons;
    for i=1:numel(Ions)
        eval(strcat('ff(i).type=Ions(i)'));
        eval(strcat('ff(i).radius_A=2*R2min_',watermodel,'(i)'));
        eval(strcat('ff(i).e_kcalmol=eps_',watermodel,'(i)'));
        
        try
            eval(strcat('ff(i).C4_kcalmolA=C4_',watermodel,'(i)'));
        catch
            
        end
        eval(strcat('ff(i).charge=2'));
    end
    
    ff_temp = mass_atom(ff); % Relies on element_atom...
    
    [ff.atnum]=ff_temp.atnum;
    [ff.mass]=ff_temp.mass;
    
    %% Constants
    eVinJ=1.60217653E-19; % SI
    Na=6.022140857E+23;
    kBinJK=1.3806485279E-23;
    
    for i=1:size(ff,2)
        %     ff(i).radius_nm=ff(i).sigma_nm*2^(1/6);
        ff(i).radius_nm=ff(i).radius_A/10;
        %    ff(i).radius_A=ff(i).radius_nm*10;
        ff(i).sigma_nm=ff(i).radius_nm/2^(1/6);
        ff(i).sigma_A=ff(i).radius_A/2^(1/6);
        
        %     ff(i).e_kcalmol=ff(i).e_kJmol/4.184;
        ff(i).e_kJmol=ff(i).e_kcalmol*4.184;
        ff(i).e_eV=ff(i).e_kJmol*1000/Na/eVinJ;
        ff(i).e_kB=ff(i).e_kJmol*1000/Na/kBinJK;
        
        try
            ff(i).C4_kJmolnm=ff(i).C4_kcalmolA*4.184/10000;
        catch
            
        end
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
    defaultAttributes={'type' 'atnum' 'mass' 'charge' 'radius_nm' 'radius_A' 'sigma_nm' 'sigma_A' 'e_kJmol' 'e_kcalmol' 'e_eV' 'e_kB' 'C4_kcalmolA' 'C4_kJmolnm' }
    ffAttributes=fieldnames(ff)';
    indDefault=find(ismember(defaultAttributes,ffAttributes));
    defaultAttributes=defaultAttributes(indDefault);
    ind_ff=find(ismember(ffAttributes,defaultAttributes));
    ffAttributes=ffAttributes(ind_ff);
    ff=orderfields(ff,unique({defaultAttributes{:} ffAttributes{:}},'stable'));
    
%      fff=load('ions_Merz_12-6-4_monovalent_spce_ff.mat')
%     
%     for i=1:size(ff,2)
%         ff(i).charge=fff.ff(i).charge;
%     end
    
     save(strcat('ions_',ffname,'_',watermodel,'_ff'),'ff')
    
end