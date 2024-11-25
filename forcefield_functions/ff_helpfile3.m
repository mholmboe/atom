format long

% ffname='ions_Zhang'
% watermodels={'opc3' 'tip3p' 'tip3pfb' 'spce' 'spceb' 'opc' 'tip4pew'  'tip4p2005' 'tip4pfb' 'tip4pd' 'a99SB-disp'}

ffname='ions_Merz_HFE'
watermodels={'opc3' 'tip3p' 'tip3pfb' 'spce' 'spceb' 'opc' 'tip4pew'  'tip4p2005' 'tip4pfb' 'tip4pd' 'a99SB-disp'}

for n=1:numel(watermodels)
    
    watermodel=char(watermodels(n));
    
    try
    ff1=load(strcat(ffname,'_monovalent_',watermodel,'_ff.mat'))
    ff2=load(strcat(ffname,'_divalent_',watermodel,'_ff.mat'))
    ff3=load(strcat(ffname,'_polyvalent_',watermodel,'_ff.mat'))
    
    ff=[ff1.ff ff2.ff ff3.ff];
    [atnum,atnum_order]=sort([ff.atnum]);
    ff_temp=ff;
    for i=1:size(ff,2)
        ff_temp(i)=ff(atnum_order(i));
    end
    ff=ff_temp;
    
    save(strcat(ffname,'_',watermodel,'_ff.mat'),'ff')
    catch
       disp('Did not find') 
       watermodel
       pause
    end
end