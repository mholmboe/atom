close all;

atomtypes1={'Li' 'Na' 'K' 'Cs' 'Ca'};
ffstrings={'spc' 'spce' 'tip3p' 'opc3'};


for i=1:numel(atomtypes1)
%     figure
    hold on;
    for j=1:numel(ffstrings)
    atomtype1=atomtypes1{i};
    atomtype2='OW';
    
    ffstring=ffstrings{j};
    atomtype2=strcat(atomtype2,'_',ffstring);
    ff1=load(strcat('ions_',ffstring,'.mat'));
    ff2=load(strcat('water_models.mat'));
    
    [r,lj,coul,Utot] = ljcoul_ff([ff1.ff ff2.ff],{atomtype1 atomtype2});
    end
end