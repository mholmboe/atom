close all;

%%
atomtypes1='OW';
ffstrings1={'opc3' 'opc3' 'opc3'};
atomtypes2={'Al'};
ffstrings2={'Merz_IOD_polyvalent_OPC3_ff' 'Merz_HFE_polyvalent_OPC3_ff' 'Merz_12_6_4_polyvalent_OPC3_ff'};

% for j=1:numel(ffstrings1)
for i=1:numel(atomtypes2)
    %     figure
    hold on;
    for j=1:numel(ffstrings2)
        atomtype2=atomtypes2{i};
        
        ffstring=ffstrings2{j}
        try
            atomtype1=strcat(atomtypes1,'_',ffstrings1{j})
        catch
            atomtype1=strcat(atomtypes1,'_',ffstrings1{end})
        end
        
        ff1=load(strcat('ions_',ffstring,'.mat'));
        ff2=load('water_models.mat');
        
        [r,lj,coul,Utot] = nonbonded_ff([ff1 ff2],{atomtype2 atomtype1},1);
    end
end
% end
%%

% tiledlayout(3,2)
% hold on;
% ff=load('clayff.mat');
% atomtypes1={'Al' 'Alt' 'Alt' 'Mgo' 'Si' 'Ob'};
% atomtypes2={'Oh'};
% for i=1:numel(atomtypes1)
%     nexttile(i)
% %     [r,lj,coul,Utot] = nonbonded_ff({ff.ff ff.ff},{atomtypes1{i} atomtypes2{1}},1);
% end
%
% %%