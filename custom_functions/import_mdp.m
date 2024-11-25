%% import_mdp.m
%% This function imports a mdp file.
%%
%% Written by MHolmboe
%% Please report bugs/issues to michael.holmboe@umu.se

%% Version
% 3.00

function mdp = import_mdp(mdp_filename)

inputfile = fopen(mdp_filename, 'r');
C = textscan(inputfile, '%s', 'Delimiter', '\n');
fclose(inputfile);
nondatarows1=find(strncmp(C{1},';',1));
nondatarows2=find(strncmp(C{1},{''},1));
C{1,1}([nondatarows1;nondatarows2])=[];

mdp_data=C{1,1};

for i=1:size(mdp_data,1)
    ind_equal=cell2mat(regexp(mdp_data(i),'='));
    ind_rm=cell2mat(regexp(mdp_data(i),';'));
    mdp_temp=char(mdp_data(i));
    mdp_param=mdp_temp(1:ind_equal-1);
    if numel(ind_rm)>0
        mdp_value=mdp_temp(ind_equal+1:ind_rm(1)-1);
    else
        mdp_value=mdp_temp(ind_equal+1:size(mdp_temp,2));
    end

    mdp_value=split(mdp_value);
    mdp_value=mdp_value(~cellfun('isempty',mdp_value));
    num_mdp_value=str2num(num2str(char(mdp_value)));
    if numel(num_mdp_value)>0
        mdp_value=num_mdp_value;
    else
        mdp_value=char(mdp_value);
    end
    field=char(mdp_param);
    field=strtrim(field);
    field=strrep(field,'-','_');
    mdp.(field)=mdp_value;
end

end