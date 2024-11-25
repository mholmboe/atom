%% write_mdp.m
%% This function writes a mdp file from an imported mdp struct.

%% Written by MHolmboe
%% Please report bugs/issues to michael.holmboe@umu.se

%% Version
% 3.00

function write_mdp(mdp,varargin)

if nargin>1
    mdp_filename_out=varargin{1};
else
    mdp_filename_out='out.mdp';
end

%% Create vars for the sections
param = fieldnames(mdp);

fid = fopen(mdp_filename_out, 'wt');

fprintf(fid, '%s % s\r\n',';','mdp file written in MATLAB. Email bugs to michael.holmboe@umu.se');

if exist('mdp','var')
    for i = 1:size(param,1)
        temp_string=getfield(mdp,param{i});
        if size(temp_string,1) > 1
            if ~isnumeric(temp_string(1,:))
                str_temp=[];
                for j=1:size(temp_string,1)
                    str_temp=strcat(str_temp,temp_string(j,:),({' '}));
                end
                fprintf(fid, '%-23s %-1s %-s \r\n',char(param(i)),'=',char(str_temp));
            elseif isnumeric(temp_string)
                str_temp=[];
                for j=1:size(temp_string,1)
                    str_temp=strcat(str_temp,num2str(temp_string(j,:),'%.2f'),({' '}));
                end
                fprintf(fid, '%-23s %-1s %-4s\r\n',char(param(i)),'=',char(str_temp));
            end
        elseif isnumeric(temp_string)
            fprintf(fid, '%-23s %-1s %-s\r\n',char(param(i)),'=',num2str(temp_string));
        elseif ischar(temp_string)>0 && length(temp_string) > 0
            fprintf(fid, '%-23s %-1s %-s\r\n',char(param(i)),'=',char(temp_string));
        end
    end
end

fprintf(fid, '\r\n');

fclose(fid);

end