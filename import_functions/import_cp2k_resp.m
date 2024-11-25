%% import_cp2k_resp.m
% * This function imports resp charges calculated in/from cp2k. Not very well tested...
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # resp_charges = import_cp2k_resp('molecule.resp')
%
function resp_charges = import_cp2k_resp(varargin)

if nargin>0
    filename=varargin{1};
else
    filename='MIN-RESP_CHARGES.resp';
end

fid = fopen(filename,'r');
data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', ''); % New addition
data=data{1}; % New addition
fclose(fid);
Total = strfind(data,'Total');
Index = find(not(cellfun('isempty',Total)));

resp_charges=[];
for i = 4:Index-1
    line = split(data{i});
    resp_charges=[resp_charges; str2double(line(end))];
end


if nargin>1
    atom=import_atom(varargin{2});
    for i=1:size(atom,2)
        atom(i).charge=resp_charges(i);
    end
    write_atom_pqr(atom,Box_dim,'resp.pqr')
end


