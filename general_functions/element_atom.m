%% element_atom.m
% * This function can replace the atomtypes names with the element names
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = element_atom(atom)
% # atom = element_atom(atom,'Ow','Hw')
%
function atom = element_atom(atom,varargin)

[atom.fftype]=atom.type;

if nargin > 1
    water_O=varargin{1}(:);
    water_H=varargin{2}(:);
else
    water_O='Ow';
    water_H='Hw';
end

if nargin > 3
    [atom.type]=atom.fftype;
end

for i=1:size(atom,2)
    
    %     [atom(i).type]=atom(i).type{1}(1:2); % Elements do not have more than two characters;
    
    if strncmpi(atom(i).type,{'Si'},2);atom(i).element={'Si'};
    elseif strncmpi(atom(i).type,{'Sr'},2);atom(i).element={'Sr'};
    elseif strncmpi(atom(i).type,{'SY'},2);atom(i).element={'Si'};
    elseif strncmpi(atom(i).type,{'SC'},2);atom(i).element={'Si'};
    elseif strncmpi(atom(i).type,{'st'},2);atom(i).element={'Si'};
    elseif strncmp(atom(i).type,{'s'},1);atom(i).element={'Si'};
    elseif strncmpi(atom(i).type,{'S'},1);atom(i).element={'S'};
    elseif strncmpi(atom(i).type,{'Al'},2);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'AC'},2);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'AY'},2);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'ao'},2);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'a'},1);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'Br'},2);atom(i).element={'Br'};
    elseif strncmpi(atom(i).type,{'B'},1);atom(i).element={'B'};
    elseif strncmpi(atom(i).type,{'I'},1);atom(i).element={'I'};
    elseif strncmpi(atom(i).type,{'Mg'},2);atom(i).element={'Mg'};
    elseif strncmpi(atom(i).type,{'Fe'},2);atom(i).element={'Fe'};
    elseif strncmpi(atom(i).type,{'F'},1);atom(i).element={'F'};
    elseif strncmpi(atom(i).type,{'U'},1);atom(i).element={'U'};
    elseif strncmpi(atom(i).type,{'V'},1);atom(i).element={'V'};
    elseif strncmpi(atom(i).type,{'Y'},1);atom(i).element={'Y'};
    elseif strncmpi(atom(i).type,{'Ow'},2);atom(i).element={water_O};
    elseif strncmpi(atom(i).type,{'Hw'},2);atom(i).element={water_H};
    elseif strncmpi(atom(i).type,{'Li'},2);atom(i).element={'Li'};
    elseif strncmpi(atom(i).type,{'Mn'},2);atom(i).element={'Mn'};
    elseif strncmpi(atom(i).type,{'Na'},2);atom(i).element={'Na'};
    elseif strncmpi(atom(i).type,{'Ni'},2);atom(i).element={'Ni'};
    elseif strncmpi(atom(i).type,{'Nh'},2);atom(i).element={'Nh'};
    elseif strncmpi(atom(i).type,{'Nb'},2);atom(i).element={'Nb'};
    elseif strncmpi(atom(i).type,{'Ne'},2);atom(i).element={'Ne'};
    elseif strncmpi(atom(i).type,{'No'},2);atom(i).element={'No'};
    elseif strncmpi(atom(i).type,{'N'},1);atom(i).element={'N'};
    elseif strncmpi(atom(i).type,{'K'},1);atom(i).element={'K'};
    elseif strncmpi(atom(i).type,{'O'},1);atom(i).element={'O'};
    elseif strncmpi(atom(i).type,{'o'},1);atom(i).element={'O'};
    elseif strncmpi(atom(i).type,{'H'},1);atom(i).element={'H'};
    elseif strncmpi(atom(i).type,{'h'},1);atom(i).element={'H'};
    elseif strncmpi(atom(i).type,{'Ca'},2);atom(i).element={'Ca'};
    elseif strncmpi(atom(i).type,{'Cl'},2);atom(i).element={'Cl'};
    elseif strncmpi(atom(i).type,{'C'},1);atom(i).element={'C'};
    else
        [atom(i).element{1}(1)]=upper(atom(i).type{1}(1));
        if size([atom(i).type{:}],2)>1
            i
            atom(1).type
            [atom(i).element{1}(2)]=lower(atom(1).type{1}(2));
        end
    end
end

[atom.type]=atom.element;
