%% element_atom.m
% * This function replaces the atomtypes names with the element names
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * atom = element_atom(atom) 
% * atom = element_atom(atom,'Ow','Hw')

function atom = element_atom(atom,varargin) 

if nargin > 1;
    water_O=varargin{1}(:);
    water_H=varargin{2}(:);
else
    water_O='Ow';
    water_H='Hw';
end

if nargin > 3;
   [atom.type]=atom.fftype; 
end

for i=1:size(atom,2);
    if strncmpi(atom(i).type,{'Si'},2);atom(i).element={'Si'};
    elseif strncmpi(atom(i).type,{'SY'},2);atom(i).element={'Si'};
    elseif strncmpi(atom(i).type,{'SC'},2);atom(i).element={'Si'};
    elseif strncmpi(atom(i).type,{'S'},1);atom(i).element={'S'};
    elseif strncmpi(atom(i).type,{'Al'},2);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'AC'},2);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'AY'},2);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'Mg'},2);atom(i).element={'Mg'};
    elseif strncmpi(atom(i).type,{'Fe'},2);atom(i).element={'Fe'};
    elseif strncmpi(atom(i).type,{'Ow'},2);atom(i).element={water_O};
    elseif strncmpi(atom(i).type,{'Hw'},2);atom(i).element={water_H};
    elseif strncmpi(atom(i).type,{'Li'},2);atom(i).element={'Li'};
    elseif strncmpi(atom(i).type,{'Na'},2);atom(i).element={'Na'};
    elseif strncmpi(atom(i).type,{'K'},1);atom(i).element={'K'};
    elseif strncmpi(atom(i).type,{'Ca'},2);atom(i).element={'Ca'};
    elseif strncmpi(atom(i).type,{'O'},1);atom(i).element={'O'};
    elseif strncmpi(atom(i).type,{'H'},1);atom(i).element={'H'};
    elseif strncmpi(atom(i).type,{'Ca'},2);atom(i).element={'Ca'};
    elseif strncmpi(atom(i).type,{'Cl'},2);atom(i).element={'Cl'};
    elseif strncmpi(atom(i).type,{'C'},1);atom(i).element={'C'};
    else
        [atom(i).element]=atom(i).type;
    end
end

[atom.type]=atom.element;
