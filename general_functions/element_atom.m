%% element_atom.m
% * This function can replace the atomtypes names with the element names
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
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
    if strncmpi(atom(i).type,{'Ag'},2);atom(i).element={'Ag'};
    elseif strncmpi(atom(i).type,{'Al'},2);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'AC'},2);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'AY'},2);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'at'},2);atom(i).element={'Alt'};
    elseif strncmpi(atom(i).type,{'ao'},2);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'a'},1);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'Br'},2);atom(i).element={'Br'};
    elseif strncmpi(atom(i).type,{'Ba'},2);atom(i).element={'Ba'};
    elseif strncmpi(atom(i).type,{'Be'},2);atom(i).element={'Be'};
    elseif strncmpi(atom(i).type,{'B'},1);atom(i).element={'B'};
    elseif strncmpi(atom(i).type,{'Ca'},2);atom(i).element={'Ca'};
    elseif strncmpi(atom(i).type,{'Ce'},2);atom(i).element={'Ce'};
    elseif strncmpi(atom(i).type,{'Cd'},2);atom(i).element={'Cd'};
    elseif strncmpi(atom(i).type,{'Co'},2);atom(i).element={'Co'};
    elseif strncmpi(atom(i).type,{'Cr'},2);atom(i).element={'Cr'};
    elseif strncmpi(atom(i).type,{'Cs'},2);atom(i).element={'Cs'};
    elseif strncmpi(atom(i).type,{'Cu'},2);atom(i).element={'Cu'};
    elseif strncmpi(atom(i).type,{'Cl'},2);atom(i).element={'Cl'};
    elseif strncmpi(atom(i).type,{'C'},1);atom(i).element={'C'};
    elseif strncmpi(atom(i).type,{'Dy'},2);atom(i).element={'Dy'};
    elseif strncmpi(atom(i).type,{'Eu'},2);atom(i).element={'Eu'};
    elseif strncmpi(atom(i).type,{'Er'},2);atom(i).element={'Er'};
    elseif strncmpi(atom(i).type,{'Fe'},2);atom(i).element={'Fe'};
    elseif strncmpi(atom(i).type,{'F'},1);atom(i).element={'F'};
    elseif strncmpi(atom(i).type,{'Gd'},2);atom(i).element={'Gd'};
    elseif strncmpi(atom(i).type,{'Hg'},2);atom(i).element={'Hg'};
    elseif strncmpi(atom(i).type,{'Hf'},2);atom(i).element={'Hf'};
    elseif strncmpi(atom(i).type,{'H'},1);atom(i).element={'H'};
    elseif strncmpi(atom(i).type,{'h'},1);atom(i).element={'H'};
    elseif strncmpi(atom(i).type,{'I'},1);atom(i).element={'I'};
    elseif strncmpi(atom(i).type,{'OW'},2);atom(i).element={water_O};
    elseif strncmpi(atom(i).type,{'Wa'},2);atom(i).element={water_O};
    elseif strncmpi(atom(i).type,{'OHH'},3);atom(i).element={water_O};
    elseif strncmpi(atom(i).type,{'Ow'},2);atom(i).element={water_O};
    elseif strncmpi(atom(i).type,{'Hw'},2);atom(i).element={water_H};
    elseif strncmpi(atom(i).type,{'In'},2);atom(i).element={'In'};
    elseif strncmpi(atom(i).type,{'K'},1);atom(i).element={'K'};
    elseif strncmpi(atom(i).type,{'Li'},2);atom(i).element={'Li'};
    elseif strncmpi(atom(i).type,{'La'},2);atom(i).element={'La'};
    elseif strncmpi(atom(i).type,{'Lu'},2);atom(i).element={'Lu'};
    elseif strncmpi(atom(i).type,{'Mg'},2);atom(i).element={'Mg'};
    elseif strncmpi(atom(i).type,{'Mn'},2);atom(i).element={'Mn'};
    elseif strncmpi(atom(i).type,{'Mo'},2);atom(i).element={'Mo'};
    elseif strncmpi(atom(i).type,{'Na'},2);atom(i).element={'Na'};
    elseif strncmpi(atom(i).type,{'Nd'},2);atom(i).element={'Nd'};
    elseif strncmpi(atom(i).type,{'Ni'},2);atom(i).element={'Ni'};
    elseif strncmpi(atom(i).type,{'Nh'},2);atom(i).element={'Nh'};
    elseif strncmpi(atom(i).type,{'Nb'},2);atom(i).element={'Nb'};
    elseif strncmpi(atom(i).type,{'Ne'},2);atom(i).element={'Ne'};
    elseif strncmpi(atom(i).type,{'No'},2);atom(i).element={'No'};
    elseif strncmpi(atom(i).type,{'N'},1);atom(i).element={'N'};
    elseif strncmpi(atom(i).type,{'O-H'},3);atom(i).element={'O'};
    elseif strncmpi(atom(i).type,{'Oh'},2);atom(i).element={'O'};
    elseif strncmpi(atom(i).type,{'O'},1);atom(i).element={'O'};
    elseif strncmpi(atom(i).type,{'o'},1);atom(i).element={'O'};
    elseif strncmpi(atom(i).type,{'W'},1);atom(i).element={'W'};
    elseif strncmpi(atom(i).type,{'Pb'},2);atom(i).element={'Pb'};
    elseif strncmpi(atom(i).type,{'Pr'},2);atom(i).element={'Pr'};
    elseif strncmpi(atom(i).type,{'Pu'},2);atom(i).element={'Pu'};
    elseif strncmpi(atom(i).type,{'Pd'},2);atom(i).element={'Pd'};
    elseif strncmpi(atom(i).type,{'Pt'},2);atom(i).element={'Pt'};
    elseif strncmpi(atom(i).type,{'P'},1);atom(i).element={'P'};
    elseif strncmpi(atom(i).type,{'RB'},2);atom(i).element={'Rb'};
    elseif strncmpi(atom(i).type,{'Si'},2);atom(i).element={'Si'};
    elseif strncmpi(atom(i).type,{'Sr'},2);atom(i).element={'Sr'};
    elseif strncmpi(atom(i).type,{'Sm'},2);atom(i).element={'Sm'};
    elseif strncmpi(atom(i).type,{'Sn'},2);atom(i).element={'Sn'};
    elseif strncmpi(atom(i).type,{'SY'},2);atom(i).element={'Si'};
    elseif strncmpi(atom(i).type,{'SC'},2);atom(i).element={'Si'};
    elseif strncmpi(atom(i).type,{'st'},2);atom(i).element={'Si'};
    elseif strncmp(atom(i).type,{'s'},1);atom(i).element={'Si'};
    elseif strcmp(atom(i).type,{'S'});atom(i).element={'S'};
    elseif strncmpi(atom(i).type,{'Tl'},2);atom(i).element={'Tl'};
    elseif strncmpi(atom(i).type,{'Ti'},2);atom(i).element={'Ti'};
    elseif strncmpi(atom(i).type,{'Tb'},2);atom(i).element={'Tb'};
    elseif strncmpi(atom(i).type,{'Th'},2);atom(i).element={'Th'};
    elseif strncmpi(atom(i).type,{'Tm'},2);atom(i).element={'Tm'};
    elseif strncmpi(atom(i).type,{'U'},1);atom(i).element={'U'};
    elseif strncmpi(atom(i).type,{'V'},1);atom(i).element={'V'};
    elseif strncmpi(atom(i).type,{'XX'},2);atom(i).element={'XX'};
    elseif strncmpi(atom(i).type,{'X'},2);atom(i).element={'X'};
    elseif strncmpi(atom(i).type,{'Yb'},2);atom(i).element={'Yb'};
    elseif strncmpi(atom(i).type,{'Y'},1);atom(i).element={'Y'};
    elseif strncmpi(atom(i).type,{'Zr'},2);atom(i).element={'Zr'};
    elseif strncmpi(atom(i).type,{'Zn'},2);atom(i).element={'Zn'};
        
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
