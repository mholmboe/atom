%% clayff_2004_param.m
% * This function holds the Clayff parameters + ion pair potentials from
% Joung and Cheatham, 2008 + some others...
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # clayff_param(Atom_label)
% # clayff_param(Atom_label,'spc')
%
function clayff_2004_param(Atom_label,varargin)

if nargin < 2
    water_model='SPC/E';
else
    water_model=varargin{1};
end

% Atom_label   = {'Hw', 'H','Ow','Oh','O', 'Omg', 'Oalt','Odsub','Ohmg','Oalsi','Oalhh','Oalh','Osih','Si','Al','Alt','Mgo','Mgh','Cao','Cah','Feo','Lio','Li','Na','K','Cs','Mg','Ca','Sr','Ba','Cl','Br'}';

% call function with find(ismember({forcefield.clayff.type},'Hw')) to fin Hw's index

%                      1     2    3    4   5     6       7      8       9     10      11      12     13    14   15   16    17    18    19    20    21    22    23   24  25   26   27   28   29   30   31   32  33  34   35 36
% From Cygan, 2004  = {'h*','ho','o*','oh','ob','obos','obts','obss', 'ohs', 'oas', 'oahhe','oahe', 'oshe','st','ao','Alt','mgo','mgh','cao','cah','feo','lio','Li','Na','K','Rb','Cs','Mg','Ca','Sr','Ba','F','Cl','Br','I'}';
% MHolmboe_CLAYFF   = {'Hw', 'H','Ow','Oh','O', 'Omg', 'Oalt','Odsub','Ohmg','Oalsi','Oalhh','Oalh','Osih','Si','Al','Alt','Mgo','Mgh','Cao','Cah','Feo','Lio','Li','Na','K','Rb','Cs','Mg','Ca','Sr','Ba','F','Cl','Br','I'}';

% Hw
forcefield.clayff(1).type =    {'Hw'}; % h* water hydrogen
forcefield.clayff(1).mass =    1.00794;
forcefield.clayff(1).charge =  0.42380;
if sum(strcmpi('SPC/E', water_model)) > 0 || sum(strcmpi('SPCE', water_model)) > 0
    forcefield.clayff(1).charge =  0.42380;
elseif sum(strcmpi('SPC', water_model)) > 0
    forcefield.clayff(1).charge =  0.41000;
elseif sum(strncmpi('TIP3P', water_model,5)) > 0
    forcefield.clayff(1).charge =  0.41700;
end
forcefield.clayff(1).radius =  0.00000; % R0 Å
forcefield.clayff(1).e_kcalmol = 0.00000; % kcal mol-1
%
% H
forcefield.clayff(2).type =    {'ho'}; % ho hydroxyl hydrogen
forcefield.clayff(2).mass =    1.00794;
forcefield.clayff(2).charge =  0.42500;
forcefield.clayff(2).radius =  0.00000;
forcefield.clayff(2).e_kcalmol = 0.00000;
%
% Ow
forcefield.clayff(3).type =    {'Ow'}; % o* water oxygen
forcefield.clayff(3).mass =    15.99410;
if sum(strcmpi('SPC/E', water_model)) > 0 || sum(strcmpi('SPCE', water_model)) > 0
    forcefield.clayff(3).charge =  -0.84760;
elseif sum(strcmpi('SPC', water_model)) > 0
    forcefield.clayff(3).charge = -0.82000;
elseif sum(strncmpi('TIP3P', water_model,5)) > 0
    forcefield.clayff(3).charge =  -0.83400;
end
forcefield.clayff(3).radius =  3.553200;
forcefield.clayff(3).e_kcalmol = 0.155400;
%
% Oh
forcefield.clayff(4).type =    {'oh'}; % oh hydroxyl oxygen
forcefield.clayff(4).mass =    15.99410;
forcefield.clayff(4).charge =  -0.95;
forcefield.clayff(4).radius =  3.553200;
forcefield.clayff(4).e_kcalmol = 0.155400;
%
% O
forcefield.clayff(5).type =    {'ob'}; % ob bridging oxygen
forcefield.clayff(5).mass =    15.99410;
forcefield.clayff(5).charge =  -1.05000;
forcefield.clayff(5).radius =  3.553200;
forcefield.clayff(5).e_kcalmol = 0.155400;
%
% Omg
forcefield.clayff(6).type =    {'obos'}; % ob bridging oxygen w. octahedral substitution
forcefield.clayff(6).mass =    15.99410;
forcefield.clayff(6).charge =  -1.18080;
forcefield.clayff(6).radius =  3.553200;
forcefield.clayff(6).e_kcalmol = 0.155400;
%
% Oalt
forcefield.clayff(7).type =    {'obts'}; % obts bridging oxygen w. tetrahedral substitution
forcefield.clayff(7).mass =    15.99410;
forcefield.clayff(7).charge =  -1.16880;%75;%-1.16880; Changed for ZSM-5
forcefield.clayff(7).radius =  3.553200;
forcefield.clayff(7).e_kcalmol = 0.155400;
%
% Odsub
forcefield.clayff(8).type =    {'obss'}; % obss bridging oxygen w. double substitution
forcefield.clayff(8).mass =    15.99410;
forcefield.clayff(8).charge =  -1.29960;
forcefield.clayff(8).radius =  3.553200;
forcefield.clayff(8).e_kcalmol = 0.155400;
%
% Ohmg
forcefield.clayff(9).type =    {'ohs'}; % ohs bridging oxygen w. substitution
forcefield.clayff(9).mass =    15.99410;
forcefield.clayff(9).charge =  -1.08090;
forcefield.clayff(9).radius =  3.553200;
forcefield.clayff(9).e_kcalmol = 0.155400;
%
% Oalh
forcefield.clayff(10).type =    {'oahe'}; % edge oxygen
forcefield.clayff(10).mass =    15.99410;
forcefield.clayff(10).charge =  -1.2375;
forcefield.clayff(10).radius =  3.553200;
forcefield.clayff(10).e_kcalmol = 0.155400;
%
% Oalhh
forcefield.clayff(11).type =    {'oahhe'}; % edge oxygen
forcefield.clayff(11).mass =    15.99410;
forcefield.clayff(11).charge =  -0.6625; %
forcefield.clayff(11).radius =  3.553200;
forcefield.clayff(11).e_kcalmol = 0.155400;
%
% Oalsi
forcefield.clayff(12).type =    {'oas'}; % 
forcefield.clayff(12).mass =    15.99410;
forcefield.clayff(12).charge =  -1.2375;
forcefield.clayff(12).radius =  3.553200;
forcefield.clayff(12).e_kcalmol = 0.155400;
%
% Osih
forcefield.clayff(13).type =    {'oshe'}; % 
forcefield.clayff(13).mass =    15.99410;
forcefield.clayff(13).charge =  -0.95;
forcefield.clayff(13).radius =  3.553200;
forcefield.clayff(13).e_kcalmol = 0.155400;
%
% Si
forcefield.clayff(14).type =    {'st'}; % st tetrahedral silicon
forcefield.clayff(14).mass =    28.08538;
forcefield.clayff(14).charge =  2.1000;
forcefield.clayff(14).radius =  3.70640;
forcefield.clayff(14).e_kcalmol = 1.8405E-06;
% Al
forcefield.clayff(15).type =    {'ao'}; % ao octahedral aluminium
forcefield.clayff(15).mass =    26.98154;
forcefield.clayff(15).charge =  1.5750;
forcefield.clayff(15).radius =  4.79430;
forcefield.clayff(15).e_kcalmol = 1.3298E-06;
%
% Alt
forcefield.clayff(16).type =    {'at'}; % at tetrahedral aluminium
forcefield.clayff(16).mass =    26.98154;
forcefield.clayff(16).charge =  1.5750;
forcefield.clayff(16).radius =  3.70640;
forcefield.clayff(16).e_kcalmol = 1.8405E-06;
% Mgo
forcefield.clayff(17).type =    {'mgo'}; % mgo octahedral magnesium
forcefield.clayff(17).mass =    24.305;
forcefield.clayff(17).charge =  1.36000;
forcefield.clayff(17).radius =  5.9090;
forcefield.clayff(17).e_kcalmol = 9.0298E-07;
%
% Mgh
forcefield.clayff(18).type =    {'mgh'}; % mgh
forcefield.clayff(18).mass =    24.305;
forcefield.clayff(18).charge =  1.0500;
forcefield.clayff(18).radius =  5.9090;
forcefield.clayff(18).e_kcalmol = 9.0298E-07;
% Cao
forcefield.clayff(19).type =    {'cao'}; % cao
forcefield.clayff(19).mass =    40.078;
forcefield.clayff(19).charge =  1.36000;
forcefield.clayff(19).radius =  6.24840;
forcefield.clayff(19).e_kcalmol = 5.02980E-06;
%
% Cah
forcefield.clayff(20).type =    {'cah'}; % cah
forcefield.clayff(20).mass =    40.078;
forcefield.clayff(20).charge =  1.05000;
forcefield.clayff(20).radius =  6.24280;
forcefield.clayff(20).e_kcalmol = 5.02980E-06;
%
% Feo
forcefield.clayff(21).type =    {'Feo'}; % feo
forcefield.clayff(21).mass =    55.845;
forcefield.clayff(21).charge =  1.5750;
forcefield.clayff(21).radius =  5.5070;
forcefield.clayff(21).e_kcalmol = 9.0298E-06;
%
% Lio
forcefield.clayff(22).type =    {'lio'}; % lio
forcefield.clayff(22).mass =    6.94;
forcefield.clayff(22).charge =  0.52500;
forcefield.clayff(22).radius =  4.72570;
forcefield.clayff(22).e_kcalmol = 9.0298E-06;
%
% Li
forcefield.clayff(23).type =    {'Li'}; %
forcefield.clayff(23).mass =    6.94;
forcefield.clayff(23).charge =  1.000;
forcefield.clayff(23).radius =  1.582;
forcefield.clayff(23).e_kcalmol = 0.3367344;
%
% Na
forcefield.clayff(24).type =    {'Na'}; %
forcefield.clayff(24).mass =    22.98977;
forcefield.clayff(24).charge =  1.000;
forcefield.clayff(24).radius =  2.424;
forcefield.clayff(24).e_kcalmol = 0.3526418;
%
% K
forcefield.clayff(25).type =    {'K'}; %
forcefield.clayff(25).mass =    39.09830;
forcefield.clayff(25).charge =  1.000;
forcefield.clayff(25).radius =  3.186;
forcefield.clayff(25).e_kcalmol = 0.4297054;
%
% Rb
forcefield.clayff(26).type =    {'Rb'}; %
forcefield.clayff(26).mass =    85.4678;
forcefield.clayff(26).charge =  1.000;
forcefield.clayff(26).radius =  3.474;
forcefield.clayff(26).e_kcalmol = 0.4451036;
%
% Cs
forcefield.clayff(27).type =    {'Cs'}; %
forcefield.clayff(27).mass =    132.90545;
forcefield.clayff(27).charge =  1.000;
forcefield.clayff(27).radius =  4.042;
forcefield.clayff(27).e_kcalmol = 0.0898565;
%
% Mg
forcefield.clayff(28).type =    {'Mg'}; %
forcefield.clayff(28).mass =    24.305;
forcefield.clayff(28).charge =  2.000;
forcefield.clayff(28).radius =  1.56880;
forcefield.clayff(28).e_kcalmol = 0.875043948;
%
% Ca
forcefield.clayff(29).type =    {'Ca'}; %
forcefield.clayff(29).mass =    40.078;
forcefield.clayff(29).charge =  2.000;
forcefield.clayff(29).radius =  2.6500;
forcefield.clayff(29).e_kcalmol = 0.449657335;
%
% Sr
forcefield.clayff(30).type =    {'Sr'}; %
forcefield.clayff(30).mass =    87.620;
forcefield.clayff(30).charge =  2.000;
forcefield.clayff(30).radius =  3.4823;
forcefield.clayff(30).e_kcalmol = 0.118225901;
%
% Ba
forcefield.clayff(31).type =    {'Ba'}; %
forcefield.clayff(31).mass =    137.330;
forcefield.clayff(31).charge =  2.000;
forcefield.clayff(31).radius =  4.24985;
forcefield.clayff(31).e_kcalmol = 0.047096425;
%
% F
forcefield.clayff(32).type =    {'F'}; %
forcefield.clayff(32).mass =    18.998403;
forcefield.clayff(32).charge =  -1.000;
forcefield.clayff(32).radius =  4.514;
forcefield.clayff(32).e_kcalmol = 0.0074005;
%
% Cl
forcefield.clayff(33).type =    {'Cl'}; %
forcefield.clayff(33).mass =    35.4530;
forcefield.clayff(33).charge =  -1.000;
forcefield.clayff(33).radius =  5.422;
forcefield.clayff(33).e_kcalmol = 0.012785;
%
% Br
forcefield.clayff(34).type =    {'Br'}; %
forcefield.clayff(34).mass =    79.9040;
forcefield.clayff(34).charge =  -1.000;
forcefield.clayff(34).radius =  5.502;
forcefield.clayff(34).e_kcalmol = 0.0269586;
%
% I
forcefield.clayff(35).type =    {'I'}; %
forcefield.clayff(35).mass =    126.90447;
forcefield.clayff(35).charge =  -1.000;
forcefield.clayff(35).radius =  5.838;
forcefield.clayff(35).e_kcalmol = 0.0427845;
%
% Oahs
forcefield.clayff(36).type =    {'oahs'}; %
forcefield.clayff(36).mass =    15.99410;
forcefield.clayff(36).charge =  -0.6625;
forcefield.clayff(36).radius =  3.553200;
forcefield.clayff(36).e_kcalmol = 0.155400;

% Du a dummy atom
forcefield.clayff(37).type =    {'Du'}; %
forcefield.clayff(37).mass =    10;
forcefield.clayff(37).charge =  0;
forcefield.clayff(37).radius =  0;
forcefield.clayff(37).e_kcalmol = 0;

% Al
forcefield.clayff(38).type =    {'aoe'}; % ao octahedral aluminium
forcefield.clayff(38).mass =    26.98154;
forcefield.clayff(38).charge =  1.8125; %1.5750+0.2375;
forcefield.clayff(38).radius =  4.79430;
forcefield.clayff(38).e_kcalmol = 1.3298E-06;

% He ( this is mostly used to account for the He atom inte clayff ff
forcefield.clayff(39).type =    {'he'}; % ho hydroxyl hydrogen
forcefield.clayff(39).mass =    1.00794;
forcefield.clayff(39).charge =  0.42500;
forcefield.clayff(39).radius =  0.00000;
forcefield.clayff(39).e_kcalmol = 0.00000;

% Ar
forcefield.clayff(40).type =    {'Ar'}; % ho hydroxyl hydrogen
forcefield.clayff(40).mass =    39.94800;
forcefield.clayff(40).charge =  0.00000;
forcefield.clayff(40).radius =  3.775333752;
forcefield.clayff(40).e_kcalmol = 0.276720296;

% Co
forcefield.clayff(41).type =    {'Co'}; %
forcefield.clayff(41).mass =    12.01073;
forcefield.clayff(41).charge =  0.6512;
forcefield.clayff(41).radius =  3.094627867;
forcefield.clayff(41).e_kcalmol = 0.05580526;

% Oc
forcefield.clayff(42).type =    {'Oc'}; %
forcefield.clayff(42).mass =    15.99941;
forcefield.clayff(42).charge =  -0.3256	;
forcefield.clayff(42).radius =  3.404427393;
forcefield.clayff(42).e_kcalmol = 0.160036571;

% Ob
forcefield.clayff(42).type =    {'o'}; % Clayffmod basal O
forcefield.clayff(42).mass =    15.99941;
forcefield.clayff(42).charge =  -1.05000;
forcefield.clayff(42).radius =  3.553200;
forcefield.clayff(42).e_kcalmol = 0.155400;

% Osi
forcefield.clayff(43).type =    {'ost'}; % Clayffmod single and lonely Si-O
forcefield.clayff(43).mass =    15.99941;
forcefield.clayff(43).charge =  -1.52500;
forcefield.clayff(43).radius =  3.553200;
forcefield.clayff(43).e_kcalmol = 0.155400;

% Oahs
forcefield.clayff(44).type =    {'oahs'}; %
forcefield.clayff(44).mass =    15.999410;
forcefield.clayff(44).charge =  -0.6625;
forcefield.clayff(44).radius =  3.553200;
forcefield.clayff(44).e_kcalmol = 0.155400;

% Oalhh
forcefield.clayff(45).type =    {'oahh'}; % edge oxygen
forcefield.clayff(45).mass =    15.999410;
forcefield.clayff(45).charge =  -0.6625; %
forcefield.clayff(45).radius =  3.553200;
forcefield.clayff(45).e_kcalmol = 0.155400;

% Osih
forcefield.clayff(46).type =    {'ohst'}; % obss bridging oxygen w. double substitution
forcefield.clayff(46).mass =    15.999410;
forcefield.clayff(46).charge =  -0.95;
forcefield.clayff(46).radius =  3.553200;
forcefield.clayff(46).e_kcalmol = 0.155400;

% Op
forcefield.clayff(46).type =    {'op'}; % op apical O, same params as ob
forcefield.clayff(46).mass =    15.99410;
forcefield.clayff(46).charge =  -1.05000;
forcefield.clayff(46).radius =  3.553200;
forcefield.clayff(46).e_kcalmol = 0.155400;

% Oal Deprotonated edge O 
forcefield.clayff(47).type =    {'oao'}; % o
forcefield.clayff(47).mass =    15.999410;
forcefield.clayff(47).charge =  -1.7650;
forcefield.clayff(47).radius =  3.553200;
forcefield.clayff(47).e_kcalmol = 0.155400;

for i=1:size(forcefield.clayff,2)
    forcefield.clayff(i).sigma=forcefield.clayff(i).radius/2^(1/6);
    forcefield.clayff(i).e_kJmol=forcefield.clayff(i).e_kcalmol*4.184;
    forcefield.clayff(i).e_eV=forcefield.clayff(i).e_kJmol*1000/6.02214149999999E23/1.60217733E-19;
end

for i = 1:length(Atom_label)
    if strncmpi(Atom_label(i),'Hw',2)
        Atom_label(i)={'Hw'};
    end
    if strncmpi(Atom_label(i),'Ow',2)
        Atom_label(i)={'Ow'};
    end
    
    ind=find(strcmpi([forcefield.clayff.type],Atom_label(i)))';
    if numel(ind)==0
        ind=find(strcmpi([forcefield.clayff.type],'Du'))';
    end
    ff.clayff(i).fftype =  [forcefield.clayff(ind).type];
    ff.clayff(i).masses =  [forcefield.clayff(ind).mass];
    ff.clayff(i).charge =  [forcefield.clayff(ind).charge];
    ff.clayff(i).epsilon =  [forcefield.clayff(ind).e_kcalmol];
    ff.clayff(i).radius =  [forcefield.clayff(ind).radius];
    ff.clayff(i).sigma =  [forcefield.clayff(ind).sigma];
end

fftype=[ff.clayff.fftype];
Masses=[ff.clayff.masses];
Charge=[ff.clayff.charge];
Epsilon=[ff.clayff.epsilon];
Radius=[ff.clayff.radius];
Sigma=[ff.clayff.sigma];

All_Bonds =  [554.1349 1.000;
    553.7600 1.000];

% All_Angles = [45.7697 109.470;
%     30.0000 109.470;
%     30.0000 109.470];

% from Teleman Jönsson, flex SPC model, eV A^2
% 24.029546 (*2) Bond
% 1.984757 % Angle

All_Angles = [45.7696 109.470];

kcalmol_to_eV = 4.184*1000/(6.0221415E+23*1.60217733E-19);

typeID = [];
e_mix = zeros(1,1);
r_mix = zeros(1,1);
k=1;

Atom_label

try
    for i=1:length(Atom_label)
        j = i;
        while j < length(Atom_label)+1 && i <= length(Epsilon) && j <= length(Epsilon)
            typeID = [typeID;[i j]];
            Epsilon(i);
            Epsilon(j);
            e_mix(k,1) = (Epsilon(i)*Epsilon(j))^.5;
            r_mix(k,1) = (Sigma(i)+Sigma(j))/2;
            j = j + 1;
            k = k + 1;
        end
    end
    assignin('caller','Epsilon',Epsilon);
    assignin('caller','Sigma',Sigma);
catch
    disp('Unable to assign all atomtypes')
    disp('We had problems with:')
    ff.clayff(i).fftype
    ff.clayff(i).masses
    ff.clayff(i).charge
    ff.clayff(i).epsilon
    ff.clayff(i).radius
    ff.clayff(i).sigma
    disp('or...')
    ff.clayff(j).fftype
    ff.clayff(j).masses
    ff.clayff(j).charge
    ff.clayff(j).epsilon
    ff.clayff(j).radius
    ff.clayff(j).sigma
end


% Send these variables to the calling script too!
assignin('caller','kcalmol_to_eV',kcalmol_to_eV);
assignin('caller','typeID',typeID);
assignin('caller','Masses',Masses);
assignin('caller','Charge',Charge);
assignin('caller','e_mix',e_mix);
assignin('caller','r_mix',r_mix);
assignin('caller','bondstr',All_Bonds);
assignin('caller','anglebend',All_Angles);
assignin('caller','Radius',Radius);
assignin('caller','ff',ff);
assignin('base','ff',ff);
assignin('caller','forcefield',forcefield);






