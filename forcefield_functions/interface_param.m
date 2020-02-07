%% interface_param.m
% * This function holds the extended interface ff parameters
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # interface_param(Atom_label)
% # interface_param(Atom_label,'tip3p')
%
function interface_param(Atom_label,varargin)

if nargin < 2
    water_model='tip3p';
else
    water_model=varargin{1};
end

water_model
%clear all;
%format short;

% Atom_label   = {'Hw', 'H','Ow','Oh','O', 'Omg', 'Oalt','Odsub','Ohmg','Oalsi','Oalhh','Oalh','Osih','Si','Al','Alt','Mgo','Mgh','Cao','Cah','Feo','Lio','Li','Na','K','Cs','Mg','Ca','Sr','Ba','Cl','Br'}';
%
% Atom_label   = {'Hw', 'Ow','Si','Al'}';

% call function with find(ismember({forcefield.interface.atom},'Hw')) to fin Hw's index

%                      1    2    3    4   5    6     7      8      9      10      11      12    13    14   15   16    17    18   19  20   21   22   23   24   25   26  27    28   29  30  31  32
% Atom_label_Heinz   = {'Hw', 'H','Ow','Oh','Ob','Op','Omg','Oalt','Ohmg','Oalsi','Oalhh','Oalh','Osih','Si','Al','Alt','Mgo','Li','Na','K','Cs','Mg','Ca','Sr','Ba','F','Cl','Cls','Br','I','He','Du'}';

% Hw
forcefield.interface(1).type =    {'HW'}; % h* water hydrogen
forcefield.interface(1).mass =    1.00794;
if strcmpi('SPC/E', water_model) > 0 || strcmpi('SPCE', water_model) > 0
    forcefield.interface(1).charge =  0.42380;
elseif strcmp('SPC', water_model) > 0
    forcefield.interface(1).charge =  0.41000;
elseif strncmpi('TIP3P', water_model,5) > 0
    forcefield.interface(1).charge =  0.41700;
end
forcefield.interface(1).radius =  0.00000; % R0 Å
forcefield.interface(1).e_kJmol = 0.00000; % kcal mol-1
%
% H
forcefield.interface(2).type =    {'H'}; % ho hydroxyl hydrogen
forcefield.interface(2).mass =    1.00794;
forcefield.interface(2).charge =  0.200;
forcefield.interface(2).radius =  1.13925;
forcefield.interface(2).e_kJmol = 0.06276;
%
% Ow
forcefield.interface(3).type =    {'OW'}; % o* water oxygen
forcefield.interface(3).mass =    15.99410;
if strcmpi('SPC/E', water_model) > 0 || strcmpi('SPCE', water_model) > 0
    forcefield.interface(3).charge =  -0.84760;
elseif strcmp('SPC', water_model) > 0
    forcefield.interface(3).charge = -0.82000;
elseif strncmpi('TIP3P', water_model,5) > 0
    forcefield.interface(3).charge =  -0.83400;
end
forcefield.interface(3).radius =  3.553200;
forcefield.interface(3).e_kJmol = 0.155400;
%
% Oh
forcefield.interface(4).type =    {'Oh'}; % Oxygen atom in octahedral aluminate sheet, hydroxl
forcefield.interface(4).mass =    15.99410;
forcefield.interface(4).charge =  -0.68333;
forcefield.interface(4).radius =  3.50000;
forcefield.interface(4).e_kJmol = 0.104600;
%
% Ob
forcefield.interface(5).type =    {'Ob'}; % Oxygen atom in silicate sheet, surface
forcefield.interface(5).mass =    15.99410;
forcefield.interface(5).charge =  -0.55000;
forcefield.interface(5).radius =  3.50000;
forcefield.interface(5).e_kJmol = 0.104600;
%
% Op
forcefield.interface(6).type =    {'Op'}; % Oxygen atom in octahedral aluminate sheet, apical
forcefield.interface(6).mass =    15.99410;
forcefield.interface(6).charge =  -0.75833;
forcefield.interface(6).radius =  3.50000;
forcefield.interface(6).e_kJmol = 0.104600;
% Omg
forcefield.interface(7).type =    {'Omg'}; % Oxygen atom in octahedral aluminate sheet, apical, if next to Mg defect
forcefield.interface(7).mass =    15.99410;
forcefield.interface(7).charge =  -0.86666;
forcefield.interface(7).radius =  3.50000;
forcefield.interface(7).e_kJmol = 0.104600;
%
% Oalt
forcefield.interface(8).type =    {'Oalt'}; % Oxygen atom in silicate sheet, surface, if next to Al defect
forcefield.interface(8).mass =    15.99410;
forcefield.interface(8).charge =  -0.78333;
forcefield.interface(8).radius =  3.50000;
forcefield.interface(8).e_kJmol = 0.104600;
%
% Ohmg
forcefield.interface(9).type =    {'Ohmg'}; % Oxygen atom in octahedral aluminate sheet, hydroxl, if next to Mg defect
forcefield.interface(9).mass =    15.99410;
forcefield.interface(9).charge =  -0.79166;
forcefield.interface(9).radius =  3.50000;
forcefield.interface(9).e_kJmol = 0.104600;
%
% Oalh
forcefield.interface(10).type =    {'Oalh'}; %
forcefield.interface(10).mass =    15.99410;
forcefield.interface(10).charge =  -0.75833; % Dummy charge
forcefield.interface(10).radius =  3.50000;
forcefield.interface(10).e_kJmol = 0.104600;
%
% Oalhh
forcefield.interface(11).type =    {'Oalhh'}; % obts bridging oxygen w. tetrahedral substitution
forcefield.interface(11).mass =    15.99410;
forcefield.interface(11).charge =  -0.4000; %??? Where does it come from? -0.07500;
forcefield.interface(11).radius =  3.50000;
forcefield.interface(11).e_kJmol = 0.104600;
%
% Oalsi
forcefield.interface(12).type =    {'Oalsi'}; % obss bridging oxygen w. double substitution
forcefield.interface(12).mass =    15.99410;
forcefield.interface(12).charge =  -0.75833; %
forcefield.interface(12).radius =  3.50000;
forcefield.interface(12).e_kJmol = 0.104600;
%
% Osih
forcefield.interface(13).type =    {'Osih'}; % obss bridging oxygen w. double substitution
forcefield.interface(13).mass =    15.99410;
forcefield.interface(13).charge =  -0.675;
forcefield.interface(13).radius =  3.50000;
forcefield.interface(13).e_kJmol = 0.104600;
%
% Si
forcefield.interface(14).type =    {'Si'}; % st tetrahedral silicon
forcefield.interface(14).mass =    28.08538;
forcefield.interface(14).charge =  1.1000;
forcefield.interface(14).radius =  4.20000;
forcefield.interface(14).e_kJmol = 0.20920;
% Al
forcefield.interface(15).type =    {'Al'}; % ao octahedral aluminium
forcefield.interface(15).mass =    26.98154;
forcefield.interface(15).charge =  1.45;%1.283333; % 1.45;% To smooth/average charge Al/Mg charge
forcefield.interface(15).radius =  4.41000;
forcefield.interface(15).e_kJmol = 0.20920;
% Alt
forcefield.interface(16).type =    {'Alt'}; % at
forcefield.interface(16).mass =    26.98154;
forcefield.interface(16).charge =  0.8;%1.45;
forcefield.interface(16).radius =  4.41000;
forcefield.interface(16).e_kJmol = 0.20920;
% Mgo
forcefield.interface(17).type =    {'Mgo'}; % mgo
forcefield.interface(17).mass =    24.305;
forcefield.interface(17).charge =  1.1000;%-0.002574/32;%-8.1250e-05;
forcefield.interface(17).radius =  4.41000;
forcefield.interface(17).e_kJmol = 0.20920;
% Li
forcefield.interface(18).type =    {'Li'}; %
forcefield.interface(18).mass =    6.94;
forcefield.interface(18).charge =  1.000;
forcefield.interface(18).radius =  1.582;
forcefield.interface(18).e_kJmol = 1.40890;
%
% Na
forcefield.interface(19).type =    {'Na'}; %
forcefield.interface(19).mass =    22.98977;
forcefield.interface(19).charge =  1.000;
forcefield.interface(19).radius =  2.424;
forcefield.interface(19).e_kJmol = 1.47545;
%
% K
forcefield.interface(20).type =    {'K'}; %
forcefield.interface(20).mass =    39.09830;
forcefield.interface(20).charge =  1.000;
forcefield.interface(20).radius =  3.186;
forcefield.interface(20).e_kJmol = 1.79789;
%
% Cs
forcefield.interface(21).type =    {'Cs'}; %
forcefield.interface(21).mass =    132.90545;
forcefield.interface(21).charge =  1.000;
forcefield.interface(21).radius =  4.042;
forcefield.interface(21).e_kJmol = 0.37596;
%
% Mg
forcefield.interface(22).type =    {'Mg'}; %
forcefield.interface(22).mass =    24.305;
forcefield.interface(22).charge =  2.000;
forcefield.interface(22).radius =  1.56880;
forcefield.interface(22).e_kJmol = 3.65874;
%
% Ca
forcefield.interface(23).type =    {'Ca'}; %
forcefield.interface(23).mass =    40.078;
forcefield.interface(23).charge =  2.000;
forcefield.interface(23).radius =  2.6500;
forcefield.interface(23).e_kJmol = 1.88011;
%
% Sr
forcefield.interface(24).type =    {'Sr'}; %
forcefield.interface(24).mass =    87.620;
forcefield.interface(24).charge =  2.000;
forcefield.interface(24).radius =  3.4823;
forcefield.interface(24).e_kJmol = 0.49433;
%
% Ba
forcefield.interface(25).type =    {'Ba'}; %
forcefield.interface(25).mass =    137.330;
forcefield.interface(25).charge =  2.000;
forcefield.interface(25).radius =  4.24985;
forcefield.interface(25).e_kJmol = 0.19692;
%
% F
forcefield.interface(26).type =    {'F'}; %
forcefield.interface(26).mass =    18.998403;
forcefield.interface(26).charge =  -1.000;
forcefield.interface(26).radius =  4.514;
forcefield.interface(26).e_kJmol = 0.0074005;
%
% Cl
forcefield.interface(27).type =    {'Cl'}; %
forcefield.interface(27).mass =    35.4530;
forcefield.interface(27).charge =  -1.000;
forcefield.interface(27).radius =  5.422;
forcefield.interface(27).e_kJmol = 0.05349;
%
% Cls
forcefield.interface(28).type =    {'Cls'}; %
forcefield.interface(28).mass =    35.4530;
forcefield.interface(28).charge =  -1.000;
forcefield.interface(28).radius =  5.422;
forcefield.interface(28).e_kJmol = 0.05349;
%
% Br
forcefield.interface(29).type =    {'Br'}; %
forcefield.interface(29).mass =    79.9040;
forcefield.interface(29).charge =  -1.000;
forcefield.interface(29).radius =  5.502;
forcefield.interface(29).e_kJmol = 0.11279;
%
% I
forcefield.interface(30).type =    {'I'}; %
forcefield.interface(30).mass =    126.90447;
forcefield.interface(30).charge =  -1.000;
forcefield.interface(30).radius =  5.838;
forcefield.interface(30).e_kJmol = 0.17901;
%
% H
forcefield.interface(31).type =    {'He'}; % Edge H
forcefield.interface(31).mass =    1.00794;
forcefield.interface(31).charge =  0.400;
forcefield.interface(31).radius =  1.13925;
forcefield.interface(31).e_kJmol = 0.06276;
%
% Du
forcefield.interface(32).type =    {'Du'}; % Edge H
forcefield.interface(32).mass =    10.0;
forcefield.interface(32).charge =  0.00;
forcefield.interface(32).radius =  0.00;
forcefield.interface(32).e_kJmol = 0.00;
% Oahs
forcefield.interface(33).type =    {'Oahs'}; %
forcefield.interface(33).mass =    15.99410;
forcefield.interface(33).charge =  -0.6625;
forcefield.interface(33).radius =  3.553200;
forcefield.interface(33).e_kcalmol = 0.155400;
% Oahs
forcefield.interface(34).type =    {'Omgsi'}; %
forcefield.interface(34).mass =    15.99410;
forcefield.interface(34).charge =  -0.9755;
forcefield.interface(34).radius =  3.553200;
forcefield.interface(34).e_kcalmol = 0.155400;
% Osi
forcefield.interface(35).type =    {'Osi'}; %
forcefield.interface(35).mass =    15.99410;
forcefield.interface(35).charge =  -1.2;
forcefield.interface(35).radius =  3.553200;
forcefield.interface(35).e_kcalmol = 0.155400;
% Al -edge
forcefield.interface(36).type =    {'Ale'}; % ao octahedral aluminium
forcefield.interface(36).mass =    26.98154;
forcefield.interface(36).charge =  1.8125; %1.5750+0.2375;
forcefield.interface(36).radius =  4.79430;
forcefield.interface(36).e_kcalmol = 1.3298E-06;
% Ar
forcefield.interface(37).type =    {'Ar'}; % Ar with SPC.. from Greeshma
forcefield.interface(37).mass =    39.94800;
forcefield.interface(37).charge =  0.00000;
forcefield.interface(37).radius =  3.775333752;
forcefield.interface(37).e_kcalmol = 0.276720296;
% Co
forcefield.interface(38).type =    {'Co'}; % C as in CO2 from Greeshma
forcefield.interface(38).mass =    12.01073;
forcefield.interface(38).charge =  0.6512;
forcefield.interface(38).radius =  3.094627867;
forcefield.interface(38).e_kcalmol = 0.05580526;
% Oc
forcefield.interface(39).type =    {'Oc'}; % O as in CO2 from Greeshma
forcefield.interface(39).mass =    15.99941;
forcefield.interface(39).charge =  -0.3256	;
forcefield.interface(39).radius =  3.404427393;
forcefield.interface(39).e_kcalmol = 0.160036571;

% Forcefield_index=find(ismember([forcefield.interface.type],Atom_label))';
% forcefield.interface=forcefield.interface(Forcefield_index);

for i=1:size(forcefield.interface,2)
    forcefield.interface(i).sigma=forcefield.interface(i).radius/2^(1/6);
    %     forcefield.interface(i).e_kJmol=forcefield.interface(i).e_kcalmol*4.184;
    forcefield.interface(i).e_eV=forcefield.interface(i).e_kJmol*1000/6.02214149999999E23/1.60217733E-19;
end

for i = 1:length(Atom_label)
    if strcmp(Atom_label(i),{'O'})
        Atom_label(i)={'Ob'};
    end
    if strncmpi(Atom_label(i),'Hw',2)
        Atom_label(i)={'Hw'};
    end
    
    ind=find(strcmpi([forcefield.interface.type],Atom_label(i)))';
    if numel(ind)==0
        ind=find(strcmpi([forcefield.interface.type],'Du'))';
    end
    
    ind=find(strcmpi([forcefield.interface.type],Atom_label(i)))';
    ff.interface(i).fftype =  [forcefield.interface(ind).type];
    ff.interface(i).masses =  [forcefield.interface(ind).mass];
    ff.interface(i).charge =  [forcefield.interface(ind).charge];
    ff.interface(i).epsilon =  [forcefield.interface(ind).e_kJmol];
    ff.interface(i).radius =  [forcefield.interface(ind).radius];
    ff.interface(i).sigma =  [forcefield.interface(ind).sigma];
end

fftype=[ff.interface.fftype];
Masses=[ff.interface.masses];
Charge=[ff.interface.charge];
Epsilon=[ff.interface.epsilon];
Radius=[ff.interface.radius];
Sigma=[ff.interface.sigma];

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

% if strcmp('SPC/E', water_model,'exact') > 0;
%     ind_Hw=find(ismember({forcefield.interface.type},'Hw')) =  0.4238;%    h*      Hw %SPCE water
%     ind_Ow=find(ismember({forcefield.interface.type},'Ow')) = -0.8476; %	o*      Ow       %SPCE water
% end


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
assignin('caller','fftype',fftype);
assignin('caller','forcefield',forcefield);

try
    for i=1:length(Atom_label)
        j = i;
        while j < length(Atom_label)+1
            %         Atom_label(i)
            %         Atom_label(j)
            %         typeID = [typeID;[i j]]
            %         e_mix
            %         Sigma
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
    ff.interface(i).fftype
    ff.interface(i).masses
    ff.interface(i).charge
    ff.interface(i).epsilon
    ff.interface(i).radius
    ff.interface(i).sigma
    disp('or...')
    ff.interface(j).fftype
    ff.interface(j).masses
    ff.interface(j).charge
    ff.interface(j).epsilon
    ff.interface(j).radius
    ff.interface(j).sigma
end

% XYZ_num=zeros(length(XYZ_labels),1);
% Atom_label=sort(unique([atom.type]));
% for i=1:length(Atom_label)
%     XYZ_num(ismember(XYZ_labels,Atom_label(i)))=i;
% end
%
% XYZ_radius=Sigma(XYZ_num(:))';
% A1=repmat(XYZ_radius,1,length(XYZ_radius));
% A2=repmat(XYZ_radius,1,length(XYZ_radius))';
% B=(A1+A2)*.3;
% D=(Dist_matrix-B)<0;
% [C1,C2]=sort(Dist_matrix(:,1));
% C3=XYZ_labels(C2(:));




