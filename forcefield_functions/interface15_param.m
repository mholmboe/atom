%% interface15_param.m
% * This function holds the extended interface ff parameters
% * Call function with find(ismember({forcefield.interface15.atom},'Hw')) to find Hw's index
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples

%% Examples
% * interface15_param({'AY1'},'tip3p')


function interface15_param(Atom_label,varargin)
%

if nargin == 1
    water_model='TIP3P'
else
    water_model=varargin{1}
end

if nargin > 2
    model_database=varargin{2}
else
    model_database='CLAY_MINERALS'
end

% Original from charmm27_interface15_v1_5.prm with comparison to clayff (MHolmobes atomnames...)
% !CLAY MINERALS (KAOLINITE, MICA, MONTMORILLONITES, PYROPHYLLITE)
% !
% Heinz, H.; Koerner, H.; Anderson, K. L.; Vaia, R. A. and Farmer, B. L., Force Field for Mica-Type
% Silicates and Dynamics of Octadecylammonium Chains Grafted to Montmorillonite, Chemistry of Materials,
% 2005, 17, 5658-5669; and extensions to kaolinite.
%                                                                                                           ; Mofied Clayff (MHolmboe)
% K+     Potassium ion in mica, mmt, and other clays 		(+1.0)											; K 			1
% NA+    Sodium ion in silica, clay, and organic salts	    (+1.0)											; Na 			1
% SY1    Silicon atom in tetrahedral silicate sheet 		(+1.1)											; Si 			2.1
% SY2    Silicon atom in tetrahedral silicate sheet 		(+1.1)											; Si 			2.1
% AYT1   Aluminum defect of SY1 							(+0.8)											; Alt 			1.575
% AYT2   Aluminum defect of SY2 							(+0.8)											; Alt 			1.575
% AY1    Aluminum atom in octahedral aluminate sheet		(+1.45) or magnesium defect	(+1.1)		 		; Al or Mgo		1.575 or 1.36
% AY2    Aluminum atom in octahedral aluminate sheet		(+1.45) or magnesium defect	(+1.1)		 		; Al or Mgo		1.575 or 1.36
% OY1    Oxygen atom in silicate sheet, surface 			(-0.55(000), -0.783(33) if next to Al defect)	; O or Omg		-1.05 or -1.1688
% OY2    Oxygen atom in silicate sheet, surface 			(-0.55(000), -0.783(33) if next to Al defect)	; O or Omg		-1.05 or -1.1688
% OY3    Oxygen atom in silicate sheet, surface 			(-0.55(000), -0.783(33) if next to Al defect)	; O or Omg		-1.05 or -1.1688
% OY4    Oxygen atom in octahedral aluminate sheet, apical 	(-0.758(33), -0.866(66) if next to Mg defect)	; O or Omg		-1.05 or -1.1808
% OY5    Oxygen atom in octahedral aluminate sheet, apical 	(-0.758(33), -0.866(66) if next to Mg defect)	; O or Omg		-1.05 or -1.1808
% OY6    Oxygen atom in octahedral aluminate sheet, hydroxl	(-0.683(33), -0.791(66) if next to Mg defect)	; Oh or Ohmg	-0.95 or -1.0808
% OY7    Oxygen atom in octahedral aluminate sheet, apical 	(-0.758(33), -0.866(66) if next to Mg defect)	; O or Omg		-1.05 or -1.1808
% OY8    Oxygen atom in octahedral aluminate sheet, apical 	(-0.758(33), -0.866(66) if next to Mg defect)	; O or Omg		-1.05 or -1.1808
% OY9    Oxygen atom in octahedral aluminate sheet, hydroxl (-0.683(33), -0.791(66) if next to Mg defect)	; Oh or Ohmg	-0.95 or -1.0808
% HOY    Hydrogen atom in (Al,Mg,Si)OH groups in clay 		(+0.2), silica (+0.4), cement, silica (excl. ettr); H		    0.425
% HOK    Hydrogen atom in AlOH surface group in kaolinite 	(+0.2)											; H				0.425
%
% SILICA
% SC4    Silicon atom in tetrahedral silicate sheet 		(+1.1)											; Si in silica
% OC23   Oxygen atom in silicate sheet, surface 			(-0.55)                                         ; Osi in silica
% OC24   Oxygen atom in silicate sheet, surface 			(-0.675)	                                    ; Osih in silica
% HOY    Hydrogen atom in silica          		            (+0.4)                                          ; H


% New additions by M Holmboe
% MY1   Octahedral Magnesium defect, replacing Aluminum atom in octahedral aluminate sheet (+1.1)		 	; Al or Mgo		1.575 or 1.36
% Edge O atoms starts here... by M Holmboe
% OAH
% OAHH
% OSH
% OAS
% OASH

% from interface15_param.m  1    2    3    4   5    6     7      8      9      10      11      12    13    14   15   16    17    18   19  20   21   22   23   24   25   26  27    28   29   30  31   32
%    Atom_label_Heinz =  {'Hw', 'H','Ow','Oh','Ob','Op','Omg','Oalt','Ohmg','Oalsi','Oalhh','Oalh','Osih','Si','Al','Alt','Mgo','Li','Na','K','Cs','Mg','Ca','Sr','Ba','F','Cl','Cls','Br','I','He','Du'}';

% Hw
forcefield.interface15(1).type =    {'HW'}; % h* water hydrogen
forcefield.interface15(1).mass =    1.00794;
if strcmpi('SPC/E', water_model) > 0 || strcmpi('SPCE', water_model) > 0
    forcefield.interface15(1).charge =  0.42380;
elseif strcmp('SPC', water_model) > 0
    forcefield.interface15(1).charge =  0.41000;
elseif strncmpi('TIP3P', water_model,5) > 0
    forcefield.interface15(1).charge =  0.41700;
end
forcefield.interface15(1).radius =  0.00000; % R0 Å
forcefield.interface15(1).e_kJmol = 0.00000; % kcal mol-1
%
% H
forcefield.interface15(2).type =    {'HOY'}; % ho hydroxyl hydrogen
if strcmp(model_database,'CLAY_MINERALS')
    forcefield.interface15(2).mass =    1.00794;
    forcefield.interface15(2).charge =  0.200;
    forcefield.interface15(2).radius =  1.13925;
    forcefield.interface15(2).e_kJmol = 0.06276;
elseif strcmpi(model_database,'SILICA')
    forcefield.interface15(2).mass =    1.00794;
    forcefield.interface15(2).charge =  0.400;
    forcefield.interface15(2).radius =  0.00000;
    forcefield.interface15(2).e_kJmol = 0.00000;
end
%
% Ow
forcefield.interface15(3).type =    {'OW'}; % o* water oxygen
forcefield.interface15(3).mass =    15.99410;
if strcmpi('SPC/E', water_model) > 0 || strcmpi('SPCE', water_model) > 0
    forcefield.interface15(3).charge =  -0.84760;
elseif strcmp('SPC', water_model) > 0
    forcefield.interface15(3).charge = -0.82000;
elseif strncmpi('TIP3P', water_model,5) > 0
    forcefield.interface15(3).charge =  -0.83400;
end
forcefield.interface15(3).radius =  3.553200;
forcefield.interface15(3).e_kJmol = 0.155400;
%
% Oh
forcefield.interface15(4).type =    {'OY6'}; % Oxygen atom in octahedral aluminate sheet, hydroxl
forcefield.interface15(4).mass =    15.99410;
forcefield.interface15(4).charge =  -0.68333;
forcefield.interface15(4).radius =  3.67500;
forcefield.interface15(4).e_kJmol = 0.104600;
%
% Ob
forcefield.interface15(5).type =    {'OY1'}; % Oxygen atom in silicate sheet, surface
forcefield.interface15(5).mass =    15.99410;
forcefield.interface15(5).charge =  -0.55000;
forcefield.interface15(5).radius =  3.67500;
forcefield.interface15(5).e_kJmol = 0.104600;
%
% Op
forcefield.interface15(6).type =    {'OY4'}; % Oxygen atom in octahedral aluminate sheet, apical
forcefield.interface15(6).mass =    15.99410;
forcefield.interface15(6).charge =  -0.75833;
forcefield.interface15(6).radius =  3.67500;
forcefield.interface15(6).e_kJmol = 0.104600;
% Omg
forcefield.interface15(7).type =    {'OY5'}; % Oxygen atom in octahedral aluminate sheet, apical, if next to Mg defect
forcefield.interface15(7).mass =    15.99410;
forcefield.interface15(7).charge =  -0.86666;
forcefield.interface15(7).radius =  3.67500;
forcefield.interface15(7).e_kJmol = 0.104600;
%
% Oalt
forcefield.interface15(8).type =    {'OY2'}; % Oxygen atom in silicate sheet, surface, if next to Al defect
forcefield.interface15(8).mass =    15.99410;
forcefield.interface15(8).charge =  -0.78333;
forcefield.interface15(8).radius =  3.67500;
forcefield.interface15(8).e_kJmol = 0.104600;
%
% % Op
% forcefield.interface15(6).type =    {'OY3'}; % Oxygen atom in octahedral aluminate sheet, apical
% forcefield.interface15(6).mass =    15.99410;
% forcefield.interface15(6).charge =  -0.75833;
% forcefield.interface15(6).radius =  3.67500;
% forcefield.interface15(6).e_kJmol = 0.104600;
% %

% % Oh
% forcefield.interface15(4).type =    {'OY7'}; % Oxygen atom in octahedral aluminate sheet, hydroxl
% forcefield.interface15(4).mass =    15.99410;
% forcefield.interface15(4).charge =  -0.68333;
% forcefield.interface15(4).radius =  3.67500;
% forcefield.interface15(4).e_kJmol = 0.104600;
% %
% % Oh
% forcefield.interface15(4).type =    {'OY8'}; % Oxygen atom in octahedral aluminate sheet, hydroxl
% forcefield.interface15(4).mass =    15.99410;
% forcefield.interface15(4).charge =  -0.68333;
% forcefield.interface15(4).radius =  3.67500;
% forcefield.interface15(4).e_kJmol = 0.104600;
% %
% Ohmg
forcefield.interface15(9).type =    {'OY9'}; % Oxygen atom in octahedral aluminate sheet, hydroxl, if next to Mg defect
forcefield.interface15(9).mass =    15.99410;
forcefield.interface15(9).charge =  -0.79166;
forcefield.interface15(9).radius =  3.67500;
forcefield.interface15(9).e_kJmol = 0.104600;
%
% Oalh
forcefield.interface15(10).type =    {'OAH'}; %
forcefield.interface15(10).mass =    15.99410;
forcefield.interface15(10).charge =  -0.75833; % Dummy charge
forcefield.interface15(10).radius =  3.67500;
forcefield.interface15(10).e_kJmol = 0.104600;
%
% Oalhh
forcefield.interface15(11).type =    {'OAHH'}; % obts bridging oxygen w. tetrahedral substitution
forcefield.interface15(11).mass =    15.99410;
forcefield.interface15(11).charge =  -0.4000; %??? Where does it come from? -0.07500;
forcefield.interface15(11).radius =  3.67500;
forcefield.interface15(11).e_kJmol = 0.104600;
%
% Oalsi
forcefield.interface15(12).type =    {'OAS'}; % obss bridging oxygen w. double substitution
forcefield.interface15(12).mass =    15.99410;
forcefield.interface15(12).charge =  -0.75833; %
forcefield.interface15(12).radius =  3.67500;
forcefield.interface15(12).e_kJmol = 0.104600;
%
% Osih
forcefield.interface15(13).type =    {'OSH'}; % obss bridging oxygen w. double substitution
forcefield.interface15(13).mass =    15.99410;
forcefield.interface15(13).charge =  -0.675;
forcefield.interface15(13).radius =  3.67500;
forcefield.interface15(13).e_kJmol = 0.104600;
%
% Si
forcefield.interface15(14).type =    {'SY1'}; % st tetrahedral silicon
forcefield.interface15(14).mass =    28.08538;
forcefield.interface15(14).charge =  1.1000;
forcefield.interface15(14).radius =  4.20000;
forcefield.interface15(14).e_kJmol = 0.20920;
%
%
% % Si
% forcefield.interface15(14).type =    {'SY2'}; % st tetrahedral silicon
% forcefield.interface15(14).mass =    28.08538;
% forcefield.interface15(14).charge =  1.1000;
% forcefield.interface15(14).radius =  4.20000;
% forcefield.interface15(14).e_kJmol = 0.20920;
%
% Al
forcefield.interface15(15).type =    {'AY1'}; % ao octahedral aluminium
forcefield.interface15(15).mass =    26.98154;
forcefield.interface15(15).charge =  1.45;% 1.283334; % To smooth/average charge Al/Mg charge
forcefield.interface15(15).radius =  4.41000;
forcefield.interface15(15).e_kJmol = 0.20920;
%
% % Al
% forcefield.interface15(15).type =    {'AY2'}; % ao octahedral aluminium
% forcefield.interface15(15).mass =    26.98154;
% forcefield.interface15(15).charge =  1.45;% 1.283334; % To smooth/average charge Al/Mg charge
% forcefield.interface15(15).radius =  4.41000;
% forcefield.interface15(15).e_kJmol = 0.20920;
%
% Alt
forcefield.interface15(16).type =    {'AYT1'}; % at tetrahedral aluminium
forcefield.interface15(16).mass =    26.98154;
forcefield.interface15(16).charge =  1.45;
forcefield.interface15(16).radius =  4.41000;
forcefield.interface15(16).e_kJmol = 0.20920;
% Mgo
forcefield.interface15(17).type =    {'MY1'}; % mgo
forcefield.interface15(17).mass =    24.305;
forcefield.interface15(17).charge =  1.1000;%-0.002574/32;%-8.1250e-05;
forcefield.interface15(17).radius =  4.41000;
forcefield.interface15(17).e_kJmol = 0.20920;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if strcmpi(model_database,'SILICA')
%     disp('Model database is SILICA')
% Si in SILICA
forcefield.interface15(18).type =    {'SC4'}; % st tetrahedral silicon
forcefield.interface15(18).mass =    28.08538;
forcefield.interface15(18).charge =  1.1000; % 0.725 when dehydrated
forcefield.interface15(18).radius =  4.20000;
forcefield.interface15(18).e_kJmol = 0.20920;

% Si in SILICA
forcefield.interface15(19).type =    {'SC5'}; % st tetrahedral silicon
forcefield.interface15(19).mass =    28.08538;
forcefield.interface15(19).charge =  0.725; % when dehydrated and bonded to OC25
forcefield.interface15(19).radius =  4.20000;
forcefield.interface15(19).e_kJmol = 0.20920;

% Osi in SILICA
forcefield.interface15(20).type =    {'OC23'}; %
forcefield.interface15(20).mass =    15.99410;
forcefield.interface15(20).charge =  -0.5500;
forcefield.interface15(20).radius =  3.67500;
forcefield.interface15(20).e_kJmol = 0.104600;

% Osih in SILICA
forcefield.interface15(21).type =    {'OC24'}; %
forcefield.interface15(21).mass =    15.99410;
forcefield.interface15(21).charge =  -0.67500; % -0.9 when dehydrated
forcefield.interface15(21).radius =  3.67500;
forcefield.interface15(21).e_kJmol = 0.104600;

% Osi- in SILICA
forcefield.interface15(22).type =    {'OC25'}; %
forcefield.interface15(22).mass =    15.99410;
forcefield.interface15(22).charge =  -0.9000; % -0.9 dehydrated as SiO-
forcefield.interface15(22).radius =  3.67500;
forcefield.interface15(22).e_kJmol = 0.104600;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ion pair-potentials from Joung-Cheatham 2008 and Babu-Lim, 20??
% Li
forcefield.interface15(23).type =    {'Li'}; %
forcefield.interface15(23).mass =    6.94;
forcefield.interface15(23).charge =  1.000;
forcefield.interface15(23).radius =  1.582;
forcefield.interface15(23).e_kJmol = 1.40890;
%
% Na
forcefield.interface15(24).type =    {'Na'}; %
forcefield.interface15(24).mass =    22.98977;
forcefield.interface15(24).charge =  1.000;
forcefield.interface15(24).radius =  2.424;
forcefield.interface15(24).e_kJmol = 1.47545;
%
% K
forcefield.interface15(25).type =    {'K'}; %
forcefield.interface15(25).mass =    39.09830;
forcefield.interface15(25).charge =  1.000;
forcefield.interface15(25).radius =  3.186;
forcefield.interface15(25).e_kJmol = 1.79789;
%
% Cs
forcefield.interface15(26).type =    {'Cs'}; %
forcefield.interface15(26).mass =    132.90545;
forcefield.interface15(26).charge =  1.000;
forcefield.interface15(26).radius =  4.042;
forcefield.interface15(26).e_kJmol = 0.37596;
%
% Mg
forcefield.interface15(27).type =    {'Mg'}; %
forcefield.interface15(27).mass =    24.305;
forcefield.interface15(27).charge =  2.000;
forcefield.interface15(27).radius =  1.56880;
forcefield.interface15(27).e_kJmol = 3.65874;
%
% Ca
forcefield.interface15(28).type =    {'Ca'}; %
forcefield.interface15(28).mass =    40.078;
forcefield.interface15(28).charge =  2.000;
forcefield.interface15(28).radius =  2.6500;
forcefield.interface15(28).e_kJmol = 1.88011;
%
% Sr
forcefield.interface15(29).type =    {'Sr'}; %
forcefield.interface15(29).mass =    87.620;
forcefield.interface15(29).charge =  2.000;
forcefield.interface15(29).radius =  3.4823;
forcefield.interface15(29).e_kJmol = 0.49433;
%
% Ba
forcefield.interface15(30).type =    {'Ba'}; %
forcefield.interface15(30).mass =    137.330;
forcefield.interface15(30).charge =  2.000;
forcefield.interface15(30).radius =  4.24985;
forcefield.interface15(30).e_kJmol = 0.19692;
%
% F
forcefield.interface15(31).type =    {'F'}; %
forcefield.interface15(31).mass =    18.998403;
forcefield.interface15(31).charge =  -1.000;
forcefield.interface15(31).radius =  4.514;
forcefield.interface15(31).e_kJmol = 0.0074005;
%
% Cl
forcefield.interface15(32).type =    {'Cl'}; %
forcefield.interface15(32).mass =    35.4530;
forcefield.interface15(32).charge =  -1.000;
forcefield.interface15(32).radius =  5.422;
forcefield.interface15(32).e_kJmol = 0.05349;
%
% Cls
forcefield.interface15(33).type =    {'Cls'}; %
forcefield.interface15(33).mass =    35.4530;
forcefield.interface15(33).charge =  -1.000;
forcefield.interface15(33).radius =  5.422;
forcefield.interface15(33).e_kJmol = 0.05349;
%
% Br
forcefield.interface15(34).type =    {'Br'}; %
forcefield.interface15(34).mass =    79.9040;
forcefield.interface15(34).charge =  -1.000;
forcefield.interface15(34).radius =  5.502;
forcefield.interface15(34).e_kJmol = 0.11279;
%
% I
forcefield.interface15(35).type =    {'I'}; %
forcefield.interface15(35).mass =    126.90447;
forcefield.interface15(35).charge =  -1.000;
forcefield.interface15(35).radius =  5.838;
forcefield.interface15(35).e_kJmol = 0.17901;
%
% H
forcefield.interface15(36).type =    {'HE'}; % Edge H
forcefield.interface15(36).mass =    1.00794;
forcefield.interface15(36).charge =  0.400;
forcefield.interface15(36).radius =  1.13925;
forcefield.interface15(36).e_kJmol = 0.06276;
%
% Du
forcefield.interface15(37).type =    {'Du'}; % Edge H
forcefield.interface15(37).mass =    10.0;
forcefield.interface15(37).charge =  0.00;
forcefield.interface15(37).radius =  0.00;
forcefield.interface15(37).e_kJmol = 0.00;
% Oahs
forcefield.interface15(38).type =    {'OAHS'}; %
forcefield.interface15(38).mass =    15.99410;
forcefield.interface15(38).charge =  -0.6625;
forcefield.interface15(38).radius =  3.553200;
forcefield.interface15(38).e_kcalmol = 0.155400;
% Oahs
forcefield.interface15(39).type =    {'OMS'}; %
forcefield.interface15(39).mass =    15.99410;
forcefield.interface15(39).charge =  -0.9755;
forcefield.interface15(39).radius =  3.553200;
forcefield.interface15(39).e_kcalmol = 0.155400;
% Al -edge
forcefield.interface15(40).type =    {'AYE'}; % ao octahedral aluminium
forcefield.interface15(40).mass =    26.98154;
forcefield.interface15(40).charge =  1.8125; %1.5750+0.2375;
forcefield.interface15(40).radius =  4.79430;
forcefield.interface15(40).e_kcalmol = 1.3298E-06;
% Ar
forcefield.interface15(41).type =    {'Ar'}; % Ar with SPC..
forcefield.interface15(41).mass =    39.94800;
forcefield.interface15(41).charge =  0.00000;
forcefield.interface15(41).radius =  3.775333752;
forcefield.interface15(41).e_kcalmol = 0.276720296;
% Co
forcefield.interface15(42).type =    {'Co'}; % C as in CO2
forcefield.interface15(42).mass =    12.01073;
forcefield.interface15(42).charge =  0.6512;
forcefield.interface15(42).radius =  3.094627867;
forcefield.interface15(42).e_kcalmol = 0.05580526;
% Oc
forcefield.interface15(43).type =    {'Oc'}; % O as in CO2
forcefield.interface15(43).mass =    15.99941;
forcefield.interface15(43).charge =  -0.3256	;
forcefield.interface15(43).radius =  3.404427393;
forcefield.interface15(43).e_kcalmol = 0.160036571;


% Forcefield_index=find(ismember([forcefield.interface15.type],Atom_label))';
% forcefield.interface15=forcefield.interface15(Forcefield_index);

for i=1:size(forcefield.interface15,2)
    forcefield.interface15(i).sigma=forcefield.interface15(i).radius/2^(1/6);
    %     forcefield.interface15(i).e_kJmol=forcefield.interface15(i).e_kcalmol*4.184;
    forcefield.interface15(i).e_eV=forcefield.interface15(i).e_kJmol*1000/6.02214149999999E23/1.60217733E-19;
end

for i = 1:length(Atom_label)
    if strcmp(Atom_label(i),{'O'})
        Atom_label(i)={'OY1'};
    end
    if strncmpi(Atom_label(i),'Hw',2)
        Atom_label(i)={'Hw'};
    end
    
    ind=find(strcmpi([forcefield.interface15.type],Atom_label(i)))';
    if numel(ind)==0
        ind=find(strcmpi([forcefield.interface15.type],'Du'))';
    end
    
    ind=find(strcmpi([forcefield.interface15.type],Atom_label(i)))';
    ff.interface15(i).fftype =  [forcefield.interface15(ind).type];
    ff.interface15(i).masses =  [forcefield.interface15(ind).mass];
    ff.interface15(i).charge =  [forcefield.interface15(ind).charge];
    ff.interface15(i).epsilon =  [forcefield.interface15(ind).e_kJmol];
    ff.interface15(i).radius =  [forcefield.interface15(ind).radius];
    ff.interface15(i).sigma =  [forcefield.interface15(ind).sigma];
end

fftype  = [ff.interface15.fftype];
Masses  = [ff.interface15.masses];
Charge  = [ff.interface15.charge];
Epsilon = [ff.interface15.epsilon];
Radius  = [ff.interface15.radius];
Sigma   = [ff.interface15.sigma];

% All_Bonds =  [554.1349 1.000;
%               553.7600 1.000];

% All_Angles = [45.7697 109.470;
%     30.0000 109.470;
%     30.0000 109.470];

% from Teleman Jönsson, flex SPC model, eV A^2
% 24.029546 (*2) Bond
% 1.984757 % Angle

% All_Angles = [45.7696 109.470];

kcalmol_to_eV = 4.184*1000/(6.0221415E+23*1.60217733E-19);

typeID = [];
e_mix = zeros(1,1);
r_mix = zeros(1,1);
k=1;

% if strcmp('SPC/E', water_model,'exact') > 0;
%     ind_Hw=find(ismember({forcefield.interface15.type},'Hw')) =  0.4238;%    h*      Hw %SPCE water
%     ind_Ow=find(ismember({forcefield.interface15.type},'Ow')) = -0.8476; %	o*      Ow       %SPCE water
% end

% Send these variables to the calling script too!
assignin('caller','kcalmol_to_eV',kcalmol_to_eV);
assignin('caller','typeID',typeID);
assignin('caller','Masses',Masses);
assignin('caller','Charge',Charge);
assignin('caller','e_mix',e_mix);
assignin('caller','r_mix',r_mix);
% assignin('caller','bondstr',All_Bonds);
% assignin('caller','anglebend',All_Angles);
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
    i
    j
    Atom_label(i)
    Atom_label(j)
    ff.interface15(i).fftype
    ff.interface15(i).masses
    ff.interface15(i).charge
    ff.interface15(i).epsilon
    ff.interface15(i).radius
    ff.interface15(i).sigma
    disp('or...')
    ff.interface15(j).fftype
    ff.interface15(j).masses
    ff.interface15(j).charge
    ff.interface15(j).epsilon
    ff.interface15(j).radius
    ff.interface15(j).sigma
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




