%% oplsaa_go_param.m - This function holds the extended oplsaa_aa ff for graphite oxide
function oplsaa_go_param(Atom_label,water_model)
%% 

%clear all;
%format short;

% Atom_label   = {'Hw', 'H','Ow','Oh','O', 'Omg', 'Oalt','Odsub','Ohmg','Oalsi','Oalhh','Oalh','Osih','Si','Al','Alt','Mgo','Mgh','Cao','Cah','Feo','Lio','Li','Na','K','Cs','Mg','Ca','Sr','Ba','Cl','Br'}';
%
% Atom_label   = {'Hw', 'Ow','Si','Al'}';

%% call function with find(ismember({forcefield.oplsaa_go.atom},'Hw')) to fin Hw's index

%%                      1    2    3    4   5    6     7      8      9      10      11      12    13    14   15   16    17    18   19  20   21   22   23   24   25   26  27    28   29  30  31  32
Atom_label_Heinz   = {'Hw', 'H','Ow','Oh','Ob','Op','Omg','Oalt','Ohmg','Oalsi','Oalhh','Oalh','Osih','Si','Al','Alt','Mgo','Li','Na','K','Cs','Mg','Ca','Sr','Ba','F','Cl','Cls','Br','I','He','Du'}';

%% Hw
forcefield.oplsaa_go(1).type =    {'H'}; % h* OH hydrogen
forcefield.oplsaa_go(1).fftype =  {'opls_155'};
forcefield.oplsaa_go(1).mass =    1.008;
forcefield.oplsaa_go(1).charge =  0.41800;
forcefield.oplsaa_go(1).radius =  0.00000; % R0 Å
forcefield.oplsaa_go(1).e_kJmol = 0.00000; % kcal mol-1
%%
%% H
forcefield.oplsaa_go(2).type =    {'Oh'}; % ho hydroxyl hydrogen
forcefield.oplsaa_go(2).fftype =  {'opls_154'};
forcefield.oplsaa_go(2).mass =    15.99410;
forcefield.oplsaa_go(2).charge = -0.6830;
forcefield.oplsaa_go(2).radius =  3.12000*2^(1/6); % 3.12000e-01  7.11280e-01
forcefield.oplsaa_go(2).e_kJmol = 7.11280e-01/4.184;
%%
%% Ow
forcefield.oplsaa_go(3).type =    {'Oe'}; % o* water oxygen
forcefield.oplsaa_go(3).fftype =  {'opls_180'};
forcefield.oplsaa_go(3).mass =    15.99410;
forcefield.oplsaa_go(3).charge = -0.4000;
forcefield.oplsaa_go(3).radius =  2.90000*2^(1/6); % 2.90000e-01  5.85760e-01
forcefield.oplsaa_go(3).e_kJmol = 5.85760e-01/4.184;
%%
%% Oh
forcefield.oplsaa_go(4).type =    {'Cen'}; % Oxygen atom in octahedral aluminate sheet, hydroxl
forcefield.oplsaa_go(4).fftype =  {'opls_141'};
forcefield.oplsaa_go(4).mass =    12.0110;
forcefield.oplsaa_go(4).charge =  0.0000;
forcefield.oplsaa_go(4).radius =  3.55000*2^(1/6); % 3.55000e-01  3.17984e-01
forcefield.oplsaa_go(4).e_kJmol = 3.17984e-01/4.184;
%%
%% Ob
forcefield.oplsaa_go(5).type =    {'Ce'}; % Oxygen atom in silicate sheet, surface
forcefield.oplsaa_go(5).fftype =  {'opls_183'};
forcefield.oplsaa_go(5).mass =    12.0110;
forcefield.oplsaa_go(5).charge =  0.0000;
forcefield.oplsaa_go(5).radius =  3.50000*2^(1/6); % 3.50000e-01  2.76144e-01
forcefield.oplsaa_go(5).e_kJmol = 2.76144e-01/4.184;
%%
%% Op
forcefield.oplsaa_go(6).type =    {'Coh'}; % Oxygen atom in octahedral aluminate sheet, apical
forcefield.oplsaa_go(6).fftype =  {'opls_159'};
forcefield.oplsaa_go(6).mass =    12.0110;
forcefield.oplsaa_go(6).charge =  0;
forcefield.oplsaa_go(6).radius =  3.50000*2^(1/6); % 3.50000e-01  2.76144e-01
forcefield.oplsaa_go(6).e_kJmol = 2.76144e-01/4.184;


% Forcefield_index=find(ismember([forcefield.oplsaa_go.type],Atom_label))';
% forcefield.oplsaa_go=forcefield.oplsaa_go(Forcefield_index);

for i=1:size(forcefield.oplsaa_go,2);
    forcefield.oplsaa_go(i).sigma=forcefield.oplsaa_go(i).radius/2^(1/6);
    %     forcefield.oplsaa_go(i).e_kJmol=forcefield.oplsaa_go(i).e_kcalmol*4.184;
    forcefield.oplsaa_go(i).e_eV=forcefield.oplsaa_go(i).e_kJmol*1000/6.02214149999999E23/1.60217733E-19;
end

for i = 1:length(Atom_label)
    if strcmp(Atom_label(i),{'O'})
        Atom_label(i)={'Ob'};
    end
    if strncmpi(Atom_label(i),'Hw',2)
        Atom_label(i)={'Hw'};
    end
    
    
    ind=find(strcmpi([forcefield.oplsaa_go.type],Atom_label(i)))';
    if numel(ind)==0;
        ind=find(strcmpi([forcefield.oplsaa_go.type],'Du'))';
    end
    
    ind=find(strcmpi([forcefield.oplsaa_go.type],Atom_label(i)))';
    ff.oplsaa_go(i).fftype =  [forcefield.oplsaa_go(ind).type];
    ff.oplsaa_go(i).masses =  [forcefield.oplsaa_go(ind).mass];
    ff.oplsaa_go(i).charge =  [forcefield.oplsaa_go(ind).charge];
    ff.oplsaa_go(i).epsilon =  [forcefield.oplsaa_go(ind).e_kJmol];
    ff.oplsaa_go(i).radius =  [forcefield.oplsaa_go(ind).radius];
    ff.oplsaa_go(i).sigma =  [forcefield.oplsaa_go(ind).sigma];
end

fftype=[ff.oplsaa_go.fftype];
Masses=[ff.oplsaa_go.masses];
Charge=[ff.oplsaa_go.charge];
Epsilon=[ff.oplsaa_go.epsilon];
Radius=[ff.oplsaa_go.radius];
Sigma=[ff.oplsaa_go.sigma];

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
%     ind_Hw=find(ismember({forcefield.oplsaa_go.type},'Hw')) =  0.4238;%    h*      Hw %SPCE water
%     ind_Ow=find(ismember({forcefield.oplsaa_go.type},'Ow')) = -0.8476; %	o*      Ow       %SPCE water
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
    ff.oplsaa_go(i).fftype
    ff.oplsaa_go(i).masses
    ff.oplsaa_go(i).charge
    ff.oplsaa_go(i).epsilon
    ff.oplsaa_go(i).radius
    ff.oplsaa_go(i).sigma
    disp('or...')
    ff.oplsaa_go(j).fftype
    ff.oplsaa_go(j).masses
    ff.oplsaa_go(j).charge
    ff.oplsaa_go(j).epsilon
    ff.oplsaa_go(j).radius
    ff.oplsaa_go(j).sigma
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




