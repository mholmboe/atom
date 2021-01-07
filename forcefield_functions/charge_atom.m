%% charge_atom.m
% * This function tries to charge the atom according to clayff or interface ff
% * If more than 4 arguments, it tries to smear out the excess charge
% * atom is the atom struct
% * Box_dim is the box dimension vector
% * ffname is 'clayff' or 'interface'
% * watermodel is not always used but should be 'spc' 'spc/e' 'tip3p'
%
%% Version
% 2.09
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * atom = charge_atom(atom,Box_dim,'clayff','spc')
% * atom = charge_atom(atom,Box_dim,'interface','tip3p','more')

function atom = charge_atom(atom,Box_dim,ffname,watermodel,varargin)

if strcmpi(ffname,'clayff')
    clayff_param(sort(unique([atom.type])),watermodel);
    for i=1:length(atom)
        if strncmpi([atom(i).type],{'Hw'},2)
            ind=strncmpi({'Hw'},[forcefield.clayff.type],2);
        else
            ind=strcmpi([atom(i).type],[forcefield.clayff.type]);
        end
        try
            atom(i).charge=[forcefield.clayff(ind).charge];
        catch
            disp('Could not set charge for:')
            atom(i)
        end
    end
    if nargin > 4
        %         Atom_label=varargin{1}(:);
        %         Charge=cell2mat(varargin(2));
        % Total_charge = charge_clayff_atom(atom,Box_dim,{'Al' 'Mgo' 'Si' 'H'},[1.575 1.36 2.1 0.425])
        Atom_label=sort(unique([atom.type]));
        clayff_param(sort(Atom_label),watermodel);
        no_adjust_labels=[Atom_label(~strncmp(Atom_label,'O',1))];
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Ow',2))];
        no_adjust_ind=ismember(Atom_label,no_adjust_labels);
        no_adjust_ind
        Atom_label
        round(Charge,5)
        try
            atom=charge_clayff_atom(atom,Box_dim,Atom_label(no_adjust_ind),Charge(no_adjust_ind));
        catch
            disp('Could not set charge for:')
            atom(i)
        end
    else
        try
            atom=check_clayff_charge(atom);
        catch
            disp('Could not set charge for:')
            atom(i)
        end
    end
    
elseif strcmpi(ffname,'clayff_2004')
    clayff_2004_param(sort(unique([atom.type])),watermodel);
    for i=1:length(atom)
        if strncmpi([atom(i).type],{'Hw'},2)
            ind=strncmpi({'Hw'},[forcefield.clayff.type],2);
        else
            ind=strcmpi([atom(i).type],[forcefield.clayff.type]);
        end
        try
            atom(i).charge=[forcefield.clayff(ind).charge];
        catch
            disp('Could not set charge for:')
            atom(i)
        end
    end
    if nargin > 4
        %         Atom_label=varargin{1}(:);
        %         Charge=cell2mat(varargin(2));
        % Total_charge = charge_clayff_atom(atom,Box_dim,{'ao' 'mgo' 'st' 'ho'},[1.575 1.36 2.1 0.425])
        Atom_label=sort(unique([atom.type]));
        clayff_2004_param(sort(Atom_label),watermodel);
        no_adjust_labels=[Atom_label(~strncmp(Atom_label,'o',1))];
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Ow',2))];
        no_adjust_ind=ismember(Atom_label,no_adjust_labels);
        no_adjust_ind
        Atom_label
        round(Charge,5)
        try
            atom=charge_clayff_2004_atom(atom,Box_dim,Atom_label(no_adjust_ind),Charge(no_adjust_ind));
        catch
            disp('Could not set charge for:')
            atom(i)
        end
    else
        try
            atom=check_clayff_2004_charge(atom);
        catch
            disp('Could not set charge for:')
            atom(i)
        end
    end
    
elseif strcmpi(ffname,'interface')
    interface_param(unique([atom.type]),watermodel);
    for i=1:length(atom)
        if strncmpi([atom(i).type],{'Hw'},2)
            ind=strncmpi({'Hw'},[forcefield.interface.type],2);
        else
            ind=strcmp([atom(i).type],[forcefield.interface.type]);
        end
        atom(i).charge=[forcefield.interface(ind).charge];
    end
    if nargin > 4
        %         Atom_label=varargin{1}(:);
        %         Charge=cell2mat(varargin(2));
        % Total_charge = charge_interface_atom(atom,Box_dim,{'Al' 'Mgo' 'Si' 'H'},[1.575 1.36 2.1 0.425])
        Atom_label=sort(unique([atom.type]));
        interface_param(sort(Atom_label),watermodel);
        no_adjust_labels=[Atom_label(~strncmp(Atom_label,'O',1))];
%        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Ow',2))];
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Ob',2))];
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Op',2))];
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Oh',2))];
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Omg',3))];
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Ohmg',4))];
%         no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Osi',3))];
%         no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Osih',4))];
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Oalsi',5))];
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Oalhh',5))];
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Oalh',4))];
%        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Oalt',4))];
        no_adjust_ind=ismember(Atom_label,no_adjust_labels);
        Atom_label
        no_adjust_ind
        no_adjust_labels
        round(Charge,5)
        atom = charge_interface_atom(atom,Box_dim,Atom_label(no_adjust_ind),Charge(no_adjust_ind));
    else
        atom=check_interface_charge(atom);
    end
elseif strcmpi(ffname,'interface15')
    interface15_param(unique([atom.fftype]),watermodel);
    for i=1:length(atom)
        if strncmpi([atom(i).fftype],{'Hw'},2)
            ind=strncmpi({'Hw'},[forcefield.interface15.type],2);
        else
            ind=strcmp([atom(i).type],[forcefield.interface15.type]);
        end
        atom(i).charge=[forcefield.interface15(ind).charge];
    end
    if nargin > 4
        %         Atom_label=varargin{1}(:);
        %         Charge=cell2mat(varargin(2));
        % Total_charge = charge_interface_atom(atom,Box_dim,{'Al' 'Mgo' 'Si' 'H'},[1.575 1.36 2.1 0.425])
        Atom_label=sort(unique([atom.fftype]));
        interface15_param(sort(Atom_label),watermodel);
        no_adjust_labels=[Atom_label(~strncmp(Atom_label,'H',1))];
        no_adjust_labels=[Atom_label(~strncmp(Atom_label,'O',1))];
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'Ow',2))];
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'OY1',3))]; % Ob
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'OY4',3))]; % Op
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'OY5',3))]; % Omg
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'OY6',3))]; % Oh
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'OY9',3))]; % Ohmg
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'OSH',3))]; % Osih
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'OAS',3))]; % Oalsi
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'OAHH',4))];% Oalhh
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'OAH',3))]; % Oalh
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'OC23',4))]; % Osisi in SILICA
        no_adjust_labels=[no_adjust_labels Atom_label(strncmp(Atom_label,'OC24',4))]; % Osih in SILICA
        no_adjust_ind=ismember(Atom_label,no_adjust_labels);
        Atom_label(no_adjust_ind)
        Charge(no_adjust_ind)
        atom = charge_interface15_atom(atom,Box_dim,Atom_label(no_adjust_ind),Charge(no_adjust_ind));
    else
        atom=check_interface15_charge(atom);
    end
end

disp('Total charge')
Total_charge=sum([atom.charge]);
round(sum(Total_charge),5)
if round(Total_charge)~=sum(Total_charge) 
   disp('Run tweak_charge_atom() to get an integer charge of the struct')
end

% atom=tweak_charge_atom(atom);
assignin('caller','Total_charge',Total_charge);
