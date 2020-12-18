%% bond_valence_atom.m
% * This function tries to calculate the bond valence values according to
% * http://www.iucr.org/resources/data/datasets/bond-valence-parameters
% * compiled by I. David Brown, McMaster University, Ontario, Canada
% * idbrown@mcmaster.ca
% * Data set bvparm2016.cif: 2016 version, (posted 2016-11-03)
% *
% * atom is the atom struct
% * Box_dim is the box dimension vector
%
%% Version
% 2.08
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = bond_valence_atom(atom,Box_dim)
% # atom = bond_valence_atom(atom,Box_dim,1.25,2.25)
%
function atom = bond_valence_atom(atom,Box_dim,varargin)

if nargin>2
    rmaxshort=varargin{1}; % the maximum O-H bond in the neighbout search
    rmaxlong=varargin{2}; % the M-O distances in the neighbout search
else
    rmaxshort=1.25;
    rmaxlong=2.25;
end

valencestates=[];
valence_ion1=-100; % Dummy value
valence_ion2=-100; % Dummy value
if nargin>4
    valence_ion1=varargin{3};
    valencestates=valence_ion1;
end

load('bond_valence_values.mat');
atom=element_atom(atom);
[atom.type]=atom.element;
[atom.fftype]=atom.element;

% replicated=0;
% if min(Box_dim)<2*rmaxlong
%    orig_Box_dim=Box_dim;
%    atom=replicate_atom(atom,Box_dim,[2 2 2]); % Will give new replicated Box_dim variable
%    replicated=1;
% end

if ~isfield(atom,'neigh')
    atom=bond_atom(atom,Box_dim,rmaxlong);
elseif numel(atom(1).neigh.type)==0
    atom=bond_atom(atom,Box_dim,rmaxlong);
end

if numel(valence_ion1)==size(atom,2)
    for i=1:size(atom,2)
        if numel(atom(i).neigh.index)>0
            for j=1:size(atom(i).neigh.index,1)
                %             if (atom(i).neigh.dist(j)< 1.25 && strncmpi([atom(i).type],'H',1)) ||...
                %                 (atom(i).neigh.dist(j)< 1.25 && strncmpi([atom(i).neigh.type(j)],'H',1) )
                %
                %                 [mean_bv,std_bv,bv,bvalue]=bond_valence_data(atom(i).type,atom(i).neigh.type(j),atom(i).neigh.dist(j),Ion_1,Ion_2,R0,b,Valence_1,Valence_2,valence_ion1,valence_ion2);
                %                 atom(i).bv(j)=bv;
                %                 atom(i).mean_bv(j)=mean_bv;
                %             elseif atom(i).neigh.dist(j)> 1.25 && ~strncmpi([atom(i).neigh.type(j)],'H',1)
                [mean_bv,std_bv,bv,bvalue]=bond_valence_data(atom(i).type,atom(i).neigh.type(j),atom(i).neigh.dist(j),...
                                           Ion_1,Ion_2,R0,b,Valence_1,Valence_2,valencestates(i),valencestates(atom(i).neigh.index(j)));
                atom(i).bv(j)=bv;
                atom(i).mean_bv(j)=mean_bv;
                %             end
            end
        end
        if size(atom(i).neigh.type,1)==0
            atom(i).bv=0;
        end
        try
            atom(i).valence=sum(atom(i).bv(:));
            size(atom(i).neigh.type,1);
            atom(i).Rdiff=bvalue*log(atom(i).valence/round([atom(i).valence])); % Rdiff calc the average valence and from R0 - R i.e. the ideal bond minus the actual bond distance
        catch
            atom(i).valence=0;
            atom(i).Rdiff=0;
            atom(i).bv=0;
            i
            %         atom(i).valence=
        end
        
        if mod(i,100)==0
            i
        end
    end
else
    for i=1:size(atom,2)
        if numel(atom(i).neigh.index)>0
            for j=1:size(atom(i).neigh.index,1)
                %             if (atom(i).neigh.dist(j)< 1.25 && strncmpi([atom(i).type],'H',1)) ||...
                %                 (atom(i).neigh.dist(j)< 1.25 && strncmpi([atom(i).neigh.type(j)],'H',1) )
                %
                %                 [mean_bv,std_bv,bv,bvalue]=bond_valence_data(atom(i).type,atom(i).neigh.type(j),atom(i).neigh.dist(j),Ion_1,Ion_2,R0,b,Valence_1,Valence_2,valence_ion1,valence_ion2);
                %                 atom(i).bv(j)=bv;
                %                 atom(i).mean_bv(j)=mean_bv;
                %             elseif atom(i).neigh.dist(j)> 1.25 && ~strncmpi([atom(i).neigh.type(j)],'H',1)
                [mean_bv,std_bv,bv,bvalue]=bond_valence_data(atom(i).type,atom(i).neigh.type(j),atom(i).neigh.dist(j),...
                                           Ion_1,Ion_2,R0,b,Valence_1,Valence_2,valence_ion1,valence_ion2);
                atom(i).bv(j)=bv;
                atom(i).mean_bv(j)=mean_bv;
                %             end
            end
        end
        if size(atom(i).neigh.type,1)==0
            atom(i).bv=0;
        end
        try
            atom(i).valence=sum(atom(i).bv(:));
            size(atom(i).neigh.type,1);
            atom(i).Rdiff=bvalue*log(atom(i).valence/round([atom(i).valence])); % Rdiff calc the average valence and from R0 - R i.e. the ideal bond distance minus the actual bond distance
        catch
            atom(i).valence=0;
            atom(i).Rdiff=0;
            atom(i).bv=0;
            i
            %         atom(i).valence=
        end
        
        if mod(i,100)==0
            i
        end
    end
end

% if replicated==1
%     nAtoms=size(atom,2);
%     Box_dim=orig_Box_dim;
%     atom=atom(1:nAtoms/(2*2*2));
% end


Atom_labels=unique([atom.type]);
for i=1:length(Atom_labels)
    label_ind=find(strcmpi([atom.type],Atom_labels(i)));
    assignin('caller',strcat(char(Atom_labels(i)),'_bv'),[atom(label_ind).bv]');
    assignin('caller',strcat(char(Atom_labels(i)),'_valence'),[atom(label_ind).valence]');
    assignin('caller',strcat(char(Atom_labels(i)),'_Rdiff'),[atom(label_ind).Rdiff]');
end


A=round([atom.valence]);
B=[atom.valence];
C(1:2:2*numel(A))=A;
C(2:2:2*numel(A))=B;
% C;

% disp('    Mean   |  Median  |  std ')
% [mean([atom.valence]-round([atom.valence])) median([atom.valence]-round([atom.valence])) std([atom.valence]-round([atom.valence]))]

