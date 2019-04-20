%% atom2mk_angndx.m
% * This function can help you print one custom gromacs angle.ndx group,
% * based on the atomtypes names. You just have to pass on the three (or
% * four) atomtypes that makes up the angle, like {'Mgo' 'Ohmg' 'H'}, or 
% * {'Al' 'Mgo' 'Oh' 'H'}, which will find the angle between all 
% * Mgo/Al - Oh - H. Note that you can use strncmpi to get a wider search.
%
%% Version
% 2.03
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% * atom2mk_angndx(filename,groupname,atomtypes,molid)

function atom2mk_angndx(filename,groupname,atomtypes)

% clear all; format compact;
% You could set the group name here manually if you want
% filename='evap_0.gro'
% groupname='MMT';
% atomtypes={'Mgo' 'Ohmg' 'H'};

atom = import_atom(filename);
atom = bond_angle_atom(atom,Box_dim,1.25,2.25);

if iscell(atomtypes)
    atomtypes=atomtypes;
else % if atomtypes
    ind=atomtypes;
end

if iscell(atomtypes)
    if numel(atomtypes)==3
        ind=[];
        str=sort(atomtypes);
        for i =1:size(Angle_index,1)
            temp_str=sort([atom(Angle_index(i,1:3)).type]);
            % Note that you can use strncmpi to get a wider search term
            if sum(strcmpi(temp_str,str))==3
                ind=[ind Angle_index(i,1:3)];
            end
        end
    elseif numel(atomtypes)==4
        ind=[];
        str1=sort(atomtypes([1 3 4]));
        str2=sort(atomtypes([2 3 4]));
        for i =1:size(Angle_index,1)
            temp_str=sort([atom(Angle_index(i,1:3)).type]);
            % Note that you can use strncmpi to get a wider search term
            if sum(strcmpi(temp_str,str1))==3
                ind=[ind Angle_index(i,1:3)];
            end
        end
        for i =1:size(Angle_index,1)
            temp_str=sort([atom(Angle_index(i,1:3)).type]);
            % Note that you can use strncmpi to get a wider search term
            if sum(strcmpi(temp_str,str2))==3
                ind=[ind Angle_index(i,1:3)];
            end
        end
    end
end


% Format the ind vector to have 15 entries per row
ext_ind=zeros(1,15*ceil(length(ind)/15));
ext_ind(1:length(ind))=ind;
ext_ind=reshape(ext_ind,15,[])';

% Print the index file
fid = fopen(strcat(groupname,'.ndx'), 'wt');
fprintf(fid, '%s %s %s\r\n','[',groupname,']');
for i = 1:size(ext_ind,1)
    row=ext_ind(i,:);
    row(row==0)=[];
    if max(row)<1000
        fprintf(fid, '%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i\r\n', row);
    elseif max(row)<10000
        fprintf(fid, '%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i\r\n', row);
    elseif max(row)<100000
        fprintf(fid, '%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i\r\n', row);
    end
end
fprintf(fid, '\r\n');
fclose(fid);
