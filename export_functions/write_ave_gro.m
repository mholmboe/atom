%% write_ave_gro.m
% * This function reads in a number of .gro files and writes out the
% average (or median, as default) atom struct and Box_dim to a file
% ave_conf.gro, as well as the variables atom_ave and Box_dim_ave below.
%
% * Note!! If chosing the 'average' option as a second argument, make sure
% your molecule does not overlap any PBC!!! If it does, just leave the
% default option which will compute the median position instead!
%
%% Version
% 2.11
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # write_ave_gro(confstring) % Basic input arguments
% # write_ave_gro(confstring,'average')
% # write_ave_gro(confstring,'median')

function write_ave_gro(confstring,varargin)

numConf=size(dir(strcat(confstring,'*')),1);

if numConf<3
    disp('Only found a few structure files, hence the average structure may not be statistically sound...')
end

atom_all = import_atom_gro(strcat(confstring,num2str(0),'.gro')); % Also gives Box_dim
nAtoms=size(atom_all,2);

atom_all = repmat(atom_all,1,numConf);
Box_dim_all = repmat(Box_dim,numConf,1);

for i=1:numConf
    [Coord,Box] = import_ave_gro(strcat(confstring,num2str(i-1),'.gro'));
    [atom_all(1+(i-1)*nAtoms:nAtoms+(i-1)*nAtoms).x]=Coord.x;
    [atom_all(1+(i-1)*nAtoms:nAtoms+(i-1)*nAtoms).y]=Coord.y;
    [atom_all(1+(i-1)*nAtoms:nAtoms+(i-1)*nAtoms).z]=Coord.z;
    Box_dim_all(i,:)=Box;
end
Box_dim_ave=mean(Box_dim_all);

atom_ave=atom_all(1:nAtoms);
if nargin>1 && strncmpi(varargin{1},'ave',3)
    disp('Computing average structure!!!!')
    for i=1:nAtoms
        [atom_ave(i).x]=mean([atom_all(i:nAtoms:end).x]);
        [atom_ave(i).y]=mean([atom_all(i:nAtoms:end).y]);
        [atom_ave(i).z]=mean([atom_all(i:nAtoms:end).z]);
    end
else
    disp('Computing median structure!!!!')
    for i=1:nAtoms
        [atom_ave(i).x]=median([atom_all(i:nAtoms:end).x]);
        [atom_ave(i).y]=median([atom_all(i:nAtoms:end).y]);
        [atom_ave(i).z]=median([atom_all(i:nAtoms:end).z]);
    end
end

assignin('caller','atom_ave',atom_ave);
assignin('caller','Box_dim_ave',Box_dim_ave);

write_atom_gro(atom_ave,Box_dim_ave,strcat('ave_',confstring,'.gro'));

end