%% import_atom_xtc.m
% * This function imports a structure file and a xtc file
% * This function relies on Gro2Mat, see
% * Gro2mat: a package to efficiently read gromacs output in MATLAB.
% * Dien H1, Deane CM, Knapp B.
% * Journal of Computational Chemistry
% * Volume 35, Issue 20, Version of Record online: 12 JUN 2014
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * [atom,traj] = import_atom_xtc('conf.gro',traj.xtc)
% * [atom,traj] = import_atom_xtc('conf.pdb',traj.xtc)

function [atom,traj] = import_atom_xtc(filenameconf,filenamextc)
% 

if regexp(filenameconf,'.gro') > 1;
    disp('Found .gro file');
    str=strcat('gmx editconf -f',{' '},filenameconf,{' '},'-o temp.pdb &>/dev/null');
    system(strcat(char({PATH2GMX()}),char(str)));
    filenamepdb='temp.pdb';
    filenamegro=filenameconf;
    scale=10;
elseif regexp(filenameconf,'.pdb') > 1;
    disp('Found .pdb file');
    str=strcat('gmx editconf -f',{' '},filenameconf,{' '},'-o temp.gro &>/dev/null');
    system(strcat(char({PATH2GMX()}),char(str)));
    filenamepdb=filenameconf;
    filenamegro='temp.gro';
    scale=1;
end

atom=import_atom_gro(filenamegro);

trj=parseTrj(filenamepdb,filenamextc);

nAtoms=size(trj.coords,1);
nFrames=size(trj.coords,3);

% This did not work...
% xcoords={reshape(trj.coords(:,1,:)/10,nAtoms,nFrames)};
% ycoords={reshape(trj.coords(:,2,:)/10,nAtoms,nFrames)};
% zcoords={reshape(trj.coords(:,3,:)/10,nAtoms,nFrames)};
%  Box_dim={trj.trajectoryData.box};
%
% [atom(:).x]=deal(xcoords{:,:});
% [atom(:).y]=deal(xcoords{:,:});
% [atom(:).z]=deal(xcoords{:,:});
% [atom(:).Box_dim]=deal(Box_dim{:});

for atomindex=1:nAtoms
    x_coords={reshape(trj.coords(atomindex,1,:)*scale,1,nFrames)};
    [atom(atomindex).x]=deal(x_coords{:})';
    y_coords={reshape(trj.coords(atomindex,2,:)*scale,1,nFrames)};
    [atom(atomindex).y]=deal(y_coords{:})';
    z_coords={reshape(trj.coords(atomindex,3,:)*scale,1,nFrames)};
    [atom(atomindex).z]=deal(z_coords{:})';
    %     Box_dim=deal({trj.trajectoryData.box(i,:)});
    %     [traj.Box_dim]=deal(Box_dim{:});
    %     step=deal({trj.trajectoryData.step(i)});
    %     [traj.step]=deal(step{:});
    %     time=deal({trj.trajectoryData.time(i)});
    %     [traj.time]=deal(time{:});
    if mod(atomindex,100)==0;
        atomindex
    end
end

trj.trajectoryData.box=trj.trajectoryData.box*scale;
traj=trj.trajectoryData;

% Result={atom,traj};

delete('./#*'); delete('./temp*');
clear xcoords ycoords zcoords
assignin('caller','trj',trj)
assignin('caller','traj',traj)
assignin('caller','atom',atom)
assignin('caller','XYZ_labels',XYZ_labels)
assignin('caller','XYZ_data',XYZ_data)




