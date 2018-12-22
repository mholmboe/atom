%% import_xtcv2.m
% * This function imports an xtc trajectory file.
% * Does it work?
% * This function imports a structure file and a xtc file
% * This function relies on Gro2Mat, see
% * Gro2mat: a package to efficiently read gromacs output in MATLAB.
% * Dien H1, Deane CM, Knapp B.
% * Journal of Computational Chemistry
% * Volume 35, Issue 20, Version of Record online: 12 JUN 2014
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * [atom,traj] = import_xtcv2('conf.gro',traj.xtc)
% * [atom,traj] = import_xtcv2('conf.pdb',traj.xtc)

function atom = import_xtcv2(filenameconf,filenamextc)
%%

if regexp(filenameconf,'.gro') > 1
    disp('Found .gro file, generating a .pdb file');
    atom=import_atom_gro(filenameconf);
    filenamepdb='temp.pdb';
    write_atom_pdb(atom,Box_dim,filenamepdb);
elseif regexp(filenameconf,'.pdb') > 1
    disp('Found .pdb file');
    filenamepdb=filenameconf;
    atom=import_atom_pdb(filenamepdb);
end

disp('Found .xtc file, will use the excellent Gro2Mat scripts');
disp('See Journal of Computational Chemistry,Volume 35, Issue 20')
disp(' ')
disp(' ')
disp('Note that in maxTol, startFrame, endFrame values can be set,')
disp('look into Gro2Mat parseTrj function for help')
disp('Note that lines 55-56 in parseTrj can give unneccesary problems')

tic
trj=parseTrj(filenamepdb,filenamextc,0.001);
toc
nAtoms=size(trj.coords,1);
nFrames=size(trj.coords,3);

%% This did not work...
xcoords=(reshape(trj.coords(:,1,:)*10,nAtoms,nFrames)).';
ycoords=(reshape(trj.coords(:,2,:)*10,nAtoms,nFrames)).';
zcoords=(reshape(trj.coords(:,3,:)*10,nAtoms,nFrames)).';
traj=zeros(size(xcoords,1),size(xcoords,2)*3);
traj(:,1:3:end)=xcoords;
traj(:,2:3:end)=ycoords;
traj(:,3:3:end)=zcoords;

trj.trajectoryData.box=trj.trajectoryData.box*10;
trajdata=trj.trajectoryData;
Box_dim=double(trajdata.box);
time=trajdata.time;
Step=trajdata.step;

delete('./#*'); delete('./temp*');
% clear xcoords ycoords zcoords
assignin('caller','trj',trj)
assignin('caller','traj',trajdata)
assignin('caller','atom',atom)
assignin('caller','Box_dim',Box_dim)
assignin('caller','traj',traj)
assignin('caller','XYZ_labels',XYZ_labels)
assignin('caller','XYZ_data',XYZ_data)
