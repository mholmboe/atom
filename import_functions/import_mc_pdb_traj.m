%% import_mc_pdb_traj.m
% * This function imports several .pdb traj files into a full traj, and can
% even handle trajs with varying number of particles
% * varargin can be used to set max frames to import, or set a value for
% nevery frame to import, i.e. a stride value
%
%% Version
% 2.0
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_mc_pdb_traj('molecule.pdb')
% # atom = import_mc_pdb_traj('molecule.pdb',1000)
% # atom = import_mc_pdb_traj('molecule.pdb',10000,10)
%
function atom = import_mc_pdb_traj(filename,varargin)

% maxframes=100;
% stride=1;
% filenames={...
%     'Movie_MMTx2_1.1.1_298.000000_101325.000000_frameworks.pdb';...
%     'Movie_MMTx2_1.1.1_298.000000_101325.000000_component_Ca_0.pdb';...
%     'Movie_MMTx2_1.1.1_298.000000_101325.000000_component_spce_1.pdb';...
%     'Movie_MMTx2_1.1.1_298.000000_101325.000000_component_glycol_2.pdb';...
%     };

for i=1:size(filename,1)
    i
    if nargin>1
        maxframes=varargin{1};
        if nargin>2
            stride=varargin{2};
            try
                temp_atom = import_pdb_traj(char(filename(i)),maxframes,stride);
            catch
                temp_atom = import_pdb_traj(filename,maxframes,stride);
            end
        else
            try
                temp_atom = import_pdb_traj(char(filename(i)),maxframes);
            catch
                temp_atom = import_pdb_traj(filename,maxframes);
            end
        end
    else
        try
            temp_atom = import_pdb_traj(char(filename(i)));
        catch
            temp_atom = import_pdb_traj(filename);
        end
    end
    
    if i==1;framework=temp_atom;framework_traj=traj;end % traj comes from import_pdb_traj() above
    if i==2;component0=temp_atom;component0_traj=traj;end
    if i==3;component1=temp_atom;component1_traj=traj;end
    if i==4;component2=temp_atom;component2_traj=traj;end
    if i==5;component3=temp_atom;component3_traj=traj;end
    if i==6;component4=temp_atom;component4_traj=traj;end
    if i>6;disp('You have to many components for this little script...'); pause;end
end

if size(filename,2)==1;atom=framework;end
if size(filename,2)==2;atom=update_atom({framework component0});end
if size(filename,2)==3;atom=update_atom({framework component0 component1});end
if size(filename,2)==4;atom=update_atom({framework component0 component1 component2});end
if size(filename,2)==5;atom=update_atom({framework component0 component1 component2 component3});end
if size(filename,2)==6;atom=update_atom({framework component0 component1 component2 component3 component4});end
if size(filename,2)==7;atom=update_atom({framework component0 component1 component2 component3 component4 component5});end
if size(filename,2)==8;atom=update_atom({framework component0 component1 component2 component3 component4 component5 component6});end

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
XYZ_labels=[atom.type]';

if size(filename,2)==1;traj=framework_traj;end
if size(filename,2)==2;traj=[framework_traj component0_traj];end
if size(filename,2)==3;traj=[framework_traj component0_traj component1_traj];end
if size(filename,2)==4;traj=[framework_traj component0_traj component1_traj component2_traj];end
if size(filename,2)==5;traj=[framework_traj component0_traj component1_traj component2_traj component3_traj];end
if size(filename,2)==6;traj=[framework_traj component0_traj component1_traj component2_traj component3_traj component4_traj];end
if size(filename,2)==7;traj=[framework_traj component0_traj component1_traj component2_traj component3_traj component4_traj component5_traj];end
if size(filename,2)==8;traj=[framework_traj component0_traj component1_traj component2_traj component3_traj component4_traj component5_traj component6_traj];end

if nargin>3
    filename_out=varargin{3};
    write_pdb_traj(atom,traj,Box_dim,filename_out)
end

assignin('caller','atom',atom);
assignin('caller','traj',traj);
assignin('caller','nAtoms',size(atom,2));
assignin('caller','Box_dim',Box_dim);
assignin('caller','XYZ_labels',XYZ_labels);
assignin('caller','XYZ_data',XYZ_data);
