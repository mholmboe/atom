%% ionize_atom.m
% * This function adds ions within a certain region defined by <limits_variable.html limits>
% * Can also add particles on a plane by setting Lx|Ly|Lz to 0 or something small
% * Compared to create_atom, this function can also add particles near a
% 'surface' or in the 'bulk', when an in_atom struct (representing a solid
% phase) is passed argument.
% * If slow, check out insert_atom or solvate_atom or grid_atom...
%
%% Function arguments
% * {type} is particle/atomtype
% * {resname} is resname
% * [limits] is a 1x3 or 1x6 volume variable
% * The number nmax is the max number of particles
% * Optional scale number, varargin{1}) is a how-many-diameters-between-the-particles-thingy
%
%% Dependencies
% * radius_ion
% * add2atom
% * replicate_atom
% * translate_atom
% * merge_atom
% * slice_atom
% * update_atom
% * distance_matrix_atom
%
%% Similar
% fuse_atom
% protonate_atom
% create_atom
% insert_atom
%
%% Version
% 2.082
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = ionize_atom('Na','Na',[10 20 30],10) % Basic input arguments
% # atom = ionize_atom('Na','Na',[10 20 30],10,2) % Nearest distance will be 2 * ionic radii
% # atom = ionize_atom('Na','Na',[10 20 30],10,2,in_atom) % Random placement
% # atom = ionize_atom('Na','Na',[10 20 30],10,2,in_atom,'surface') % Preferred placement at the 'surface' or 'bulk'
% # atom = ionize_atom('Na','Na',[10 20 30],10,2,in_atom,'surface'|'bulk','x'|'y'|'z'|20) % Preferred placement at the 'surface' or 'bulk' within the x|y|z or [value] range
%
function atom = ionize_atom(type,resname,limits,nmax,varargin)

sorted_distratio=[];

if iscell(type)==0;type={type}; end
if iscell(resname)==0;resname={resname};end

radii = abs(radius_ion(type));
if nargin > 4
    distance_factor=varargin{1};
    if numel(distance_factor)>1
        rmin=distance_factor(2);
        distance_factor=distance_factor(1);
    else
        rmin=radii*distance_factor;
    end
else
    distance_factor=2;
    rmin=radii*distance_factor;
end

if distance_factor<2
    disp('distance_factor should be >=2 !!!')
    pause(2)
end
Box_dim_temp=distance_factor*2*[radii radii radii];
atom = add2atom(type,[0 0 0],resname,[]);

if numel(limits)==1
    Lx=limits(1);
    Ly=limits(1);
    Lz=limits(1);
    limits(4)=limits(1);
    limits(5)=limits(1);
    limits(6)=limits(1);
    limits(1:3)=0;
elseif numel(limits)==3
    Lx=limits(1);
    Ly=limits(2);
    Lz=limits(3);
    limits(4)=limits(1);
    limits(5)=limits(2);
    limits(6)=limits(3);
    limits(1:3)=0;
elseif numel(limits)==6
    Lx=limits(4)-limits(1);
    Ly=limits(5)-limits(2);
    Lz=limits(6)-limits(3);
elseif numel(limits)==9
    Lx=limits(1);
    Ly=limits(2);
    Lz=limits(3);
    limits(4)=limits(1);
    limits(5)=limits(2);
    limits(6)=limits(3);
    limits(1:3)=0;
end

nx=ceil(Lx/Box_dim_temp(1));
ny=ceil(Ly/Box_dim_temp(2));
nz=ceil(Lz/Box_dim_temp(3));

if nx==0;nx=1;end
if ny==0;ny=1;end
if nz==0;nz=1;end

atom=replicate_atom(atom,Box_dim_temp,[nx ny nz]); % nx/ny/nz==0 is set to 1 in replicate_atom

molid=num2cell(1:size(atom,2));
[atom.molid]=deal(molid{:});

% Move the particles around a little bit
for i=1:size(atom,2)
    if nx>0 && (limits(4)-limits(1))>5;atom(i).x=atom(i).x-rmin*(rand(1)-0.5);end
    if ny>0 && (limits(5)-limits(2))>5;atom(i).y=atom(i).y-rmin*(rand(1)-0.5);end
    if nz>0 && (limits(6)-limits(3))>5;atom(i).z=atom(i).z-rmin*(rand(1)-0.5);end
end


if (limits(1)+limits(2)+limits(3)) ~= 0
    disp('Translating the water box');
    atom=translate_atom(atom,[limits(1) limits(2) limits(3)],'all');
end

disp('nAtom before merge')
size(atom,2)

if nargin>5 && size(varargin{2},2) > 0
    in_atom=varargin{2};
    if size(atom,2) > 10000 || size(in_atom,2) > 20000
        natom_block=size(atom,2)/(nx*ny*nz);
        atom_count=1;atom_merged=[];count=1;
        while atom_count< size(atom,2)
            atom_block= atom(atom_count:atom_count+natom_block-1);
            atom_block = merge_atom(in_atom,limits(4:6),atom_block,'type',type,rmin);
            atom_merged = [atom_merged atom_block];
            atom_count=atom_count+natom_block;
            disp('box number...')
            count=count+1
        end
        atom=atom_merged;
    else
        atom = merge_atom(in_atom,limits(4:6),atom,'type',type,rmin);
    end
else
    atom = slice_atom(atom,limits,0);
end

% assignin('base','atom3',atom);
atom=update_atom(atom);

% Set nmax
if iscellstr({nmax}) == 1
    nmax=size(atom,2);
end

% Randomize order of the particles
nAtoms=size(atom,2);
ind_rand=randperm(nAtoms);
ind_sel=ismember(ind_rand,1:nAtoms);
atom_ind=ind_rand(ind_sel);
atom(atom_ind)=atom;

% Check that no added particles are too close
distmatrix=dist_matrix_atom(atom,Box_dim);
distmatrix(distmatrix==0)=1000000; % Dummy distance in the distance matrix
i=1;ind=[];
while i<size(atom,2)+1
    [minvalue,ind]=min(distmatrix(i,:));
    if minvalue<distance_factor*radii
        ind=[ind i];
        %         atom(ind)=[];         v2.06
        %         distmatrix(ind,:)=[]; v2.06
        %         distmatrix(:,ind)=[]; v2.06
    else
        i=i+1;
    end
    %     if mod(i,1000)==1
    %         if i > 1
    %             [i-1 size(atom,2)]
    %         end
    %     end
    i=i+1;
end
if numel(ind)>0
    atom(ind)=[];
end

if nargin>6
    if nargin>7
        if strncmpi(varargin{4},'x',1)
            D=limits(4)-limits(1);
        elseif strncmpi(varargin{4},'y',1)
            D=limits(5)-limits(2);
        elseif strncmpi(varargin{4},'z',1)
            D=limits(6)-limits(3);
        else
            D=varargin{4}; % Set to a number in Ångström
        end
    else
        D=limits(6)-limits(3); % i.e. along Z
    end
    
    % We do this to limit the number of atoms in in_atom-->surface_atom before calculating the distance_matrix
    indxlo=find([in_atom.x]<limits(1));
    indxhi=find([in_atom.x]>limits(4));
    indylo=find([in_atom.y]<limits(2));
    indyhi=find([in_atom.y]>limits(5));
    indzlo=find([in_atom.z]<limits(3));
    indzhi=find([in_atom.z]>limits(6));
    ind = unique([indxlo indxhi indylo indyhi indzlo indzhi]);
    surface_atom=in_atom;
    surface_atom(ind)=[];
    
    distmatrix=dist_matrix_atom(atom,surface_atom,limits);
    i=1;distratio=zeros(1,size(atom,2));distmatrix(distmatrix==0)=1000000; % Dummy distance in the distance matrix
    while i<size(atom,2)+1 %>nmax+1
        if strncmpi(varargin{3},'surface',1)
            dist=min(distmatrix(i,:));
            distratio(i)=exp(2*dist/D)/rand(1);
        else % if strncmpi(varargin{3},'bulk')
            dist=min(distmatrix(i,:));
            distratio(i)=exp(-dist/(D*2))/rand(1);
        end
        i=i+1;
    end
    [sorted_distratio,distratio_order]=sort(distratio);
    atom=atom(distratio_order(1:nmax));
else
    %     i=1;
    distmatrix=dist_matrix_atom(atom,Box_dim);
    distmatrix(distmatrix==0)=1000000; % Dummy distance in the distance matrix
    while size(atom,2)>nmax+1
        [row,col]=find(distmatrix==min(min(distmatrix)));
        ind_rm=max([row(1) col(1)]);
        %         if ind_rm>i
        %             i=i+1;
        %         end
        atom(row(1))=[];
        distmatrix(row(1),:)=[];
        distmatrix(:,col(1))=[];
        
        if mod(size(atom,2),100)==1
            if size(atom,2) > 1
                [size(atom,2) nmax]
            end
        end
        
        
    end
end

% Do we need this? Delete particles if not using the <maxion> option
if iscellstr({nmax}) == 0
    if nmax > size(atom,2)
        disp('Ooops, you asked for too many particles...')
        disp('Max number of particles allowed without changing scale is:')
        size(atom,2)
        atom=atom(1:nmax);
    else
        try
            atom=atom(1:nmax);
        catch
            disp('ionize_atom didd not manage to add enough particles!!');
        end
    end
end

disp('nIon after merge')
size(atom,2)

atom=update_atom(atom);

assignin('caller','limits',limits);
assignin('caller','distratio',sorted_distratio);



