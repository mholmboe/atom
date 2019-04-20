%% cell_list_dist_matrix_atom.m
% * This function calculates the distance matrix from the atom struct,
% using a cell list algorithm adapted from the Matlab MDtoolbox by
% Yasuhiro Matsunaga, link to github
% * http://github.com/ymatsunaga/mdtoolbox/
%
% Original reference to the cell list algorithm:
% Heinz, T.N. & Hunenberger, P.H.
% "A fast pairlist-construction algorithm for molecular simulations
% under periodic boundary conditions."
% J Comput Chem 25,  Sep;25(12):1474-1486 (2004).
%
%% Version
% 2.03
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # dist_matrix = cell_list_dist_matrix_atom(atom,Box_dim)
% # dist_matrix = cell_list_dist_matrix_atom(atom,Box_dim,1.25) % Setting the max distance for bonds with H's
% # dist_matrix = cell_list_dist_matrix_atom(atom,Box_dim,1.25,2.25) % Setting the max distance for bonds, otherwise default is 12 Å
% # dist_matrix = cell_list_dist_matrix_atom(atom,Box_dim,1.25,2.25,2.25) % Setting the max distance for angles, default is not to calculate any angles
%
function dist_matrix = cell_list_dist_matrix_atom(atom,Box_dim,varargin) % ,atom2,Box_dim); % or % ,Box_dim);

nAtoms=size(atom,2);
dist_matrix=zeros(nAtoms);
X_dist=dist_matrix;
Y_dist=dist_matrix;
Z_dist=dist_matrix;

if nargin > 2
    rHmax=varargin{1};
else
    rHmax=1.25;
end

if nargin >3
    calc_bond_list=1;
    rmax=varargin{2};
else
    calc_bond_list=0;
    rmax=15;
end

if nargin > 4
    calc_angle_list=1;
    rmax_angle=varargin{3};
else
    calc_angle_list=0;
    %   rmax_angle=2.25;
end

rmax2=rmax.^2;

if numel(Box_dim)==1
    Box_dim(1)=Box_dim(1);
    Box_dim(2)=Box_dim(1);
    Box_dim(3)=Box_dim(1);
    xy=0;xz=0;yz=0;
elseif numel(Box_dim)==3
    xy=0;xz=0;yz=0;
else
    xy=Box_dim(6);xz=Box_dim(8);yz=Box_dim(9);
end

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
if numel(Box_dim)==3
    XYZ_data(:,1) = XYZ_data(:,1) - Box_dim(1)*floor(XYZ_data(:,1)./Box_dim(1));
    XYZ_data(:,2) = XYZ_data(:,2) - Box_dim(2)*floor(XYZ_data(:,2)./Box_dim(2));
    XYZ_data(:,3) = XYZ_data(:,3) - Box_dim(3)*floor(XYZ_data(:,3)./Box_dim(3));
elseif numel(Box_dim)==9 && numel(nonzeros(Box_dim(1,4:end))) > 0 %&& sum(abs(Box_dim([6 8 9]))) > 0
    % Untested
    XYZ_data(:,1) = XYZ_data(:,1) - xy*floor(XYZ_data(:,2)./Box_dim(2));
    XYZ_data(:,1) = XYZ_data(:,1) - xz*floor(XYZ_data(:,3)./Box_dim(3));
    XYZ_data(:,2) = XYZ_data(:,2) - yz*floor(XYZ_data(:,3)./Box_dim(3));
    XYZ_data(:,1) = XYZ_data(:,1) - Box_dim(1)*floor(XYZ_data(:,1)./Box_dim(1));
    XYZ_data(:,2) = XYZ_data(:,2) - Box_dim(2)*floor(XYZ_data(:,2)./Box_dim(2));
    XYZ_data(:,3) = XYZ_data(:,3) - Box_dim(3)*floor(XYZ_data(:,3)./Box_dim(3));
end

if any(Box_dim(1:3) < (2*rmax))
    disp('Decrease rmax!!! The Box_dim are too small...')
    disp('Assuming you have a small system, reverting to ')
    disp('the basic bond_angle_atom() function... ')
    pause(1)
    dist_matrix = dist_matrix_atom(atom,Box_dim);
    if nargin>5
        atom=bond_angle_atom(atom,Box_dim,1.25,min(Box_dim(1:3))/2,'more');
        assignin('caller',char(varargin{4}),atom);
    end
    assignin('caller','dist_matrix',dist_matrix)
    return
end

rcell = max([8.0 8.0 8.0], rmax);
ncell = floor(Box_dim(1:3)./rcell); % [Nx Ny Nz] discretized grid cells
rcell = Box_dim(1:3)./ncell; % dimensions of the discretized grid cells

nx = floor(XYZ_data(:,1)./rcell(1))'; % x positions index
ny = floor(XYZ_data(:,2)./rcell(2))'; % y positions index
nz = floor(XYZ_data(:,3)./rcell(3))'; % z positions index

% Number of grid cells
M = prod(ncell);

% Grid cell indexes
m = ncell(1)*ncell(2)*nz + ncell(1)*ny + nx + 1;

mindex = cell(M, 1);
mindex3 = cell(M, 1);
for i = 1:M
    mindex{i} = find(m == i);
    ind3 = zeros(1, numel(mindex{i})*3);
    ind3(1:3:end) = 3.*(mindex{i}-1)+1;
    ind3(2:3:end) = 3.*(mindex{i}-1)+2;
    ind3(3:3:end) = 3.*(mindex{i}-1)+3;
    mindex3{i} = ind3;%to3(mindex{i});
end

% calculate cell mask
dm = 1:(M - 1);
dmx = mod(dm, ncell(1));
dmy = floor(mod(dm, ncell(1)*ncell(2)) ./ ncell(1));
dmz = floor(dm ./ (ncell(1)*ncell(2)));

dnx = abs(minimum_image(dmx, ncell(1)));

dny = zeros(size(dmy));
logical_index = (dmx == 0) | (dmy == (ncell(2)-1) & dmz == (ncell(3)-1));
dny(logical_index) = abs(minimum_image(dmy(logical_index), ncell(2)));
dny(~logical_index) = min(abs(minimum_image(dmy(~logical_index), ncell(2))), abs(minimum_image(dmy(~logical_index) + 1, ncell(2))));

dnz = zeros(size(dmz));
logical_index = (dmz == (ncell(3)-1)) | (dmx == 0 & dmy == 0);
dnz(logical_index) = abs(minimum_image(dmz(logical_index), ncell(3)));
dnz(~logical_index) = min(abs(minimum_image(dmz(~logical_index), ncell(3))), abs(minimum_image(dmz(~logical_index) + 1, ncell(3))));

mask = zeros(size(dm));
mask = (max(dnx, 1) - 1).^2 * rcell(1).^2 + (max(dny, 1) - 1).^2 * rcell(2).^2 + (max(dnz, 1) - 1).^2 * rcell(3).^2 <= rmax2;
%mask = true(1, M - 1); % bug check of mask
mask_index = find(mask);

% calculate the distances of the atoms in masked cells
pair = zeros(nAtoms*1000, 2);
dist = zeros(nAtoms*1000, 1);
icount = 1;
for m1 = 1:M
    m1index = mindex{m1};
    if numel(m1index) == 0
        continue
    elseif numel(m1index) > 1
        [lpair,ldist,num,rx,ry,rz] = calcpair_atom_pbc(XYZ_data(m1index,:), rmax, Box_dim); % New
        if num > 0
            X_dist(m1index,m1index)=-rx;
            Y_dist(m1index,m1index)=-ry;
            Z_dist(m1index,m1index)=-rz;
            pair(icount:(icount+num-1), :) = [m1index(lpair(:,1)')' m1index(lpair(:,2)')'];
            dist(icount:(icount+num-1)) = ldist;
            icount = icount + num;
        end
    end
    
    m2 = m1 + mask_index;
    m2 = m2(m2 <= M);
    if isempty(m2)
        continue
    end
    m2index = [mindex{m2}];
    [lpair,ldist,num,rx,ry,rz] = calcpair2_atom_pbc(XYZ_data(m1index,:), XYZ_data(m2index,:), rmax, Box_dim); % New
    if num > 0
        X_dist(m1index,m2index)=-rx;
        Y_dist(m1index,m2index)=-ry;
        Z_dist(m1index,m2index)=-rz;
        X_dist(m2index,m1index)=rx';
        Y_dist(m2index,m1index)=ry';
        Z_dist(m2index,m1index)=rz';
        pair(icount:(icount+num-1),:) = [m1index(lpair(:,1)')' m2index(lpair(:,2)')'];
        dist(icount:(icount+num-1)) = ldist;
        icount = icount + num;
    end
end

pair(icount:end,:) = [];
dist(icount:end) = [];

% sort pair index (upper triangular form)
for i = 1:size(pair,1)
    if pair(i,1) > pair(i,2)
        tmp = pair(i,1);
        pair(i,1) = pair(i,2);
        pair(i,2) = tmp;
    end
    dist_matrix(pair(i,1),pair(i,2))=dist(i);
    dist_matrix(pair(i,2),pair(i,1))=dist(i); % Do we need this?
end

%         assignin('caller','pair',pair)
%         assignin('caller','dist',dist)

if calc_bond_list==1
    Bond_index = [pair dist];
    % Special treatment for the short H-distances
    H_ind=find(strncmpi([atom.type],'H',1)); % Find all H indexes
    Hcol1=find(ismember(Bond_index(:,1),H_ind)); % Find all H indexes in first column
    Hcol2=find(ismember(Bond_index(:,2),H_ind)); % Find all H indexes in second column
    H_ind=sort([Hcol1; Hcol2]);
    Hrmax_ind=find(Bond_index(:,3)>rHmax); % Find all too long H interactions
    H2_ind=intersect(Hcol1,Hcol2); % Find all H-H interactions
    H_rm_ind=unique(sort([Hrmax_ind; H2_ind])); % Clean up all H indexes in the Bond_index
    
    Bond_index(intersect(H_ind,H_rm_ind),:)=[]; % Remove all H bonds that are either too long or H-H
    Bond_index(Bond_index(:,3)>rmax,:)=[]; % Remove all other bonds that are too long
    [Y,I] = sort(Bond_index(:,2));
    Bond_index = Bond_index(I,:);
    Bond_index = unique(Bond_index,'rows','stable');
    Bond_index(~any(Bond_index,2),:) = [];
    [Y,I] = sort(Bond_index(:,1));
    Bond_index = Bond_index(I,:);
    
    i=1;
    while i<size(Bond_index,1)
        if atom(Bond_index(i,1)).molid==atom(Bond_index(i,2)).molid
            i=i+1;
        else
            Bond_index(i,:)=[];
        end
    end
    assignin('caller','Bond_index',Bond_index)
end

if calc_angle_list==1
    a=1;
    Angle_dist_index=Bond_index(Bond_index(:,3)<rmax_angle,:);
    Angle_dist_index=[Angle_dist_index;[Angle_dist_index(:,2) Angle_dist_index(:,1) Angle_dist_index(:,3)]];
    Angle_index=zeros(6*size(XYZ_data,1),10);
    mid_angle_ind=unique(Angle_dist_index(:,1));
    for u=1:length(mid_angle_ind)
        angle_ind=find(Angle_dist_index(:,1)==mid_angle_ind(u));
        neigh_ind=Angle_dist_index(angle_ind,2);
        [pair, dist, num,rx,ry,rz] = calcpair2_atom_pbc(XYZ_data(mid_angle_ind(u),:),XYZ_data(neigh_ind,:),rmax_angle,Box_dim);
        neigh_vec=[rx' ry' rz'];
        for v=1:size(neigh_ind,1)
            for w=1:size(neigh_ind,1)
                angle=rad2deg(atan2(norm(cross(neigh_vec(v,:),neigh_vec(w,:))),dot(neigh_vec(v,:),neigh_vec(w,:))));
                if angle > 60 && angle <= 150
                    if v < w
                        Angle_index(a,1)= neigh_ind(v,1);
                        Angle_index(a,2)= mid_angle_ind(u);
                        Angle_index(a,3)= neigh_ind(w,1);
                        Angle_index(a,4)= angle;
                        Angle_index(a,5:7)= neigh_vec(v,:);
                        Angle_index(a,8:10)= neigh_vec(w,:);
                        a=a+1;
                    else
                        Angle_index(a,1)= neigh_ind(w,1);
                        Angle_index(a,2)= mid_angle_ind(u);
                        Angle_index(a,3)= neigh_ind(v,1);
                        Angle_index(a,4)= angle;
                        Angle_index(a,5:7)= neigh_vec(w,:);
                        Angle_index(a,8:10)= neigh_vec(v,:);
                        a=a+1;
                    end
                end
            end
        end
    end
    [Y,I]=sort(Angle_index(:,2));
    Angle_index=Angle_index(I,:);
    Angle_index = unique(Angle_index,'rows','stable');
    Angle_index(~any(Angle_index,2),:) = [];
    assignin('caller','Angle_index',Angle_index)
end

if nargin > 5
    radius=repmat(rmax,nAtoms,1);
    % radius([find(strncmpi([atom.type],'Ow',2)) find(strncmpi([atom.type],'H',1))])=max_short_dist;
    % This is a test
    % radius_ion(find(strncmpi([atom.type],'H',1)))=max_short_dist;
    % dist_matrix = dist_matrix_atom(atom,Box_dim);
    for i=1:size(atom,2)
        Neigh_ind=intersect(find(dist_matrix(:,i)>0),find(dist_matrix(:,i)<radius(i)));%radius_ion([atom.type])));
        Neigh_dist=dist_matrix(Neigh_ind,i);
        atom(i).neigh.type = [atom(Neigh_ind).type]';
        atom(i).neigh.index = Neigh_ind;
        atom(i).neigh.dist = Neigh_dist;
        atom(i).neigh.coords = [[atom(Neigh_ind).x]' [atom(Neigh_ind).y]' [atom(Neigh_ind).z]']; % PBC NOT TAKEN INTO ACCOUNT!!!
        atom(i).bond.type = [];
        atom(i).bond.index = [];
        atom(i).bond.dist = [];
        atom(i).angle.type = [];
        atom(i).angle.index = [];
        atom(i).angle.angle = [];
        atom(i).angle.vec1 = [];
        atom(i).angle.vec2 = [];
        
        if ismember(i,Bond_index(:,1:2))
            [A,B]=find(Bond_index(:,1:2)==i);
            atom(i).bond.type = 1;
            atom(i).bond.index = Bond_index(A,1:2);
            atom(i).bond.dist = Bond_index(A,3);
            if ismember(i,Angle_index(:,1:3))
                [C,D]=find(Angle_index(:,1:3)==i);
                atom(i).angle.type = 1;
                atom(i).angle.index = Angle_index(C,1:3);
                atom(i).angle.angle = Angle_index(C,4);
                atom(i).angle.vec1 = Angle_index(C,5:7);
                atom(i).angle.vec2 = Angle_index(C,8:10);
            end
        end
        if mod(i,100)==1
            i-1
        end
    end
    
    %%%%%%%%%%
    
    Atom_labels=unique([atom.type]);
    for i=1:length(Atom_labels)
        label_ind=find(strcmpi([atom.type],Atom_labels(i)));
        Tot_dist=[];Tot_type=[];Tot_index=[];Tot_angleindex=[];Tot_bondindex=[];Tot_neighindex=[];Tot_coords=[];Tot_bonds=[];Tot_angles=[];
        for j=label_ind
            if numel([atom(j).neigh])>0
                Tot_index=[Tot_index; repmat(j,numel([atom(j).neigh.index]),1)];
                Tot_dist=[Tot_dist; [atom(j).neigh.dist]];
                Tot_type=[Tot_type; [atom(j).neigh.type]];
                Tot_neighindex=[Tot_neighindex; [atom(j).neigh.index]];
                Tot_coords=[Tot_coords; [atom(j).neigh.coords]];
            end
            
            if numel([atom(j).bond])>0
                Tot_bondindex=[Tot_bondindex; [atom(j).bond.index]];
                Tot_bonds=[Tot_bonds; [atom(j).bond.dist]];
            end
            if numel([atom(j).angle])>0
                Tot_angleindex=[Tot_angleindex; [atom(j).angle.index]];
                Tot_angles=[Tot_angles; [atom(j).angle.angle]];
            end
        end
        assignin('caller',strcat(char(Atom_labels(i)),'_dist')',[num2cell(Tot_index) num2cell(Tot_neighindex) Tot_type num2cell(Tot_dist)]);
        assignin('caller',strcat(char(Atom_labels(i)),'_coords')',[[atom(Tot_neighindex).x]' [atom(Tot_neighindex).y]' [atom(Tot_neighindex).z]']);
        assignin('caller',strcat(char(Atom_labels(i)),'_bonds')',[Tot_bondindex Tot_bonds]);
        assignin('caller',strcat(char(Atom_labels(i)),'_angles')',[Tot_angleindex Tot_angles]);
        assignin('caller',strcat(char(Atom_labels(i)),'_atom')',atom(ismember([atom.type],Atom_labels(i))));
    end
    assignin('caller',char(varargin{4}),atom);
end

if nAtoms < 20000
%     assignin('caller','dist_matrix',dist_matrix);
    assignin('caller','X_dist',(X_dist)');
    assignin('caller','Y_dist',(Y_dist)');
    assignin('caller','Z_dist',(Z_dist)');
end


%%%%%%%%%%%%% End of main function %%%%%%%%%%%%%%

    function mi_n = minimum_image(n, N)
        mi_n = n - N .* sign(n) .* floor((abs(n) + 0.5*N) ./ N);
    end

    function [pair, dist, num, rx, ry, rz] = calcpair_atom_pbc(XYZ_data, rmax, Box_dim)
        if numel(Box_dim)==3
            rx = bsxfun(@minus,XYZ_data(:,1),XYZ_data(:,1)');
            rx = rx - Box_dim(1)*round(rx./Box_dim(1));
            ry = bsxfun(@minus,XYZ_data(:,2),XYZ_data(:,2)');
            ry = ry - Box_dim(2)*round(ry./Box_dim(2));
            rz = bsxfun(@minus,XYZ_data(:,3),XYZ_data(:,3)');
            rz = rz - Box_dim(3)*round(rz./Box_dim(3));
        else %if numel(Box_dim)==9 && numel(nonzeros(Box_dim(1,4:end))) > 0 % Untested
            rx = bsxfun(@minus,XYZ_data(:,1),XYZ_data(:,1)');
            ry = bsxfun(@minus,XYZ_data(:,2),XYZ_data(:,2)');
            rz = bsxfun(@minus,XYZ_data(:,3),XYZ_data(:,3)');
            
            rx = rx - xz*round(rz./Box_dim(3));
            ry = ry - yz*round(rz./Box_dim(3));
            rz = rz - Box_dim(3)*round(rz./Box_dim(3));
            
            rx = rx - xy*round(ry./Box_dim(2));
            ry = ry - Box_dim(2)*round(ry./Box_dim(2));
            
            rx = rx - Box_dim(1)*round(rx./Box_dim(1));
        end
        dist = sqrt(rx.^2 + ry.^2 + rz.^2);
        index = find(tril(dist < rmax, -1));
        if isempty(index)
            pair = [];
            dist = [];
            num = 0;
            rx = [];
            ry = [];
            rz = [];
        else
            [pair1, pair2] = ind2sub(size(dist), index);
            pair = zeros(numel(pair1), 2);
            pair(:,1) = pair1;
            pair(:,2) = pair2;
            dist = dist(index);
            num = numel(dist);
        end
    end

    function [pair, dist, num, rx, ry, rz] = calcpair2_atom_pbc(XYZ_data1, XYZ_data2, rmax, Box_dim)
        if numel(Box_dim)==3
            rx = bsxfun(@minus,XYZ_data1(:,1),XYZ_data2(:,1)');
            rx = rx - Box_dim(1)*round(rx./Box_dim(1));
            ry = bsxfun(@minus,XYZ_data1(:,2),XYZ_data2(:,2)');
            ry = ry - Box_dim(2)*round(ry./Box_dim(2));
            rz = bsxfun(@minus,XYZ_data1(:,3),XYZ_data2(:,3)');
            rz = rz - Box_dim(3)*round(rz./Box_dim(3));
        else %if numel(Box_dim)==9 && numel(nonzeros(Box_dim(1,4:end))) > 0 % Untested
            rx = bsxfun(@minus,XYZ_data1(:,1),XYZ_data2(:,1)');
            ry = bsxfun(@minus,XYZ_data1(:,2),XYZ_data2(:,2)');
            rz = bsxfun(@minus,XYZ_data1(:,3),XYZ_data2(:,3)');
            
            rx = rx - xz*round(rz./Box_dim(3));
            ry = ry - yz*round(rz./Box_dim(3));
            rz = rz - Box_dim(3)*round(rz./Box_dim(3));
            
            rx = rx - xy*round(ry./Box_dim(2));
            ry = ry - Box_dim(2)*round(ry./Box_dim(2));
            
            rx = rx - Box_dim(1)*round(rx./Box_dim(1));
        end
        dist = sqrt(rx.^2 + ry.^2 + rz.^2);
        index = find(dist < rmax);
        if isempty(index)
            pair = [];
            dist = [];
            num = 0;
            rx = [];
            ry = [];
            rz = [];
        else
            [pair1, pair2] = ind2sub(size(dist), index);
            pair = zeros(numel(pair1), 2);
            pair(:,1) = pair1;
            pair(:,2) = pair2;
            dist = dist(index);
            num = numel(dist);
        end
    end
end % end function

