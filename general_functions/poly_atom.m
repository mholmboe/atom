%% poly_atom.m
% * This function tries to plot pretty polyhedras, similar to the show_atom
% * function
%
%% Version
% 2.082
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=poly_atom(atom,Box_dim)                         % Plots polyhedras for all atomtypes except O and H,
% # atom=poly_atom(atom,Box_dim,{'Al' 'Mgo' 'Si'})       % Plots polyhedras for the stated atomtypes
% # atom=poly_atom(atom,Box_dim,0,2.25)                  % 2.25 is a bond cutoff (sort of)
% # atom=poly_atom(atom,Box_dim,0,2.25,1)                % Last argument (1/0) to plot the Box
% # atom=poly_atom(atom,Box_dim,0,2.25,1,0.3)            % Last argument sets the transperancy

function atom = poly_atom(atom,Box_dim,varargin)

if nargin>2
    polytype=varargin{1}; % Dummy value
else
    polytype=unique([atom.type]);
    polytype(strncmpi(polytype,'H',1))=[];
    polytype(strncmpi(polytype,'O',1))=[];
end

if isnumeric(polytype)
    polytype=unique([atom.type]);
    polytype(strncmpi(polytype,'H',1))=[];
    polytype(strncmpi(polytype,'O',1))=[];
end

polytype_ind=find(ismember([atom.type],polytype));

if nargin>3
    rmaxlong=varargin{2}; % Dummy value
else
    rmaxlong=2.25;
end

if nargin>5
    alpha=1-varargin{4};
else
    alpha=1; % Transperacy
end

atom=element_atom(atom);
XYZ_labels=[atom.type]';
nAtoms = size(XYZ_labels,1);

% if strncmpi(style,'crystal',3)
%     radii = 1/4*abs(radius_crystal(XYZ_labels));
% else
%     radii = 1/4*abs(radius_vdw(XYZ_labels));
% end
color =  1*element_color(XYZ_labels);

if nargin>6
    trans_vec=varargin{5};
    if numel(trans_vec)==3
        atom=translate_atom(atom,trans_vec(1:3));
        if numel(trans_vec)==4
            atom=wrap_atom(atom,Box_dim);
        end
    end
end

if nargin>7
    color=varargin{6};
    color=repmat(color,nAtoms,1);
end

% XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];

assignin('caller','color',color)

% if nargin>1

%    Box_dim=varargin{1};

if numel(Box_dim)>0
    
    if numel(Box_dim)==1
        Box_dim(1)=Box_dim(1);
        Box_dim(2)=Box_dim(1);
        Box_dim(3)=Box_dim(1);
    end
    
    atom = bond_atom(atom,Box_dim,rmaxlong,.6);
    
    if max([atom.molid])<max([atom.index])
        i=1;
        while numel(Bond_index)==0 && i < 10
            atom = bond_atom(atom,Box_dim,rmaxlong+i/4,.6+i/10);
            i=i+1;
        end
    end
    
end

% Sets plot limits for the data
xlo = floor(min([-5 min([atom.x])-2])); xhi = ceil(max([max([atom.x])+2 Box_dim(1)])/5)*5;
ylo = floor(min([-5 min([atom.y])-2])); yhi = ceil(max([max([atom.y])+2 Box_dim(2)])/5)*5;
zlo = floor(min([-5 min([atom.z])-2])); zhi = ceil(max([max([atom.z])+2 Box_dim(3)])/5)*5;
% else
%     xlo = floor(min([-5 min([atom.x])-max(radii)])); xhi = ceil(max(max([atom.x]))/5)*5;
%     ylo = floor(min([-5 min([atom.y])-max(radii)])); yhi = ceil(max(max([atom.y]))/5)*5;
%     zlo = floor(min([-5 min([atom.z])-max(radii)])); zhi = ceil(max(max([atom.z]))/5)*5;
% end

if xhi < 5
    xhi=5;
end

if yhi < 5
    yhi=5;
end

if zhi < 5
    zhi=5;
end

% xhi=35
% yhi=40
% zhi=45
hold on;
cameratoolbar
rotate3d on;
camlight(220,210,'infinite');
set(gcf,'Visible','on','Color',[1 1 1]);
set(gca,'Color',[1 1 1],'PlotBoxAspectRatio',[(xhi-xlo)/(zhi-zlo) (yhi-ylo)/(zhi-zlo) (zhi-zlo)/(zhi-zlo)],'FontSize',24); %

fig = gcf;
fig.Color = [1 1 1];
ax = fig.CurrentAxes;
ax.XLim = [xlo xhi];
ax.YLim = [ylo yhi];
ax.ZLim = [zlo zhi];
xlabel('X [Å]'); ylabel('Y [Å]'); zlabel('Z [Å]');
view([0,0]);

for ip=1:numel(polytype_ind)
    
    i=polytype_ind(ip);
    
    color_temp = color(i,:);
    
    if numel(atom(i).neigh.index)>0
        
        r_vec = [[atom(i).neigh.r_vec(:,1)] [atom(i).neigh.r_vec(:,2)] [atom(i).neigh.r_vec(:,3)]];
        
        dist=[];
        for j=1:size(r_vec,1)
            dist(j,:)=([ (r_vec(j,1)-r_vec(:,1)).^2 + (r_vec(j,2)-r_vec(:,2)).^2 + ([r_vec(j,3)-r_vec(:,3)]).^2] ).^.5;
        end
        
        dist(dist>1.5*rmaxlong)=0;
        
        PolyInd=[];n=0;
        for ik=1:size(r_vec,1)
            for il=ik+1:size(r_vec,1)
                for im=il+1:size(r_vec,1)
                    if dist(il,ik)>0 && dist(im,ik)>0 && dist(im,il)>0
                        n=n+1;
                        PolyInd(n,:)=[ik il im];
                    end
                end
            end
        end
        
        poly_color=color_temp.^.5;
        poly_color_edge=color_temp./5;
        patch('Faces',PolyInd,'Vertices',[atom(i).x atom(i).y atom(i).z]+[atom(i).neigh.r_vec],'FaceColor',poly_color,...
            'EdgeColor',poly_color_edge,'FaceLighting','gouraud','Facealpha',alpha,...
            'AmbientStrength',0.8,'DiffuseStrength',0.6,'SpecularStrength',0.2,'LineWidth',1);
    end
    
    if mod(i,1000)==1 || i==size(atom,2)
        if i > 1
            i-1
            drawnow limitrate
        end
    end
end


if nargin>4
    if varargin{3}>0
        Simbox = draw_box_atom(Box_dim,[0 0 0.8],2);
    end
end

hold off;

