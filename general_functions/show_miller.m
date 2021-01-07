%%  show_miller.m
% * This function draws up the h,k,l Miller plane within the Box_dim/Cell
% * This function is inspired by miller.m, written by Lucas Ennes:
% * Miller index (https://www.mathworks.com/matlabcentral/fileexchange/68887-miller-index)

%% Version
% 2.09
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # show_miller(h,k,l) % h k l designates the Miller plane of a ortohogonal 1x1x1 cell
% # show_miller(h,k,l,Box_dim) % the 1x3 or the 1x9 Box_dim variable
% # show_miller(h,k,l,Box_dim,[0 0 1]) % Color
% # show_miller(h,k,l,Box_dim,[0 0 1],0.5) % Transperancy

function show_miller(h,k,l,varargin)

if nargin==3
    Box_dim=[1 1 1];
    color='b';
elseif nargin==4
    Box_dim=varargin{1};
    color='b';
elseif nargin>4
    Box_dim=varargin{1};
    color=varargin{2};
end

if nargin>5
    alpha=1-varargin{3};
else
    alpha=.5; % Transperacy
end

if size(Box_dim,2) == 1
    Box_dim(2)=Box_dim(1);
    Box_dim(3)=Box_dim(2);
end

%% Set the unit cell parameters
if size(Box_dim,2) == 9
    lx=Box_dim(1);
    ly=Box_dim(2);
    lz=Box_dim(3);
    xy=Box_dim(6);
    xz=Box_dim(8);
    yz=Box_dim(9);
    a=lx;
    b=(ly^2+xy^2)^.5;
    c=(lz^2+xz^2+yz^2)^.5;
    alfa=rad2deg(acos((ly*yz+xy*xz)/(b*c)));
    beta=rad2deg(acos(xz/c));
    gamma=rad2deg(acos(xy/b));
    %% Set the unit cell parameters
elseif size(Box_dim,2) == 6
    a = Box_dim(1);
    b = Box_dim(2);
    c = Box_dim(3);
    alfa = Box_dim(4);
    beta = Box_dim(5);
    gamma = Box_dim(6);
    lx = a;
    xy = b * cos(deg2rad(gamma));
    ly = (b^2-xy^2)^.5;
    xz = c*cos(deg2rad(beta));
    yz = (b*c*cos(deg2rad(alfa))-xy*xz)/ly;
    lz = (c^2 - xz^2 - yz^2)^0.5;
    Box_dim=[lx ly lz 0 0 xy 0 xz yz];
    Box_dim(Box_dim<0.00001&Box_dim>-0.00001)=0;
    if sum(find(Box_dim(4:end)))<0.0001
        Box_dim=Box_dim(1:3);
    end
elseif size(Box_dim,2) == 3
    lx=Box_dim(1);
    ly=Box_dim(2);
    lz=Box_dim(3);
    xy=0;
    xz=0;
    yz=0;
    a=lx;
    b=ly;
    c=lz;
    alfa=90;
    beta=90;
    gamma=90;
end

%% Converting to radians
alfa_rad=alfa*pi/180;
beta_rad=beta*pi/180;
gamma_rad=gamma*pi/180;

%% Determine the two theta values for each h,k,l
Volume=a*b*c*(1-cos(alfa_rad)^2-cos(beta_rad)^2-cos(gamma_rad)^2+2*cos(alfa_rad)*cos(beta_rad)*cos(gamma_rad))^0.5;
recip_d_hkl=1/Volume.*...
    (h.^2*b^2*c^2*sin(alfa_rad)^2+...
    k.^2*a^2*c^2*sin(beta_rad)^2+...
    l.^2*a^2*b^2*sin(gamma_rad)^2+...
    2*h.*k*a*b*c^2*(cos(alfa_rad)*cos(beta_rad)-cos(gamma_rad))+...
    2*k.*l*a^2*b*c*(cos(beta_rad)*cos(gamma_rad)-cos(alfa_rad))+...
    2*h.*l*a*b^2*c*(cos(alfa_rad)*cos(gamma_rad)-cos(beta_rad))).^(0.5);
recip_d_hkl=real(recip_d_hkl);

d_hkl=1./recip_d_hkl;

%%

H=[];
K=[];
L=[];
%% Set the plane
if h ~= 0 && k ~= 0 && l ~= 0
    
    % Triangular plane
    h = 1 / h;    k = 1 / k;    l = 1 / l;
    
    %     high = max([abs(h),abs(k),abs(l)]);
    %
    %     h = h / high;
    %     k = k / high;
    %     l = l / high;
    
    subh = 0;    subk = 0;    subl = 0;
    if h < 0
        subh = - h;
    end
    if k < 0
        subk = - k;
    end
    if l < 0
        subl = - l;
    end
    H = [h+subh, subh, subh];
    K = [subk, k+subk, subk];
    L = [subl, subl, l+subl];
    
    if h < 0
        H(H>0)=1;
        H(H==0)=1+h;
    end
    if k < 0
        K(K>0)=1;
        K(K==0)=1+k;
    end
    if l < 0
        L(L>0)=1;
        L(L==0)=1+l;
    end
    
    
else
    % Rectangular plane
    % Count zeros
    zeros = 0;
    if h == 0
        zeros = zeros + 1;
    end
    if k == 0
        zeros = zeros + 1;
    end
    if l == 0
        zeros = zeros + 1;
    end
    % 2 zeros
    if zeros == 2
        if h ~= 0
            for i=1:abs(h)
                H = [H;i/h*[1 1 1 1]];
                K = [K;    [0 0 1 1]];
                L = [L;    [0 1 1 0]];
                if h < 0
                    H(H<0)=H(H<0)+1;
                end
            end
        elseif k ~= 0
            for i=1:abs(k)
                H = [H;    [0 0 1 1]];
                K = [K;i/k*[1 1 1 1]];
                L = [L;    [0 1 1 0]];
                if k < 0
                    K(K<0)=K(K<0)+1;
                end
            end
        elseif l ~= 0
            for i=1:abs(l)
                H = [H;    [0 1 1 0]];
                K = [K;    [0 0 1 1]];
                L = [L;i/l*[1 1 1 1]];
                if l < 0
                    L(L<0)=L(L<0)+1;
                end
            end
        end
        % 1 zero
    elseif zeros == 1
        
        subh = 0;        subk = 0;        subl = 0;
        
        if h == 0
            k = (1/k);
            l = (1/l);
            
            %             high = max([k, l]);
            %             k = k/high;
            %             l = l/high;
            
            if k < 0
                subk = -k;
            end
            if l < 0
                subl = -l;
            end
            %
            H = [0 1 1 0];
            K = [subk subk k+subk k+subk];
            L = [l+subl l+subl subl subl];
            
        end
        if k == 0
            h = (1/h);
            l = (1/l);
            
            %             high = max([h, l]);
            %             h = h/high;
            %             l = l/high;
            
            if h < 0
                subh = -h;
            end
            if l < 0
                subl = -l;
            end
            K = [0 1 1 0];
            H = [h+subh h+subh subh subh];
            L = [subl subl l+subl l+subl];
        end
        if l == 0
            k = (1/k);
            h = (1/h);
            
            %             high = max([k, h]);
            %             k = k/high;
            %             h = h/high;
            
            if k < 0
                subk = -k;
            end
            if h < 0
                subh = -h;
            end
            L = [0 1 1 0];
            K = [subk subk k+subk k+subk];
            H = [h+subh h+subh subh subh];
        end
        
        if h < 0
            H(H>0)=1;
            H(H==0)=1+h;
        end
        if k < 0
            K(K>0)=1;
            K(K==0)=1+k;
        end
        if l < 0
            L(L>0)=1;
            L(L==0)=1+l;
        end
        
    else
        error('Weird plane')
    end
end

Simbox = [a b c alfa beta gamma Volume];
assignin('base','Simbox',Simbox);

v=(1 - cos(deg2rad(alfa))^2 - cos(deg2rad(beta))^2 - cos(deg2rad(gamma))^2 + 2*cos(deg2rad(alfa))*cos(deg2rad(beta))*cos(deg2rad(gamma)))^.5;

FromFrac=[a b*cos(deg2rad(gamma))  c*cos(deg2rad(beta));...
    0 b*sin(deg2rad(gamma))  c*(cos(deg2rad(alfa))-cos(deg2rad(beta))*cos(deg2rad(gamma)))/sin(deg2rad(gamma));...
    0 0                      c*v/sin(deg2rad(gamma))];

max_order=max([size(H,1) size(K,1) size(L,1)]);
max_dim=max([size(H,2) size(K,2) size(L,2)]);
for i=1:max_order
    for j=1:max_dim
        HKL=FromFrac*[H(i,j) K(i,j) L(i,j)]';
        H(i,j)=HKL(1);
        K(i,j)=HKL(2);
        L(i,j)=HKL(3);
    end
end

hold on;
cameratoolbar
rotate3d on;
camlight(220,210,'infinite');
set(gcf,'Visible','on','Color',[1 1 1]);
set(gca,'Color',[1 1 1],'FontSize',24); % 'PlotBoxAspectRatio',[(xhi-xlo)/(zhi-zlo) (yhi-ylo)/(zhi-zlo) (zhi-zlo)/(zhi-zlo)]

fig = gcf;
fig.Color = [1 1 1];
xlabel('a [Å]'); ylabel('b [Å]'); zlabel('c [Å]');
% view([0,0]);

% if nargin>5
%     if varargin{3}==0
draw_box_atom(Box_dim,[0 0 0.8],2);
%     end
% end

hold on;
assignin('caller','H',H);
assignin('caller','K',K);
assignin('caller','L',L);

for i=1:size(H,1)
    patch(H(i,:),K(i,:),L(i,:),color,'FaceAlpha',alpha)
end
view(3)
% xlim([0 1]);ylim([0 1]);zlim([0 1])

