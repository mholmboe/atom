%% wrap_atom.m
% * This function wraps the atoms into the box. It can also wrap only along
% the x and y dimensions, ie neglecting wrapping along the z-direction
% * Which one is fastest? Ortogonal or triclinic version?
%
%% Version
% 2.11
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = wrap_atom(atom,Box_dim)
% # atom = wrap_atom(atom,Box_dim,'xy')
%
function atom = wrap_atom(atom,Box_dim,varargin)

attribut_in=fieldnames(atom);

if nargin==2
    if size(Box_dim(1,:),2)==3
        % disp('assuming orthogonal box when wrapping!!!')
%         ind_hiz=find([atom.z]>=Box_dim(3));
%         z_shift=num2cell([[atom(ind_hiz).z]-Box_dim(3)]');
%         [atom((ind_hiz)).z]=deal(z_shift{:});
%         ind_loz=find([atom.z]<0);
%         z_shift=num2cell([[atom(ind_loz).z]+Box_dim(3)]');
%         [atom((ind_loz)).z]=deal(z_shift{:});
%         
%         ind_hiy=find([atom.y]>=Box_dim(2));
%         y_shift=num2cell([[atom(ind_hiy).y]-Box_dim(2)]');
%         [atom((ind_hiy)).y]=deal(y_shift{:});
%         ind_loy=find([atom.y]<0);
%         y_shift=num2cell([[atom(ind_loy).y]+Box_dim(2)]');
%         [atom((ind_loy)).y]=deal(y_shift{:});
%         
%         ind_hix=find([atom.x]>=Box_dim(1));
%         x_shift=num2cell([[atom(ind_hix).x]-Box_dim(1)]');
%         [atom((ind_hix)).x]=deal(x_shift{:});
%         ind_lox=find([atom.x]<0);
%         x_shift=num2cell([[atom(ind_lox).x]+Box_dim(1)]');
%         [atom((ind_lox)).x]=deal(x_shift{:});
     
        X_data = num2cell([atom.x]' - Box_dim(1)*floor([atom.x]'./Box_dim(1)));
        Y_data = num2cell([atom.y]' - Box_dim(2)*floor([atom.y]'./Box_dim(2)));
        Z_data = num2cell([atom.z]' - Box_dim(3)*floor([atom.z]'./Box_dim(3)));
        [atom.x]=deal(X_data{:});
        [atom.y]=deal(Y_data{:});
        [atom.z]=deal(Z_data{:});  
        
    else
 %       disp('will try to wrap triclinic Box_dim!!!')
        xy=Box_dim(6); xz=Box_dim(8); yz=Box_dim(9);
        orto=orto_atom(atom,Box_dim);
        X_data = num2cell([orto.xfrac]' - floor([orto.xfrac]'));
        Y_data = num2cell([orto.yfrac]' - floor([orto.yfrac]'));
        Z_data = num2cell([orto.zfrac]' - floor([orto.zfrac]'));
        [orto.x]=deal(X_data{:});
        [orto.y]=deal(Y_data{:});
        [orto.z]=deal(Z_data{:});
        orto = scale_atom(orto,[1 1 1],orto_Box_dim,'ALL');
        atom = triclinic_atom(orto,orto_Box_dim,[xy xz yz],'tilt');
    end
    
elseif nargin == 3
    if size(Box_dim(1,:),2)==3
 %       disp('assuming orthogonal box when wrapping!!!')
        % ind_hiz=find([atom.z]>=Box_dim(3));
        % z_shift=num2cell([[atom(ind_hiz).z]-Box_dim(3)]');
        % [atom((ind_hiz)).z]=deal(z_shift{:});
        % ind_loz=find([atom.z]<0);
        % z_shift=num2cell([[atom(ind_loz).z]+Box_dim(3)]');
        % [atom((ind_loz)).z]=deal(z_shift{:});
        
%         ind_hiy=find([atom.y]>=Box_dim(2));
%         y_shift=num2cell([[atom(ind_hiy).y]-Box_dim(2)]');
%         [atom((ind_hiy)).y]=deal(y_shift{:});
%         ind_loy=find([atom.y]<0);
%         y_shift=num2cell([[atom(ind_loy).y]+Box_dim(2)]');
%         [atom((ind_loy)).y]=deal(y_shift{:});
%         
%         ind_hix=find([atom.x]>=Box_dim(1));
%         x_shift=num2cell([[atom(ind_hix).x]-Box_dim(1)]');
%         [atom((ind_hix)).x]=deal(x_shift{:});
%         ind_lox=find([atom.x]<0);
%         x_shift=num2cell([[atom(ind_lox).x]+Box_dim(1)]');
%         [atom((ind_lox)).x]=deal(x_shift{:});
        
        X_data = num2cell([atom.x]' - Box_dim(1)*floor([atom.x]'./Box_dim(1)));
        Y_data = num2cell([atom.y]' - Box_dim(2)*floor([atom.y]'./Box_dim(2)));
%         Z_data = num2cell([atom.z]' - Box_dim(3)*floor([atom.z]'./Box_dim(3)));
        [atom.x]=deal(X_data{:});
        [atom.y]=deal(Y_data{:});
%         [atom.z]=deal(Z_data{:});  
        
    else
 %       disp('will try to wrap triclinic Box_dim!!!')
        xy=Box_dim(6); xz=Box_dim(8); yz=Box_dim(9);
        orto=orto_atom(atom,Box_dim);
        X_data = num2cell([orto.xfrac]' - floor([orto.xfrac]'));
        Y_data = num2cell([orto.yfrac]' - floor([orto.yfrac]'));
        Z_data = num2cell([orto.zfrac]'); % - floor([orto.zfrac]'));
        [orto.x]=deal(X_data{:});
        [orto.y]=deal(Y_data{:});
        [orto.z]=deal(Z_data{:});
        orto = scale_atom(orto,[1 1 1],orto_Box_dim,'ALL');
        atom = triclinic_atom(orto,orto_Box_dim,[xy xz yz],'tilt');
    end
end

try
    attribut_out=fieldnames(atom);
    different_ind=find(~ismember(attribut_out,attribut_in));
    atom=rmfield(atom,attribut_out(different_ind));
catch
    
end

%assignin('caller','atom',atom);