%% round_atom.m
% * This function rounds the coordinates in the atom struct
%
%% Version
% 2.08
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = round_atom(atom,Box_dim,varargin)
%
function atom = round_atom(atom,Box_dim,varargin)

disp('Rounding the coordinates')

if nargin>2
    precision=varargin{1};
else
    precision=3;
end

if nargin>3 % will only round the fractional coordinates
    x_coord=num2cell(round([atom.xfrac],precision)); [atom.xfrac]=deal(x_coord{:});
    
    y_coord=num2cell(round([atom.yfrac],precision)); [atom.yfrac]=deal(y_coord{:});
    
    z_coord=num2cell(round([atom.zfrac],precision)); [atom.zfrac]=deal(z_coord{:});
else
    x_coord=num2cell(round([atom.x],precision)); [atom.x]=deal(x_coord{:});
    
    y_coord=num2cell(round([atom.y],precision)); [atom.y]=deal(y_coord{:});
    
    z_coord=num2cell(round([atom.z],precision)); [atom.z]=deal(z_coord{:});
    try
        x_coord=num2cell(round([atom.xfrac],precision)); [atom.xfrac]=deal(x_coord{:});
        
        y_coord=num2cell(round([atom.yfrac],precision)); [atom.yfrac]=deal(y_coord{:});
        
        z_coord=num2cell(round([atom.zfrac],precision)); [atom.zfrac]=deal(z_coord{:});
    catch
        disp('Found no fractional coordinates to round...')
    end
end



atom=update_atom(atom);
