%% COM_func.m
% * This super  old function calculates the center of mass for water. Slow due to pbc...
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = COM_atom(atom,Box_dim)
%
function AtomCoords_COM = COM_SOL(MolID,XYZ_data,Atom_label,XYZ_labels,Box_dim)

All_index = MolID{1};

Lx = Box_dim(1);
Ly = Box_dim(2);
Lz = Box_dim(3);
xy = Box_dim(6);
xz = Box_dim(8);
yz = Box_dim(9);

mO = 15.99941;
mH = 1.00794;
mH2O = mO + 2*mH;

AtomCoords_COM = zeros(size(XYZ_data,1),size(XYZ_data,2)/3);

if strncmp(XYZ_labels(All_index(1)+2),Atom_label{1},1)==1
    for i = 1:9:size(XYZ_data,2)
        for j = 1:size(XYZ_data,1)
            
            % H_1
            if (XYZ_data(j,i+2+3)-XYZ_data(j,i+2)) > Lz/2
                XYZ_data(j,i+2+3) = XYZ_data(j,i+2+3) - Lz;
                XYZ_data(j,i+3) = XYZ_data(j,i+3) - xz;
                XYZ_data(j,i+1+3) = XYZ_data(j,i+1+3) - yz;
            elseif (XYZ_data(j,i+2+3)-XYZ_data(j,i+2)) < -Lz/2
                XYZ_data(j,i+2+3) = XYZ_data(j,i+2+3) + Lz;
                XYZ_data(j,i+3) = XYZ_data(j,i+3) + xz;
                XYZ_data(j,i+1+3) = XYZ_data(j,i+1+3) + yz;
            end
            
            if (XYZ_data(j,i+1+3)-XYZ_data(j,i+1)) > Ly/2
                XYZ_data(j,i+1+3) = XYZ_data(j,i+1+3) - Ly;
                XYZ_data(j,i+3) = XYZ_data(j,i+3) - xy;
            elseif (XYZ_data(j,i+1+3)-XYZ_data(j,i+1)) < -Ly/2
                XYZ_data(j,i+1+3) = XYZ_data(j,i+1+3) + Ly;
                XYZ_data(j,i+3) = XYZ_data(j,i+3) + xy;
            end
            
            if (XYZ_data(j,i+3)-XYZ_data(j,i)) > Lx/2
                XYZ_data(j,i+3) = XYZ_data(j,i+3) - Lx;
            elseif (XYZ_data(j,i+3)-XYZ_data(j,i)) < -Lx/2
                XYZ_data(j,i+3) = XYZ_data(j,i+3) + Lx;
            end
            
            % H_2
            if (XYZ_data(j,i+2+6)-XYZ_data(j,i+2)) > Lz/2
                XYZ_data(j,i+2+6) = XYZ_data(j,i+2+6) - Lz;
                XYZ_data(j,i+6) = XYZ_data(j,i+6) - xz;
                XYZ_data(j,i+1+6) = XYZ_data(j,i+1+6) - yz;
            elseif (XYZ_data(j,i+2+6)-XYZ_data(j,i+2)) < -Lz/2
                XYZ_data(j,i+2+6) = XYZ_data(j,i+2+6) + Lz;
                XYZ_data(j,i+6) = XYZ_data(j,i+6) + xz;
                XYZ_data(j,i+1+6) = XYZ_data(j,i+1+6) + yz;
            end
            
            if (XYZ_data(j,i+1+6)-XYZ_data(j,i+1)) > Ly/2
                XYZ_data(j,i+1+6) = XYZ_data(j,i+1+6) - Ly;
                XYZ_data(j,i+6) = XYZ_data(j,i+6) - xy;
            elseif (XYZ_data(j,i+1+6)-XYZ_data(j,i+1)) < -Ly/2
                XYZ_data(j,i+1+6) = XYZ_data(j,i+1+6) + Ly;
                XYZ_data(j,i+6) = XYZ_data(j,i+6) + xy;
            end
            
            if (XYZ_data(j,i+6)-XYZ_data(j,i)) > Lx/2
                XYZ_data(j,i+6) = XYZ_data(j,i+6) - Lx;
            elseif (XYZ_data(j,i+6)-XYZ_data(j,i)) < -Lx/2
                XYZ_data(j,i+6) = XYZ_data(j,i+6) + Lx;
            end
            
            
        end
        
    end
    
    
    for k=1:3;
        AtomCoords_COM(:,k:3:end) = (mO/mH2O*XYZ_data(:,k:9:end) + mH/mH2O*XYZ_data(:,k+3:9:end) + mH/mH2O*XYZ_data(:,k+6:9:end));
        k
    end
else
    disp('Wrong index for O and H?')
end
