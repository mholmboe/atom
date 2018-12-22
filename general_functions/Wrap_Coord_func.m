%% Wrap_Coord_func.m 
% * This is an old function that wraps atoms 'sticking out' back into the box.
% * Important notice: Untested for triclinic boxes... use wrap_atom instead
%
%% Similar
% wrap_atom
%
%% Version
% 2.0
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # XYZ_data = Wrap_Coord_func(XYZ_data,Box_dim)

function XYZ_data = Wrap_Coord_func(XYZ_data,Box_dim)
%% 

if size(Box_dim,2) > 3
    Lx = Box_dim(1);
    Ly = Box_dim(2);
    Lz = Box_dim(3);
    xy = Box_dim(6);
    xz = Box_dim(8);
    yz = Box_dim(9);
    for i = 1:3:size(XYZ_data,2)
        for j = 1:size(XYZ_data,1) %    for j = 1:size(Atomcoord,1)-1;
            
            if XYZ_data(j,i+2) > Lz
                XYZ_data(j,i+2) = XYZ_data(j,i+2) - Lz;
                XYZ_data(j,i) = XYZ_data(j,i) - xz;
                XYZ_data(j,i+1) = XYZ_data(j,i+1) - yz;
            elseif XYZ_data(j,i+2) < 0
                XYZ_data(j,i+2) = XYZ_data(j,i+2) + Lz;
                XYZ_data(j,i) = XYZ_data(j,i) + xz;
                XYZ_data(j,i+1) = XYZ_data(j,i+1) + yz;
            end
            
            if XYZ_data(j,i+1) > Ly
                XYZ_data(j,i+1) = XYZ_data(j,i+1) - Ly;
                XYZ_data(j,i) = XYZ_data(j,i) - xy;
            elseif XYZ_data(j,i+1) < 0
                XYZ_data(j,i+1) = XYZ_data(j,i+1) + Ly;
                XYZ_data(j,i) = XYZ_data(j,i) + xy;
            end
            
            if XYZ_data(j,i) > Lx
                XYZ_data(j,i) = XYZ_data(j,i) - Lx;
            elseif XYZ_data(j,i) < 0
                XYZ_data(j,i) = XYZ_data(j,i) + Lx;
            end
            
        end
    end
    
    
elseif length(Box_dim) == 3
%     Lx = Box_dim(1);
%     Ly = Box_dim(2);
%     Lz = Box_dim(3);
%     
%     for i = 1:3:size(XYZ_data,2);
%         for j = 1:size(XYZ_data,1); %    for j = 1:size(Atomcoord,1)-1;
%             
%             if XYZ_data(j,i+2) > Lz;
%                 XYZ_data(j,i+2) = XYZ_data(j,i+2) - Lz;
%                
%             elseif XYZ_data(j,i+2) < 0;
%                 XYZ_data(j,i+2) = XYZ_data(j,i+2) + Lz;
%                
%             end
%             
%             if XYZ_data(j,i+1) > Ly;
%                 XYZ_data(j,i+1) = XYZ_data(j,i+1) - Ly;
%                
%             elseif XYZ_data(j,i+1) < 0;
%                 XYZ_data(j,i+1) = XYZ_data(j,i+1) + Ly;
%                 
%             end
%             
%             if XYZ_data(j,i) > Lx;
%                 XYZ_data(j,i) = XYZ_data(j,i) - Lx;
%             elseif XYZ_data(j,i) < 0;
%                 XYZ_data(j,i) = XYZ_data(j,i) + Lx;
%             end
%             
%         end
%     end
    %% Does this work? 
    for i = 1:3:size(XYZ_data,2)
        indxlo=find(XYZ_data(:,i)<0);
        XYZ_data(indxlo,i)=XYZ_data(indxlo,i)+Box_dim(1);
        indxhi=find(XYZ_data(:,i)>Box_dim(1));
        XYZ_data(indxhi,i)=XYZ_data(indxhi,i)-Box_dim(1);
        
        indylo=find(XYZ_data(:,i+1)<0);
        XYZ_data(indylo,i+1)=XYZ_data(indylo,i+1)+Box_dim(2);
        indyhi=find(XYZ_data(:,i+1)>Box_dim(2));
        XYZ_data(indyhi,i+1)=XYZ_data(indyhi,i+1)-Box_dim(2);
        
        indzlo=find(XYZ_data(:,i+2)<0);
        XYZ_data(indzlo,i+2)=XYZ_data(indzlo,i+2)+Box_dim(3);
        indzhi=find(XYZ_data(:,i+2)>Box_dim(3));
        XYZ_data(indzhi,i+2)=XYZ_data(indzhi,i+2)-Box_dim(3);
    end
else
    indlo=find(XYZ_data<0);
    XYZ_data(indlo)=XYZ_data(indlo)+Box_dim;
    indhi=find(XYZ_data(:)>Box_dim);
    XYZ_data(indhi)=XYZ_data(indhi)-Box_dim;
end


