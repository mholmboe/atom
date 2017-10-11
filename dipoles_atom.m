%% dipole_atom.m
% * This function calculates the dipole vector of water. Similar to the COM_SOL.
% * Is it correct?
% * Tested when?
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * dipole_vec = dipoles_atom(traj,Box_dim)


function dipole_vec = dipoles_atom(traj,Box_dim)
 

% All_index = MolID{1};

Lx = Box_dim(1);
Ly = Box_dim(2);
Lz = Box_dim(3);
xy = Box_dim(6);
xz = Box_dim(8);
yz = Box_dim(9);

% if strncmp(textdata(All_index(1)+2),Atom_label{1},1)==1;
for i = 1:9:size(traj,2);
    for j = 1:size(traj,1); %
        if isnan(traj(j,i))==0 % In case the data is set to NaN...
            % H_1
            if (traj(j,i+2+3)-traj(j,i+2)) > Lz/2;
                traj(j,i+2+3) = traj(j,i+2+3) - Lz;
                %                 Elements(j,i+3) = Elements(j,i+3) - xz;
                %                 Elements(j,i+1+3) = Elements(j,i+1+3) - yz;
            elseif (traj(j,i+2+3)-traj(j,i+2)) < -Lz/2;
                traj(j,i+2+3) = traj(j,i+2+3) + Lz;
                %                 Elements(j,i+3) = Elements(j,i+3) + xz;
                %                 Elements(j,i+1+3) = Elements(j,i+1+3) + yz;
            end
            
            if (traj(j,i+1+3)-traj(j,i+1)) > Ly/2;
                traj(j,i+1+3) = traj(j,i+1+3) - Ly;
                %                 Elements(j,i+3) = Elements(j,i+3) - xy;
            elseif (traj(j,i+1+3)-traj(j,i+1)) < -Ly/2;
                traj(j,i+1+3) = traj(j,i+1+3) + Ly;
                %                 Elements(j,i+3) = Elements(j,i+3) + xy;
            end
            
            if (traj(j,i+3)-traj(j,i)) > Lx/2;
                traj(j,i+3) = traj(j,i+3) - Lx;
            elseif (traj(j,i+3)-traj(j,i)) < -Lx/2;
                traj(j,i+3) = traj(j,i+3) + Lx;
            end
            
            % H_2
            if (traj(j,i+2+6)-traj(j,i+2)) > Lz/2;
                traj(j,i+2+6) = traj(j,i+2+6) - Lz;
                %                 Elements(j,i+6) = Elements(j,i+6) - xz;
                %                 Elements(j,i+1+6) = Elements(j,i+1+6) - yz;
            elseif (traj(j,i+2+6)-traj(j,i+2)) < -Lz/2;
                traj(j,i+2+6) = traj(j,i+2+6) + Lz;
                %                 Elements(j,i+6) = Elements(j,i+6) + xz;
                %                 Elements(j,i+1+6) = Elements(j,i+1+6) + yz;
            end
            
            if (traj(j,i+1+6)-traj(j,i+1)) > Ly/2;
                traj(j,i+1+6) = traj(j,i+1+6) - Ly;
                %                 Elements(j,i+6) = Elements(j,i+6) - xy;
            elseif (traj(j,i+1+6)-traj(j,i+1)) < -Ly/2;
                traj(j,i+1+6) = traj(j,i+1+6) + Ly;
                %                 Elements(j,i+6) = Elements(j,i+6) + xy;
            end
            
            if (traj(j,i+6)-traj(j,i)) > Lx/2;
                traj(j,i+6) = traj(j,i+6) - Lx;
            elseif (traj(j,i+6)-traj(j,i)) < -Lx/2;
                traj(j,i+6) = traj(j,i+6) + Lx;
            end
        else
            traj(j,i:i+8)=0;
        end
    end
end

dipole_vec = zeros(size(traj,1),size(traj,2)/3);

for k=1:3;
    dipole_vec(:,k:3:end) = traj(:,k:9:end) - (traj(:,k+3:9:end) + traj(:,k+6:9:end))/2;
end

dipole_vec=dipole_vec.*2.3506/0.58;

assignin('caller','Element_all',traj);


