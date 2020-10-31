%% Cell2Box_dim.m
% * This function transforms the 1x6 Cell variable containing the a, b, c 
% cell values and  the alfa, beta, gamma angle values as used in a typical 
% .pdb file, into a 1x3 or the 1x9  Box_dim variable  
%

%% Version
% 2.08
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # Box_dim = Cell2Box_dim(Cell)

function Box_dim = Cell2Box_dim(Cell)
    a=Cell(1);
    b=Cell(2);
    c=Cell(3);
    alfa=Cell(4);
    beta=Cell(5);
    gamma=Cell(6);
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
end