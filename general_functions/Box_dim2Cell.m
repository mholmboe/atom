%% Box_dim2Cell.m
% * This function transforms the 1x3 or the 1x9 Box_dim variable to the 1x6
% Cell variable, containing the a, b, c cell values and  the alfa, beta,
% gamma angle values as used in a typical .pdb file
%

%% Version
% 2.08
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # Cell=Box_dim2Cell(Box_dim)

function Cell=Box_dim2Cell(Box_dim)

if length(Box_dim)==9
    Box_dim(Box_dim<0.00001&Box_dim>-0.00001)=0;
    if sum(find(Box_dim(4:end)))<0.0001
        Box_dim=Box_dim(1:3);
    end
end

if length(Box_dim)==3
    
    lx=Box_dim(1);
    ly=Box_dim(2);
    lz=Box_dim(3);
    xy=0;
    xz=0;
    yz=0;
    
    a=lx;
    b=ly;
    c=lz;
    alfa=90.00;
    beta=90.00;
    gamma=90.00;
    
    Cell=[a b c alfa beta gamma];
    
elseif length(Box_dim)==6
    
    a=Box_dim(1);
    b=Box_dim(2);
    c=Box_dim(3);
    alfa=Box_dim(4);
    beta=Box_dim(5);
    gamma=Box_dim(6);
    
    Cell=[a b c alfa beta gamma];
    
elseif length(Box_dim)==9
    
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
    
    Cell=[a b c alfa beta gamma];
    
else
    Cell=[];
    disp('No proper box_dim information')
end



end