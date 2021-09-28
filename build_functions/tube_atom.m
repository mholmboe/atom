%% tube_atom.m
% * This quirky function can be used to create a nano-tube or nano-roll of the 
% coordinates from an atom struct. It works best if the the input atom struct 
% consists of one centered unit cell (to keep the number of atoms down). 
% You can always the replicate_atom() function later to build the entire roll/tube.
% * See needed variables/parameters below on line 22-27 and 37 to play around with.
% For non-centrosymmetric layers, chosing +R or -R on line 45 allow you to
% choose the type of inner-surface atoms. The spiral_vector and Rshift on 
% 37 can be used to skew the spiral in the x and/or the y direction/s.
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = tube_atom(atom,Box_dim,Radii)

function atom = tube_atom(atom,Box_dim,Dim,nUC,AngularRange,UCaccuracy,deltaaccuracy,Rshift)

% Dim=2; % 1 for X-direction, 2 for the Y-direction
% nUC=36; % Ideal number of unit cells per revolution
% AngularRange=2*360; % Ideal number of revolutions, ie n*360deg. The final angular range might be different..
% UCaccuracy=0.05; % Connected to UC placement accuracy, lower is better, but harder.. will affect the actual angular range
% deltaaccuracy=0.4; % Connected to UC placement speed, higher is faster, but worse.. will affect the actual angular range
% Rshift=0.3; % Lateral shift in Å between each UC, resulting in a roll rather than a tube (Rshift=0)

if Dim==1
    L=Box_dim(1);
elseif Dim==2
    L=Box_dim(2);
elseif Dim==3
    L=Box_dim(3);
end

spiral_vector = [0   0   Rshift*(Rshift/(Rshift+i/AngularRange/2))] % [0 0 0] makes a tube, anything else creates a spiral roll perpendicular to [x y z]
R = L*nUC/(2*pi) % Target inner radii
Cf = 2*pi*R % Target inner circumference, not used

atom=translate_atom(atom,-[min([atom.x]) min([atom.y]) max([atom.z])]);
atom=replicate_atom(atom,Box_dim,[1 1 1]);

atom=translate_atom(atom,-Box_dim./2);
atom=translate_atom(atom,[0 0 R]); % Try also -R, do you see any difference

System=atom;
nUCreal=nUC*AngularRange/360;
prev_temp=atom;AngleList(1)=0;DeltaList(1)=0;

for i=1:nUCreal
    A=0;delta=0;angle=0;MajorAtomicOverlap=0;
    while abs(A-L)>(UCaccuracy*L)
        
        if A-L > UCaccuracy*L
            angle=(i-1)*AngularRange/nUCreal-delta;
        elseif A-L < UCaccuracy*L
            angle=(i-1)*AngularRange/nUCreal+delta;
        end
        
        if angle<AngleList(end)
            angle=AngleList(end);
        end
        
        %         if angle>AngularRange
        %             disp('angle larger than the AngularRange')
        %             [angle AngularRange]
        %
        %         end
        
        if Dim==1
            trans_temp = spiral_atom(atom,Box_dim,[angle 0 0],i*spiral_vector);
        elseif Dim==2
            trans_temp = spiral_atom(atom,Box_dim,[0 angle 0],i*spiral_vector);
        elseif Dim==3
            disp('Not supported in Z')
            pause
            %             trans_temp = spiral_atom(atom,Box_dim,[0 0 angle],i*[0 0 Rshift]);
        end
        
        D=dist_matrix_noPBC_atom(trans_temp,prev_temp);
        A=mean(diag(D));
        delta=delta+deltaaccuracy;
        
    end
    
    angle
    
    DeltaList(i)=delta;
    AngleList(i)=angle;
    System=update_atom({System trans_temp});
    prev_temp=trans_temp;
    
    %     if mod(i,10)==1
    %         if i > 1
    %             i-1
    %         end
    %     end
    
end

if Rshift==0
    D=dist_matrix_noPBC_atom(atom,trans_temp);
    A=mean(diag(D));
    disp('last-to-first distance         ideal distance')
    [A L]
    if A-L > UCaccuracy*L
        disp('Play around with AngularRange and UCaccuracy to optimize the sealing of the tube')
    elseif A-L < UCaccuracy*L
        disp('Play around with AngularRange and UCaccuracy to optimize the sealing of the tube')
    end
end

atom=System;

assignin('caller','DeltaList',DeltaList);
assignin('caller','AngleList',AngleList);
disp('.')
disp('.')
disp('.')
disp('Note that this function offers no error-checking in terms of overlapping atoms')
disp('..no clever math can be found here..')

