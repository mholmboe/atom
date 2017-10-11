%% hist_atom.m
% * This function is used to calculate density profiles in the Z-direction
% * Tested when? Does it work?
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * Hist = hist_atom(atom,0.02)

function Hist = hist_atom(atom,s)

ind_H=find(strncmp([atom.fftype],'H',1));
ind_O=find(strncmp([atom.fftype],'O',1));
ind_rm = [ind_H ind_O];
atom(ind_H)=[];

% ind_S=find(strncmp([atom.fftype],'S',1));
% atom=[atom atom(ind_S) atom(ind_S) atom(ind_S) atom(ind_S) atom(ind_S)];

ind_C=find(strncmp([atom.fftype],'C',1));
atom=[atom atom(ind_C) atom(ind_C)];

Coord = [atom.x; atom.y; atom.z]';

% xmin=min(Coord(:,1));
% ymin=min(Coord(:,2));
% zmin=min(Coord(:,3));

Coord = [Coord(:,1)-mean(Coord(:,1))+50 Coord(:,2)-mean(Coord(:,2))+50 Coord(:,3)-mean(Coord(:,3))+50];

Bins = (0:s:s*1000)';
Hist=hist(Coord,Bins);

sigma = 4;
window = 100;
x = linspace(-window / 2, window / 2, window);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum(gaussFilter); % normalize
%Hist(:,1)=filter(gaussFilter,1, Hist(:,1));

for i=1:3
    Hist(:,i)=conv(Hist(:,i), gaussFilter, 'same');
    Hist(:,i)=Hist(:,i)/sum(Hist(:,i));
end

% assignin('caller','Bins',Bins);
% assignin('caller','Coord',Coord);
% assignin('caller','xmin',xmin);
% assignin('caller','ymin',ymin);
% assignin('caller','zmin',zmin);


