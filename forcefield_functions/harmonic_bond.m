function [r,Ubond] = harmonic_bond(dist,kb) % ff and atomtype can be single variables or cell 'tuplets'

% if ischar(atomtype)
%     disp('atomtype must be a {1x1} or {1x2} cell')
% end
% 
% if size(atomtype,2)==1
%     atomtype1=atomtype;
%     atomtype2=atomtype;
% elseif size(atomtype,2)==2
%     atomtype1=atomtype(1);
%     atomtype2=atomtype(2);
% end
% 
% mass1=mass_atom(Atom_label1);
% mass2=mass_atom(Atom_label2);
% 
% reduced_mass=(mass1*mass2)/(mass1+mass2);
r=dist(1)-.5:.0005:dist(1)+.5;
Ubond=kb*(r-dist(1)).^2;

hold on
plotmin=1000*(round(min(Ubond*1.5)/1000));
if plotmin>=0
    plotmin=10000;
end

if numel(dist)>1
    plot(r,Ubond+dist(2));
else
    plot(r,Ubond);
end
xlabel('r [nm]');
ylabel('U [kJ/mol]');
xlim([-.2,1.2]);
ylim(sort([-plotmin plotmin]))
% 
% titlestring=strcat(char(atomtype1.type),'-',char(atomtype2.type));
% titlestring=strrep(titlestring,'_','-');
% title(titlestring)
