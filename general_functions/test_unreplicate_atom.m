%% test_unreplicate_atom.m
% * This function checks that unreplicate_atom folds a supercell back onto
% * its unit cell, for orthogonal as well as triclinic boxes.
% * A supercell replicated by condense_factors and then condensed by the same
% * factors must put every replica of a site on top of that site, so the
% * folded structure holds exactly as many distinct positions as the unit
% * cell it came from. The box the function reports back must also match the
% * unit cell box it started from.
% * The triclinic cases are the ones worth keeping: folding the input
% * coordinates rather than the orthogonalized ones leaves the box tilt in
% * the in-plane positions, which splits each site into one position per
% * layer. On a 2x2x2 pyrophyllite cell that misplaced atoms by up to 3.46 A,
% * exactly the xz tilt, while every orthogonal box still passed.
%
%% Version
% 3.0.2
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # test_unreplicate_atom
% # ok = test_unreplicate_atom

%% See also
% # <unreplicate_atom.html unreplicate_atom>
% # <replicate_atom.html replicate_atom>
% # <orto_atom.html orto_atom>
% # <triclinic_atom.html triclinic_atom>
% # <wrap_atom.html wrap_atom>
%
function ok = test_unreplicate_atom

tol=0.01; % A, well above the 4 decimals the transforms round to

% name, Cell as [a b c alfa beta gamma]
cases={ 'orthogonal',            [5.1600 8.9658 9.3450 90.00 90.00 90.00];...
        'monoclinic pyrophyllite',[5.1600 8.9658 9.3450 90.00 100.46 90.00];...
        'triclinic sheared',     [5.1600 8.9658 9.3450 95.00 100.00 80.00]};

factors={[2 2 2],[3 2 1]};

% fractional sites, kept clear of the cell faces and of each other
frac=[0.1234 0.3456 0.5678;...
      0.7777 0.1111 0.2500;...
      0.4321 0.6789 0.8765;...
      0.2500 0.8888 0.3333;...
      0.6111 0.5000 0.7000;...
      0.9000 0.2222 0.1500];

ok=true;
fprintf('\n%-26s %-8s %7s %9s %8s %s\n','case','factors','sites','expected','spread','box');
for c=1:size(cases,1)
    name=cases{c,1};
    uc_box=Cell2Box_dim(cases{c,2});
    uc=sites_atom(frac,uc_box);
    n=numel(uc);

    for f=1:numel(factors)
        cf=factors{f};
        [sc,sc_box]=replicate_atom(uc,uc_box,cf);
        rec=unreplicate_atom(sc,sc_box,cf);
        rec_box=Box_dim; % unreplicate_atom assigns this into its caller

        R=[[rec.x]' [rec.y]' [rec.z]'];

        % how many distinct positions survive the fold, and how far apart the
        % replicas of a site ended up; both independent of the atom order
        [nsites,spread]=count_sites(R,tol);

        boxok=numel(rec_box)>=3 && max(abs(rec_box(1:3)-uc_box(1:3)))<tol;
        pass=(nsites==n) && (spread<tol) && boxok;
        ok=ok && pass;

        fprintf('%-26s %-8s %7d %9d %8.4f %s   %s\n',name,mat2str(cf),nsites,n,spread,...
            ternary(boxok,'ok','WRONG'),ternary(pass,'pass','FAIL'));
    end
end

fprintf('\n');
if ~ok
    error('test_unreplicate_atom:failed','unreplicate_atom did not fold a supercell back onto its unit cell');
end
disp('test_unreplicate_atom passed');

%% Group positions that sit within tol of each other, and report the count
% and the widest spread inside a group
function [k,spread] = count_sites(R,tol)
k=0; spread=0; left=R;
while ~isempty(left)
    k=k+1;
    d=sqrt(sum((left-left(1,:)).^2,2));
    sel=d<=tol;
    g=left(sel,:);
    spread=max(spread,max(max(g,[],1)-min(g,[],1)));
    left=left(~sel,:);
end

%% Build an atom struct at the given fractional sites of a Box_dim cell
function atom = sites_atom(frac,Box_dim)
lx=Box_dim(1); ly=Box_dim(2); lz=Box_dim(3);
if numel(Box_dim)>3
    xy=Box_dim(6); xz=Box_dim(8); yz=Box_dim(9);
else
    xy=0; xz=0; yz=0;
end
H=[lx 0 0; xy ly 0; xz yz lz]; % rows are the a, b and c box vectors
XYZ=frac*H;
for i=1:size(XYZ,1)
    atom(i).molid=1;
    atom(i).resname={'TST'};
    atom(i).type={'O'};
    atom(i).fftype={'O'};
    atom(i).index=i;
    atom(i).x=XYZ(i,1);
    atom(i).y=XYZ(i,2);
    atom(i).z=XYZ(i,3);
end

function s = ternary(cond,a,b)
if cond, s=a; else, s=b; end
