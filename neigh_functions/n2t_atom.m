%% n2t_atom.m — write GROMACS atomname2type.n2t files
%
% Columns: Element  AtomType  Charge  Mass[u]  CN  [NeighborElement  Distance(nm)] x CN
% - Element <- element_atom(atom) -> element(i).type
% - AtomType <- original atom(i).type
% - CN = unique neighbor atoms
% - Neighbors = exactly CN entries (no aggregation)
% - Distances = mean per neighbor position (Å→nm)
% - Charge = mean per environment if atom(i).charge exists, else 0.0
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% txt = n2t_atom(atom,Box_dim,'atomtype2name.n2t'); % Basic input arguments
%

function txt = n2t_atom(atom, Box_dim, outfile, varargin)

if nargin < 3 || isempty(outfile), outfile = 'atomname2type.n2t'; end

% % % Optional explicit element struct
% % if ~isempty(varargin) && isstruct(varargin{1})
% %     element = varargin{1}; varargin(1) = [];
% % else
% %     element = [];
% % end

% Options
useCharges = true; autoBond = true; rmaxlong = 2.45; verbose = true; to_nm = 0.1;
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k})); if k+1>numel(varargin), break; end
    val = varargin{k+1};
    switch key
        case "usecharges",  useCharges = logical(val);
        case "autobond",    autoBond   = logical(val);
        case "rmaxlong",    rmaxlong   = val;
        case "verbose",     verbose    = logical(val);
    end
end

% Ensure masses
if ~isfield(atom,'mass') || isempty([atom.mass]), atom = mass_atom(atom); end

% Elements 
 element = element_atom(atom);

% if isempty(element), element = element_atom(atom); end
% if ~isfield(element,'type') || numel(element) ~= numel(atom)
%     error('element_atom(atom) must return element(i).type for each atom.');
% end
getEl = @(i) get_field_cell(element,i,'type','X');
getTy = @(i) get_field_str(atom,i,'type','X');

% % % Ensure neighbors/bonds
% % need_bonds = (~isfield(atom,'bond')) || all(arrayfun(@(a) isempty(a.bond), atom));
% % need_neigh = (~isfield(atom,'neigh')) || all(arrayfun(@(a) isempty(a.neigh), atom));
% % if autoBond && (need_bonds || need_neigh)
% %     if nargin<2 || isempty(Box_dim)
% %         warning('No bonds/neighbors and Box_dim missing; proceeding without distances.');
% %     else
         atom = bond_atom(atom, Box_dim, rmaxlong);
% %     end
% % end

% Build per-site neighbor sequences
% ----------------- per-site neighbor sequences -----------------
Natom = numel(atom);
site = struct('centerEl',{},'centerTy',{},'CN',{}, ...
              'neighElSeq',{},'neighDistsA',{},'hasCharge',{});

for i = 1:Natom
    centerEl = string(getEl(i)); if strlength(centerEl)==0, centerEl="X"; end
    centerTy = string(getTy(i)); if strlength(centerTy)==0, centerTy=centerEl; end

    % --- neighbors by index (unique atoms) ---
    neigh_idx  = get_field_vec(atom,i,'neigh.index',[]);
    if ~isempty(neigh_idx)
        neigh_idx = unique(neigh_idx,'stable');
    end

    % --- build distance map from bonds (robust), fallback to neigh.dist, then coords ---
    distMap = containers.Map('KeyType','double','ValueType','double');
    if isfield(atom,'bond') && ~isempty(atom(i).bond) && isfield(atom(i).bond,'index') ...
                              && ~isempty(atom(i).bond.index)
        bidx = atom(i).bond.index;   % [i j]
        bdis = atom(i).bond.dist;    % Å
        for b = 1:size(bidx,1)
            j = bidx(b,2);
            if j>=1 && j<=Natom
                distMap(double(j)) = double(bdis(b));
            end
        end
    end
    neigh_dist = get_field_vec(atom,i,'neigh.dist',[]); % may be empty/misaligned

    % --- assemble canonical neighbor element + distance arrays (Å) ---
    nEl = strings(0,1); nD = [];
    for k = 1:numel(neigh_idx)
        j = neigh_idx(k);
        if j>=1 && j<=Natom
            % element of neighbor j
            nej_el = string(getEl(j)); if strlength(nej_el)==0, nej_el="X"; end
            nEl(end+1,1) = nej_el; %#ok<AGROW>

            % distance preference: bond map -> neigh.dist(k) -> from coordinates
            if isKey(distMap,double(j))
                d = distMap(double(j));
            elseif ~isempty(neigh_dist) && k<=numel(neigh_dist) && isfinite(neigh_dist(k))
                d = double(neigh_dist(k));
            elseif all(isfield(atom, {'x','y','z'}))
                d = sqrt( (atom(i).x-atom(j).x)^2 + (atom(i).y-atom(j).y)^2 + (atom(i).z-atom(j).z)^2 );
            else
                d = NaN;
            end
            nD(end+1,1) = d; %#ok<AGROW>
        end
    end

    % --- sort neighbors by (element, then distance) for canonical ordering ---
    if ~isempty(nEl)
        % numeric secondary key via formatted strings for stable sort
        [~, ord] = sortrows([nEl, string(num2str_vector(nD))]);
        nEl = nEl(ord); nD = nD(ord);
    end

    site(end+1) = struct( ...
        'centerEl',    centerEl, ...
        'centerTy',    centerTy, ...
        'CN',          numel(nEl), ...
        'neighElSeq',  nEl, ...
        'neighDistsA', nD, ...
        'hasCharge',   isfield(atom,'charge') && ~isempty(atom(i).charge) ); %#ok<AGROW>
end
% ----------------- group identical environments -----------------

% Group identical environments: key = El|Ty|CN|ElSeq
envMap = containers.Map();
for s = 1:numel(site)
    el  = site(s).centerEl; ty = site(s).centerTy; cn = site(s).CN;
    seq = strjoin(site(s).neighElSeq,'|');
    key = sprintf('%s|%s|%d|%s', el, ty, cn, seq);
    if ~isKey(envMap,key)
        pos = poscell(site(s).neighDistsA); posDists = cell(1,cn);
        for p = 1:cn
            posDists{p} = []; if p<=numel(pos) && isfinite(pos{p}), posDists{p}(end+1,1) = pos{p}; end
        end
        envMap(key) = struct('centerEl',el,'centerTy',ty,'CN',cn,'seqEl',site(s).neighElSeq, ...
                             'posDists',{posDists},'members',1,'memberIdx',s,'charges',get_charge(atom,s));
    else
        E = envMap(key);
        pos = poscell(site(s).neighDistsA);
        for p = 1:max(numel(E.posDists), numel(pos))
            if p>numel(E.posDists) || isempty(E.posDists{p}), E.posDists{p} = []; end
            if p<=numel(pos) && isfinite(pos{p}), E.posDists{p}(end+1,1) = pos{p}; end
        end
        E.members   = E.members + 1;
        E.memberIdx = [E.memberIdx s]; %#ok<AGROW>
        qc = get_charge(atom,s); if ~isempty(qc), E.charges(end+1,1) = qc; end
        envMap(key) = E;
    end
end

% Render lines
keysEnv = envMap.keys;
records = struct('Element',{},'AtomType',{},'Charge',{},'Mass',{},'CN',{});
neighborBlocks = {};
for kk = 1:numel(keysEnv)
    E = envMap(keysEnv{kk});
    % Charge
    q = 0.0; if useCharges && ~isempty(E.charges), q = mean(E.charges); end
    % MASS: use masses of THIS atomtype in THIS environment (members), not element-average
    masses = [atom(E.memberIdx).mass];
    if isempty(masses)
        m = NaN;
    else
        m = masses(1);
        if any(abs(masses - m) > 1e-6)
            m = median(masses);
            if verbose
                warning('n2t_atom:massInconsistency', ...
                    'Multiple masses for AtomType %s (%s); using median %.5f', ...
                    E.centerTy, E.centerEl, m);
            end
        end
    end
    % Distances per position
    CN  = E.CN; seq = E.seqEl; dnm = zeros(CN,1);
    for p = 1:CN
        if p<=numel(E.posDists) && ~isempty(E.posDists{p}), dnm(p) = mean(E.posDists{p})*to_nm; else, dnm(p)=0.0; end
    end
    rec.Element = char(E.centerEl); rec.AtomType = char(E.centerTy);
    rec.Charge = q; rec.Mass = m; rec.CN = CN; records(end+1) = rec; %#ok<AGROW>
    block = cell(CN,2);
    for p=1:CN, block{p,1}=char(seq(p)); block{p,2}=dnm(p); end
    neighborBlocks{end+1} = block; %#ok<AGROW>
end

% Order: descending CN, then Element, then AtomType
if ~isempty(records)
    CNlist = arrayfun(@(r) r.CN, records);
    elNames = string({records.Element}).'; tyNames = string({records.AtomType}).';
    [~, order] = sortrows([ -CNlist(:), elNames, tyNames ]);
else
    order = [];
end

header = [
    "; atomname2type.n2t generated by n2t_atom()"
    "; Columns: Element  AtomType  Charge  Mass[u]  CN  [NeighborElement  Distance(nm)]"
    "; Element from element_atom(atom) -> element(i).type ; AtomType from original atom(i).type"
    "; Mass from atom(i).mass for the members of each (Element,AtomType,CN,NeighborSeq) environment."
];

lines = strings(0,1);
for ix = order'
    r = records(ix);
    line = sprintf('%-2s %-16s % .6f %10.5f %2d', r.Element, r.AtomType, r.Charge, r.Mass, r.CN);
    block = neighborBlocks{ix};
    for k = 1:r.CN
        line = sprintf('%s  %-2s %6.3f', line, block{k,1}, block{k,2});
    end
    lines(end+1,1) = string(line); %#ok<AGROW>
end

txt = strjoin([header; lines], newline);
fid = fopen(outfile,'w');
if fid<0
    warning('Could not open %s for writing. Returning text.', outfile);
else
    fprintf(fid,'%s\n',txt); fclose(fid);
    if verbose, fprintf('Wrote %d entries to %s\n', numel(lines), outfile); end
end
end

% ---------------- helpers ----------------
function val = get_field_vec(S,i,field,default)
try
    parts = strsplit(field,'.'); v = S(i).(parts{1});
    for t = 2:numel(parts), v = v.(parts{t}); end
    if isempty(v), val = default; else, val = v; end
catch, val = default; end
end

function val = get_field_cell(S,i,field,default)
try, val = S(i).(field); if isempty(val), val = default; end
catch, val = default; end
end

function val = get_field_str(S,i,field,default)
try, val = S(i).(field); if isempty(val), val = default; end
    if isstring(val), val = char(val); end
catch, val = default; end
end

function out = poscell(v)
if isempty(v), out = {}; else, out = arrayfun(@(x) x, v, 'UniformOutput', false); end
end

function s = num2str_vector(v)
if isempty(v), s = ""; else, s = string(arrayfun(@(x) sprintf('%12.6f',x), v, 'UniformOutput', false)); end
end

function q = get_charge(atom, idx_site)
i = idx_site;
if isfield(atom,'charge') && ~isempty(atom(i).charge), q = atom(i).charge; else, q = []; end
end
