%% G2_atom.m
% * This function calculates the continuous G2 factor fromthe cos and sin
% * terms and also saves a struct variable for G2_calc_func()
% * You might want to edit the atomtype names below to fit your needs...
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% * [G2,G_cos,G_sin,Gtwotheta,structure] = G2_atom(atom,Box_dim)
% * [G2,G_cos,G_sin,Gtwotheta,structure] = G2_atom(atom,Box_dim,filename)
%
function [G2,G_cos,G_sin,Gtwotheta,structure,old_structure] = G2_atom(atom,Box_dim,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format compact;

if nargin==2
    filename='outfile';
elseif nargin>2
    filename=varargin{1};
end
Inorganic_atoms = {'Si','Al','Alt','Feo','Mgo','Ob','Op','O','Oh','Omg','Ohmg','Oalt','Odsub','H','Ow','Hw','OW','HW','HW1','HW2','Li','Na','K','Cs','Mg','Ca','Sr','Ba','Cl','Br','Cc','Hc','Oc','Nc'}; % Octahedral and tetrahedral substitusions
Atoms_noWater  = {'Si','Al','Alt','Feo','Mgo','Mg','Ob','Op','O','Oh','Omg','Ohmg','Oalt','Odsub','H','Na','Cc','Hc','Oc','Nc'};%,'Li','Na','K','Cs','Mg','Ca','Sr','Ba','Cl'};
Atoms_clay     = {'Si','Al','Alt','Feo','Mgo','Mg','Ob','Op','O','Oh','Omg','Ohmg','Oalt','Odsub','H'};%,'Li','Na','K','Cs','Mg','Ca','Sr','Ba','Cl'};
System         = {'all'};%

atom=translate_atom(atom,-[atom.x atom.y atom.z]);
atom=wrap_atom(atom,Box_dim);

Elements(1,1:3:3*length([atom.x]))=[atom.x];
Elements(1,2:3:3*length([atom.x]))=[atom.y];
Elements(1,3:3:3*length([atom.x]))=[atom.z];

Elements(:,1:3:end)=Elements(:,1:3:end) - min(min(Elements(:,1:3:end)));
Elements(:,2:3:end)=Elements(:,2:3:end) - min(min(Elements(:,2:3:end)));
Elements(:,3:3:end)=Elements(:,3:3:end) - min(min(Elements(:,3:3:end)));

% Element = Wrap_Coord_func(Elements, Box_dim);

% ind_Feo=find(ismember([atom.type],'Al'));
% ind_Feo=ind_Feo(1:9:end);
% [atom(ind_Feo).type]=deal({'Feo'});

for push=1:1
    
    scale_water=1;
    nUC=6*4*3;
    dim='z';
    stride=1;
    pushSOL=0;%(-5+push)/10; %Å
    nLayers=3;
    lambda=1.54187;
    step=0.01; % Do not change
    twotheta=[2:step:60]';
    Gtwotheta=[2:.05:60]; % To be used in the modelling
    Linecolor = [0 0 1];% [1-MooreReynolds-CLAYFF 0 MooreReynolds];
    Linewidth = 2;
    Fontsize = 18;
    Add_Atom1 = {'Feo'}; Replace_Atom1={'Al'}; Atom1_frac=0;% 11in percent, normally < 4 Al per UC
    Add_Atom2 = {'Alt'}; Replace_Atom2={'Si'}; Atom2_frac=0;% 1in percent, normally < 8 Si per UC
    centerAtomtype='Al';
    center_ind=find(strcmp([atom.type],centerAtomtype))';
    %     center_ind=center_ind(length(center_ind)/nLayers+1:2*length(center_ind)/nLayers);%         center_ind=center_ind(1:length(center_ind)/nLayers);
    scattering_factors='default';%'WaasmaierKirfel';%'Wright'; % 'Wright' or 'CromerMann'
    
    for i=1:size(atom,2)
        temp_label=[atom(i).type];
        if size(temp_label{1},2)
            temp_label{1}(2:end)=lower(temp_label{1}(2:end));
            [atom(i).type]=temp_label;
        end
    end
    
    inorganic_ind=find(ismember([atom.type],Inorganic_atoms));
    Cc_ind=find(strncmpi([atom.type],'C',1));
    Cc_ind=setdiff(Cc_ind,inorganic_ind);
    [atom(Cc_ind).type]=deal({'Cc'})
    
    Hc_ind=find(strncmpi([atom.type],'H',1));
    Hc_ind=setdiff(Hc_ind,inorganic_ind);
    [atom(Hc_ind).type]=deal({'Hc'})
    
    Oc_ind=find(strncmpi([atom.type],'O',1));
    Oc_ind=setdiff(Oc_ind,inorganic_ind);
    [atom(Oc_ind).type]=deal({'Oc'})
    
    Nc_ind=find(strncmpi([atom.type],'N',1));
    Nc_ind=setdiff(Nc_ind,inorganic_ind);
    [atom(Nc_ind).type]=deal({'Nc'})
    
    if strcmpi(System,{'all'})
        Atom_label = unique([atom.type]);
    else
        Atom_label = System(ismember(System,unique([atom.type])))
    end
    
    Ow_ind=find(strncmpi([atom.type],'Ow',2));[atom(Ow_ind).type]=deal('Ow');
    Hw_ind=find(strncmpi([atom.type],'Hw',2));[atom(Hw_ind).type]=deal('Hw');
    Atom_label(find(strncmpi(Atom_label,{'MW'},2)))=[];
    Atom_label(find(strncmpi(Atom_label,{'HW2'},3)))=[];
    Atom_label(find(strncmpi(Atom_label,{'Ow'},2)))={'Ow'};
    Atom_label(find(strncmpi(Atom_label,{'Hw'},2)))={'Hw'};
    
    if strcmp(dim,'x')
        dim_column=2;
        L = Box_dim(1);
        d001=L/nLayers;
        Distance=0:step:L;
    elseif strcmp(dim,'y')
        dim_column=1;
        L = Box_dim(2);
        d001=L/nLayers;
        Distance=0:step:L;
    else%if strcmp(dim,'z');
        dim_column=0;
        L = Box_dim(3);
        d001=L/nLayers;
        Distance=0:step:L;
    end
    
    ind_center_atom=3*center_ind-dim_column;
    center_atom=reshape(Elements(:,ind_center_atom),[],1);
    center_atom_hist=hist(center_atom,Distance)';
    [Max_center,Max_row]=max(center_atom_hist);
    shift=Max_row*step;
    Distance_symm=Distance-d001/2+3/2*step;
    
    Total_density=zeros(ceil(d001/step),1);
    Total_density_disc=zeros(ceil(d001/step),1);
    G_cos=zeros(size(twotheta,1),1);
    G_sin=zeros(size(twotheta,1),1);
    
    %hold on;
    
    % for h=1:length(Atom_label)
    %     ind_atom=3*find(strcmp([atom.type],Atom_label(h)))'-dim_column;
    %     Element = reshape(Elements(:,ind_atom),[],1);
    %     Element = Element - shift + d001/2;
    %     Element = Wrap_Coord_func(Element, L);
    %     Element_density = histcounts(Element,Distance)';
    %     Element_density=(Element_density+flipud(Element_density))/2;
    %     Element_density=Element_density(1:ceil(d001/step),1);
    %     Total_density=Total_density+Element_density;
    % end
    
    SUMNDISCRETE=0;
    norm=size(Elements,1)*nUC;
    Total_density(:)=0;Total_Electron_density=0;
    Element_count=0;
    
    for h=1:length(Atom_label)
        
        Element_ave_sum=0;
        Element_count=Element_count+1;
        ind_atom=3*find(strcmpi([atom.type],Atom_label(h)))'-dim_column;
        disp('Atom_label')
        Atom_label(h)
        size(ind_atom)
        disp('-----------')
        Element = reshape(Elements(:,ind_atom),[],1);
        Element = Element - shift + d001/2;
        Element = Wrap_Coord_func(Element, L);
        Element_density = histcounts(Element,Distance)'/norm;
        Element_density = Element_density-floor(Element_density./d001)*d001;
        
        Element_density = Element_density(1:ceil(d001/step),1);
        Element_density = smooth((Element_density+flipud(Element_density))/2,11);
        
        DW=0;
        %% Set the Debye-Waller factor
        if ismember(Atom_label(h),{'Si','Al','Alt','Mgo','Mg','H','Hc','Fe','Feo'})
            DW=1.5
        elseif ismember(Atom_label(h),{'Li','Na','K','Cs','Mg','Ca','Sr','Ba'})
            DW=2
        elseif ismember(Atom_label(h),{'Op','Ob','O','Oh','Oapical','Obasal','Omg','Ohmg','Oalt','Odsub','Br','Cl'})
            DW=2
        elseif ismember(Atom_label(h),{'Ow','Hw'})
            disp('atomtype --> water')
            DW=2
        else
            disp('Do not recognise the atomtype!!!')
            DW=2
        end
        
        if ismember(Atom_label(h),{'Ow','Hw','OW','HW','HW1','HW2'})
            Element_density=Element_density*scale_water;
            if pushSOL>0
                extendedElement_density=[Element_density(1:length(Element_density)/2);zeros(1,floor(2*pushSOL/step))';Element_density(length(Element_density)/2+1:end)];
                Element_density=interp1(1:length(extendedElement_density),extendedElement_density,1:length(Element_density))';
                Element_density(1:length(Element_density)/2)=0;
                Element_density = sum(extendedElement_density)/sum(Element_density)*(Element_density+flipud(Element_density))/2;
            elseif pushSOL<0
                truncatedElement_density=Element_density(1:length(Element_density)/2-floor(-pushSOL/step));
                truncatedElement_density=[truncatedElement_density;flipud(truncatedElement_density)];
                Element_density=interp1(1:length(truncatedElement_density),truncatedElement_density,1:length(Element_density))';
                Element_density(isnan(Element_density))=truncatedElement_density(end);
                Element_density(1:length(Element_density)/2)=0;
                Element_density = sum(truncatedElement_density)/sum(Element_density)*(Element_density+flipud(Element_density))/2;
                if ismember(Atom_label(h),{'Ow' 'OW'})
                    disp('water content')
                    sum(Element_density)
                end
            end
        end
        
        Element_f = atomic_scattering_factors(Atom_label(h),lambda,twotheta,DW);
        Atomtype
        Element_Electron_density =  Element_density*nElectrons;
        
        Element_cos=zeros(size(twotheta,1),1); Element_sin=zeros(size(twotheta,1),1);
        for i=1:size(twotheta,1)
            PhaseAngle=2*pi()*Distance_symm(1:length(Element_density)).*(2*sin(twotheta(i)/2*pi()/180)/lambda);
            Element_cos(i,1) =  Element_f(i)*Element_density.'*cos(PhaseAngle)';
            Element_sin(i,1) =  Element_f(i)*Element_density.'*sin(PhaseAngle)';
        end
        %                 Element_cos=zeros(size(twotheta,1),1); Element_sin=zeros(size(twotheta,1),1);
        %                 for i=1:size(twotheta,1);
        %                     for j=1:ceil(d001/step);
        %                         cos_term(i,j) = Element_density(j) * Element_f(i) * cos(2*pi()*Distance(j)*(2*sin(twotheta(i)/2*pi()/180)/lambda));
        %                         sin_term(i,j) = Element_density(j) * Element_f(i) * sin(2*pi()*Distance(j)*(2*sin(twotheta(i)/2*pi()/180)/lambda));
        %                     end
        %                     Element_cos(i,1)=sum(cos_term(i,:),2);
        %                     Element_sin(i,1)=sum(sin_term(i,:),2);
        %                 end
        Element_G2 = Element_cos.^2+Element_sin.^2;
        Total_density=Total_density+Element_density;
        %         Total_density_disc=Total_density_disc+Element_density_disc;
        Total_Electron_density = Total_Electron_density+Element_Electron_density;
        G_cos=G_cos+Element_cos;
        G_sin=G_sin+Element_sin;
        
        if sum(Element_density) > 0
            SUMNDISCRETE=SUMNDISCRETE+Element_ave_sum;
            %                 assignin('base',strcat(char(Atom_label(h)),'_discr_sum'),Element_ave_sum); % new
            %                 assignin('base',strcat(char(Atom_label(h)),'_discr'),Element_ave); % new
            assignin('base',strcat(char(Atom_label(h)),'_cos'),Element_cos); % new
            assignin('base',strcat(char(Atom_label(h)),'_sin'),Element_sin); % new
            assignin('base',strcat(char(Atom_label(h)),'_G2'),Element_G2);
            assignin('base',strcat(char(Atom_label(h)),'_density'),Element_density);
            assignin('base',strcat(char(Atom_label(h)),'_electron_density'),Element_Electron_density);
            assignin('base',strcat(char(Atom_label(h)),'_f'),Element_f);
            %plot(twotheta,Element_f)
        end
        
        %             hold on;
        %             plot([UC(Element_count).Zdist],[UC(Element_count).Pdist])
        %             findpeaks([UC(Element_count).Pdist],[UC(Element_count).Zdist]','MinPeakDistance',0.1,'MinPeakProminence',.0005)
        %             plot(Distance_symm(1:length(Total_density)),circshift(Element_density,[0 0]))
        %             stem(Distance(1:length(Element_density)),circshift(Element_density,[0 0]),'Marker','none')
        %             plot(Distance(1:length(Element_density)),circshift(Element_density,[ceil(d001/2/step) 0]))
        %             stem(Distance(1:length(Element_density_disc)),circshift(Element_density_disc,[ceil(d001/2/step) 0]),'Marker','none')
        
        
        %% Collect the structural data
        structure(h).type=Atom_label(h);
        structure(h).d001=d001;
        structure(h).P=single(Element_density)';%UC(h).P;;%single(UC(h).P);
        structure(h).e=nElectrons;
        structure(h).Z=Distance_symm(1:length(Element_density));%(find(Element_density));;%Distance_symm(find(Element_density));
        structure(h).B=DW;
        
        
    end
    
    %% Clear up structure
    k=1;i=1;
    newstructure=structure;
    while i < size(structure,2)+1
        for j=1:numel([structure(i).P])
            if numel([structure(i).P]) > 1
                newstructure(k)=structure(i);
                newstructure(k).Z=structure(i).Z(j);
                newstructure(k).P=structure(i).P(j);%/numel([structure(i).Z]);
                newstructure(k).e=[structure(i).P(j)]*[structure(i).e];%/numel([structure(i).Z]);
            end
            k=k+1;
        end
        i=i+1;
    end
    
    ind_zero=find([newstructure.P]==0);
    newstructure(ind_zero)=[];
    old_structure=structure;
    structure=newstructure;
    
    G_cos=interp1(twotheta,G_cos,Gtwotheta);
    G_sin=interp1(twotheta,G_sin,Gtwotheta);
    
    G2=G_cos.^2+G_sin.^2;
    save(strcat('G2md_',filename,'.mat'),'G2','G_cos','G_sin','G2','Gtwotheta','structure');
    
    figure
    hold on
    plot(Gtwotheta,G2,'-','color',Linecolor,'LineWidth',Linewidth);set(gca,'PlotBoxAspectRatio',[1 1 1],'FontSize',Fontsize)
    plot(Gtwotheta,G_cos,'-','color',Linecolor,'LineWidth',Linewidth);set(gca,'PlotBoxAspectRatio',[1 1 1],'FontSize',Fontsize)
    plot(Gtwotheta,G_sin,'-','color',Linecolor,'LineWidth',Linewidth);set(gca,'PlotBoxAspectRatio',[1 1 1],'FontSize',Fontsize)
    
    figure
    hold on
    plot(Distance(1:length(Total_density)),circshift(Total_density,[0 0]))
%     plot(Distance_symm(1:length(Total_density)),circshift(Total_Electron_density,[0 0]))
    %     xlabel('2Theta','FontSize',Fontsize); ylabel('a.u.','FontSize',Fontsize);
end

