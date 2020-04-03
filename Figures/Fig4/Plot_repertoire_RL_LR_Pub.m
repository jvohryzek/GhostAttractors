% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)
% %
% % Function to plot Figure 4
% % Plots the repertoire of Centroids obrained for the LR session
% % Plots the probabilities, Lifetimes for the 2 sessions.
%
% % It creates a new Figure with all K states, their probabilities,
% % lifetimes and list of brain areas.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Adapted version
% Jakub Vohryzek and Joana Cabral
% jakub.vohryzek@queens.ox.ac.uk and joanacabral@med.uminho.pt
% Ghost Attractors in Spontaneous Brain Activity: Recurrent Excursions Into
% Functionally-Relevant BOLD Phase-Locking States. (Vohryzek et al. 2020)
% doi: 10.3389/fnsys.2020.00020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DIRECTORY
statTest = 'correlation'% 'ttest'

Directory='/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/';

addpath(genpath(Directory))
Extension = '_unfiltered_100unrelated'; % '_filtered_100unrelated'

load([Directory 'LEiDA_HCP/K5_results_LR_RL' Extension],'Centroids_K5_LR',...
    'P_K5_LR','P_K5_RL','LT_K5_LR','LT_K5_RL')
%%
K=5;
Parcellation= 'AAL116';
N_areas=90;
n_Subjects=99; % 95

[cc_V_yeo7,p_V_yeo7,rangeK] = Overlap_LEiDA_Yeo (Parcellation,N_areas,Extension,0);
cc_k5_V_yeo7=squeeze(cc_V_yeo7(rangeK==K,1:K,:));
p_k5_V_yeo7=squeeze(p_V_yeo7(rangeK==K,1:K,:));

color_vecs = [120 18 134; 70 30 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78]./256;

Yeo_names={'VIS','SMT','DAT','VAT','LBC','FPA','DMN'};
%Yeo_names={'Visual','Somatomotor','Dorsal Att.','Ventral Att.','Limbic','Frontoparietal','Default Mode'};

Order=[1:2:N_areas N_areas:-2:2];

% LEiDA networks colored according to closest RSN
Volume=struct2array(load('ParcelsMNI2mm',['V_' Parcellation])); 
Brain_Mask_tmp=load_nifti('MNI152_T1_2mm_brain_mask.nii');
Brain_Mask=Brain_Mask_tmp.vol;

scortex=smooth3(Brain_Mask>0);
%%
figure('Name',['LEiDA 90 K=5 Centroids V1 ' Extension])
colormap(jet)
% Pannel A - Plot the BOLD PL patterns on the cortex (top and side)
% Pannel B - Plot the FC patterns in matrix format
% Pannel C - Plot the probability of each state in each condition

for c=1:K
    
    disp(['c=' num2str(c)])
    
    %% First Plot the arrows view from top
    subplot_tight(6,2*K,c*2-1,0.01)
    plot_arrows_in_cortex_general(Parcellation,Centroids_K5_LR(c,:))
    material dull; lighting gouraud
    daspect([1 1 1])
    view(-90,90)    % Top view    Side view:   view(0,0)     
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    %% Then render the areas on a transparent cortex view from top
    subplot_tight(6,2*K,c*2,0.01) 
    hold on
    cortexpatch=patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    
    % Find Yeo network with maximal overlap
    [~, net]= max(cc_k5_V_yeo7(c,:)); 
    Centroid_Vol=zeros(size(Volume)); 
    % Color all areas with positive V
    for n=find(Centroids_K5_LR(c,:)>0)
        Centroid_Vol(Volume==n)=1;
    end
    sregion=smooth3(Centroid_Vol>0);
    if p_k5_V_yeo7(c,net)<0.05/K
            patch(isosurface(sregion,0.3), 'FaceColor', color_vecs(net,:), 'EdgeColor', 'none') %'FaceAlpha', V(n))
    end
    %title({['Centroid #' num2str(c)]},'FontSize',6)
    %set(get(gca,'title'),'Position',[50.8314 98.5068 39.3312])
    material dull; lighting gouraud
    daspect([1 1 1])
    view(-90,90)    % Top view    Side view:   view(0,0)     
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    %% Same but view from the side
    % First Plot the arrows view from top
    subplot_tight(6,2*K,c*2-1+2*K,0.01)
    plot_arrows_in_cortex_general(Parcellation,Centroids_K5_LR(c,:))
    material dull; lighting gouraud
    daspect([1 1 1])
    view(0,0)    % Top view    Side view:   view(0,0)     
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    %% Then render the areas on a transparent cortex view from top
    subplot_tight(6,2*K,c*2+2*K,0.01) 
    hold on
    cortexpatch=patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    if p_k5_V_yeo7(c,net)<0.05/K
            patch(isosurface(sregion,0.3), 'FaceColor', color_vecs(net,:), 'EdgeColor', 'none') %'FaceAlpha', V(n))
    end
    material dull; lighting gouraud
    daspect([1 1 1])
    view(0,0)    % Top view    Side view:   view(0,0)     
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    %% Plot the FC matrices (most people studying FC like to see matrices
    
    subplot_tight(6,K,2*K+c,0.05)
    V=Centroids_K5_LR(c,:);
    V=V/max(abs(V));
    FC_V=V(Order)'*V(Order);
    %li=max(abs(FC_V(:)));
    imagesc(FC_V,[-1 1])
    axis square
    title(['FC state ' num2str(c)])
    if c==1
        ylabel('Brain area #')
        set(gca,'YTick',20:20:80)

    end
    xlabel('Brain area #')
    %set(gca,'XTick',20:20:80,'Fontsize',12)

    
    subplot_tight(6,K,3*K+c,0.06)
    switch statTest
        case 'ttest'
    hold on
    % P=[squeeze(P_K5_LR(:,c)) squeeze(P_K5_RL(:,c))];
    h1 = violinplot([P_K5_LR(:,c) P_K5_RL(:,c)]);
    h1(1,1).ViolinColor = [1 .7 .1];h1(1,2).ViolinColor = [.8 .3 .1];
    h1(1,1).BoxColor = [0.5 0.5 0.5];h1(1,1).BoxWidth = 0.005;
    h1(1,2).BoxColor = [0.5 0.5 0.5];h1(1,2).BoxWidth = 0.005;
    h1(1,1).MedianColor = [1 1 1];h1(1,2).MedianColor = [1 1 1];
    xticks(1:2);xticklabels({'LR','RL'});ylim([0 1]);
    hold on
    [hP(c) pP(c)]=ttest(P_K5_LR(:,c),P_K5_RL(:,c))
    if pP(c)<(0.05/K)
     % To plot an asterisk when significant
        plot(1.5,max([max(P_K5_LR(:,c)),max(P_K5_LR(:,c))]),'*k')
    end
    %histogram(P_K5_LR(:,c),0:0.05:1,'FaceColor',[1 .7 .1]);
    %histogram(P_K5_RL(:,c),0:0.05:1,'FaceColor',[.8 .3 .1]);
    
%     [t, P_val_2sessions]=ttest2(squeeze(P_K5_LR(:,c)),squeeze(P_K5_RL(:,c)),0.05/K)
% %     title(['p=' num2str(P_val_2sessions*K)])
% %     disp(['P p=' num2str(P_val_2sessions*K)])
    
    %plot([mean(P(:,rangeK==K,c)) mean(P(:,rangeK==K,c))],[0 100],':k')
    xlabel('Fractional Occupancy')
    if c==1
        ylabel('Probability','Fontsize',16)
    end
    %xlim([0 1])
    %ylim([0 40])
    box off
    if c==5
        %legend({'LR fMRI session','RL fMRI session'})
        lP = legend([h1(1).ScatterPlot;h1(2).ScatterPlot],'LR fMRI session','RL fMRI session')
        set(lP,'Location','northeast', 'FontSize',10)
    end
    
    %%
        case 'correlation'
            hold on
            if c==1
                plot(squeeze(P_K5_LR(:,c)),squeeze(P_K5_RL(:,c)),'.','Color',[0.5 0.5 0.5],'MarkerSize',8)
            else
                plot(squeeze(P_K5_LR(:,c)),squeeze(P_K5_RL(:,c)),'.','Color',color_vecs(net,:),'MarkerSize',8)
            end
            hold on
            plot(0:0.1:1,0:0.1:1,'r','LineWidth',2)
            ylim([0 1]);xlim([0 1])
    if c==1
        ylabel('RL: Fractional Occupancy');xlabel('LR: Fractional Occupancy')
    end
    axis square
    box off
    end
    
    subplot_tight(6,K,4*K+c,0.06)
    %%
    switch statTest
        case 'ttest' 
    hold on
    max_Lt=max(max([LT_K5_LR(:) LT_K5_RL(:)]));
    h1 = violinplot([LT_K5_LR(:,c) LT_K5_RL(:,c)]);
    h1(1,1).ViolinColor = [.1 .7 1];h1(1,2).ViolinColor = [.1 .3 .8];
    h1(1,1).BoxColor = [0.5 0.5 0.5];h1(1,1).BoxWidth = 0.005;
    h1(1,2).BoxColor = [0.5 0.5 0.5];h1(1,2).BoxWidth = 0.005;
    h1(1,1).MedianColor = [1 1 1];h1(1,2).MedianColor = [1 1 1];
    xticks(1:2);xticklabels({'LR','RL'});ylim([0 10])
    %histogram(LT_K5_LR(:,c),0:round(max_Lt)/20:max_Lt,'FaceColor',[.1 .7 1]);
    %histogram(LT_K5_RL(:,c),0:round(max_Lt)/20:max_Lt,'FaceColor',[.1 .3 .8]);
    
    hold on
    [hLT(c) pLT(c)]=ttest(LT_K5_LR(:,c),LT_K5_RL(:,c))
    if pLT(c)<(0.05/K)
     % To plot an asterisk when significant
        plot(1.5,max([max(LT_K5_LR(:,c)),max(LT_K5_LR(:,c))]),'*k')
    end
    
    if c==5
        lLT = legend([h1(1).ScatterPlot;h1(2).ScatterPlot],'LR fMRI session','RL fMRI session')
        %lLT.position = 
        set(lLT,'Location','northeast', 'FontSize',10)
    end
    %    plot([mean(LT(:,rangeK==K,c)) mean(LT(:,rangeK==K,c))],[0 100],':k')
    xlabel('Dwell Time (secs)')
    if c==1
        ylabel('# of TRs','Fontsize',16)
    end
    
%     [t, P_val_2sessions]=ttest2(squeeze(LT_K5_LR(:,c)),squeeze(LT_K5_RL(:,c)),0.05/K)
%     title(['p=' num2str(P_val_2sessions*K)])
%     disp(['LT p=' num2str(P_val_2sessions*K)])
    
    %xlim([0 max_Lt])
    %ylim([0 45])
    box off
    
        case 'correlation'
           hold on
           if c ==1
            plot(squeeze(LT_K5_LR(:,c)),squeeze(LT_K5_RL(:,c)),'.','Color',[0.5 0.5 0.5],'MarkerSize',8)
           else
            plot(squeeze(LT_K5_LR(:,c)),squeeze(LT_K5_RL(:,c)),'.','Color',color_vecs(net,:),'MarkerSize',8)
           end
            hold on
            plot(0:1:10,0:1:10,'r','LineWidth',2)
            ylim([0 10]);xlim([0 10])
    if c==1
        ylabel('RL: Dwell Time');xlabel('LR: Dwell Time')
    end
    axis square
    box off   
            
    end
   
    if c ==1
        display('state 1 no correlation')
    else
    subplot_tight(6,K,5*K+c,0.06)

    hold on
    for net=1:7
        bar(net,cc_k5_V_yeo7(c,net),'FaceColor',color_vecs(net,:))
    end
    set(gca,'XTick',1:7)
    set(gca,'XTickLabel',Yeo_names,'Fontsize',12)
    set(gca,'XTickLabelRotation',45)
    xlabel('Functional Networks','Fontsize',16)
    if c==2
        ylabel('Correlation','Fontsize',16)
    end
    ylim([-.4 1.1])
    box off
    for net=1:7
        if p_k5_V_yeo7(c,net)<(0.05/K)
            % To plot an asterisk when significant
            plot(net,cc_k5_V_yeo7(c,net)*1.1,'*k')
        end
    end
    end
    set(gca,'YGrid','on')
end
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

%%
saveas(gcf,[Directory 'Figures/Fig3/Figure3_Centroids_FNsv3Reviews' statTest Extension],'jpg')
saveas(gcf,[Directory 'Figures/Fig3/Figure3_Centroids_FNsv3Reviews' statTest Extension],'fig')



