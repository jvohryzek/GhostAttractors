function Plot_repertoire_RL_LR

%%%%%%%%%
%
%  Function to plot Figure 4
%  Plots the repertoire of Centroids obrained for the LR session
%  Plots the probabilities, Lifetimes for the 2 sessions.
%
%  It creates a new Figure with all K states, their probabilities,
%  lifetimes and list of brain areas.
%
%  Joana Cabral November 2019
%
%%%%%%%%%%%

Directory='/Users/joana/Documents/Work/LEiDA general/LEiDA_HCP/';
addpath(genpath(Directory))

load([Directory 'AAL90_filtered_K5_results_LR_RL'],'Centroids_K5_LR',...
    'P_K5_LR','P_K5_RL','LT_K5_LR','LT_K5_RL')
K=5;
Parcellation= 'AAL116';
N_areas=90;
n_Subjects=95;

[cc_V_yeo7,p_V_yeo7,rangeK] = Overlap_LEiDA_Yeo (Parcellation,N_areas,'_filtered',0);
cc_k5_V_yeo7=squeeze(cc_V_yeo7(rangeK==K,1:K,:));
p_k5_V_yeo7=squeeze(p_V_yeo7(rangeK==K,1:K,:));

color_vecs = [120 18 134; 70 30 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78]./256;

%Yeo_names={'Visual','Somatomotor','Dorsal Att.','Ventral Att.','Limbic','Frontoparietal','Default Mode'};
Yeo_names={'VIS','SMT','DAT','VAT','LBC','FPA','DMN'};


Order=[1:2:N_areas N_areas:-2:2];

% LEiDA networks colored according to closest RSN
Volume=struct2array(load('ParcelsMNI2mm',['V_' Parcellation])); 
Brain_Mask=niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex=smooth3(Brain_Mask>0);

figure('Name','LEiDA 90 K=5 Centroids V1 filtered')
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
    axis on
    
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
    title({['Centroid #' num2str(c)]})
    material dull; lighting gouraud
    daspect([1 1 1])
    view(-90,90)    % Top view    Side view:   view(0,0)     
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis on
    
    %% Same but view from the side
    % First Plot the arrows view from top
    subplot_tight(6,2*K,c*2-1+2*K,0.01)
    plot_arrows_in_cortex_general(Parcellation,Centroids_K5_LR(c,:))
    material dull; lighting gouraud
    daspect([1 1 1])
    view(0,0)    % Top view    Side view:   view(0,0)     
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis on
    
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
    axis on
    
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
        set(gca,'YTick',20:20:80,'Fontsize',8)

    end
    xlabel('Brain area #')
    set(gca,'XTick',20:20:80,'Fontsize',8)

    
    subplot_tight(6,K,3*K+c,0.06)
    hold on
    % P=[squeeze(P_K5_LR(:,c)) squeeze(P_K5_RL(:,c))];
    histogram(P_K5_LR(:,c),0:0.05:1,'FaceColor',[1 .7 .1]);
    histogram(P_K5_RL(:,c),0:0.05:1,'FaceColor',[.8 .3 .1]);
    
%     [t, P_val_2sessions]=ttest2(squeeze(P_K5_LR(:,c)),squeeze(P_K5_RL(:,c)),0.05/K)
% %     title(['p=' num2str(P_val_2sessions*K)])
% %     disp(['P p=' num2str(P_val_2sessions*K)])
    
    %plot([mean(P(:,rangeK==K,c)) mean(P(:,rangeK==K,c))],[0 100],':k')
    xlabel('Probability')
    if c==1
        ylabel('Subjects #','Fontsize',10)
    end
    xlim([0 1])
    ylim([0 40])
    box off
    if c==1
        legend({'LR fMRI session','RL fMRI session'})
    end
    
    subplot_tight(6,K,4*K+c,0.06)
    hold on
    max_Lt=max(max([LT_K5_LR(:) LT_K5_RL(:)]));
    histogram(LT_K5_LR(:,c),0:round(max_Lt)/20:max_Lt,'FaceColor',[.1 .7 1]);
    histogram(LT_K5_RL(:,c),0:round(max_Lt)/20:max_Lt,'FaceColor',[.1 .3 .8]);
    if c==1
        legend({'LR fMRI session','RL fMRI session'})
    end
    %    plot([mean(LT(:,rangeK==K,c)) mean(LT(:,rangeK==K,c))],[0 100],':k')
    xlabel('Dwell time (secs)')
    if c==1
        ylabel('Subjects #','Fontsize',10)
    end
    
%     [t, P_val_2sessions]=ttest2(squeeze(LT_K5_LR(:,c)),squeeze(LT_K5_RL(:,c)),0.05/K)
%     title(['p=' num2str(P_val_2sessions*K)])
%     disp(['LT p=' num2str(P_val_2sessions*K)])
    
    xlim([0 max_Lt])
    ylim([0 45])
    box off
    
    subplot_tight(6,K,5*K+c,0.06)
    hold on
    for net=1:7
        bar(net,cc_k5_V_yeo7(c,net),'FaceColor',color_vecs(net,:))
    end
    set(gca,'XTick',1:7)
    set(gca,'XTickLabel',Yeo_names,'Fontsize',8)
    set(gca,'XTickLabelRotation',45)
    xlabel('Functional Networks','Fontsize',10)
    if c==2
        ylabel('Correlation','Fontsize',10)
    end
    ylim([-.4 1.1])
    box off
    for net=1:7
        if p_k5_V_yeo7(c,net)<(0.05/K)
            % To plot an asterisk when significant
            plot(net,cc_k5_V_yeo7(c,net)*1.1,'*k')
        end
    end
    set(gca,'YGrid','on')
end


