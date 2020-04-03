% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % INTER-CLASS CORRELATION (ICC) ANALYSIS
% %
% % This function compute ICC for different measures of LEiDA% %
% % 
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

%% LOAD DATA
% LR
Directory='/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/';

Extension= '_unfiltered_100unrelated';% '_filtered_100unrelated';

load([Directory '/LEiDA_HCP/K5_results_LR_RL' Extension '.mat'])
n_Subjects = 99;
%% ICC
for p =1:5
    %P_ICC(p) = ICC(2,'single',[P_K5_LR(:,p), P_K5_RL(:,p)]);
    P_ICC(p) = ICC(1,'single',[P_K5_LR(:,p), P_K5_RL(:,p)]);

    %LT_ICC(p) = ICC(2,'single',[LT_K5_LR(:,p), LT_K5_RL(:,p)]);
    LT_ICC(p) = ICC(1,'single',[LT_K5_LR(:,p), LT_K5_RL(:,p)]);

    for r = 1:5
        %PTnorm_ICC(p,r) = ICC(2,'single',[PTnorm_K5_LR(:,p,r), PTnorm_K5_RL(:,p,r)]);
        %PT_ICC(p,r) = ICC(2,'single',[PT_K5_LR(:,p,r), PT_K5_RL(:,p,r)]);
        
        PTnorm_ICC(p,r) = ICC(1,'single',[PTnorm_K5_LR(:,p,r), PTnorm_K5_RL(:,p,r)]);
        PT_ICC(p,r) = ICC(1,'single',[PT_K5_LR(:,p,r), PT_K5_RL(:,p,r)]);

    end
end

%% PLOTTING PRnorm
k1 = figure,
imagesc(PTnorm_ICC)
colormap('parula')
for j1 = 1:5
    for j2 = 1:5
        caption = sprintf('%.2f', PTnorm_ICC(j1,j2));
        % setting all negative values to 0
        if PTnorm_ICC(j1,j2) <0
            caption = sprintf('%.2f', 0);
        end
        text(j2-0.2, j1, caption, 'FontSize', 15, 'Color', [1, 1, 1],'FontWeight','bold');
    end
end
axis square
colorbar
title('ICC: Probability of Transition')
xticks(1:1:5); xticklabels({'State 1','State 2','State 3','State 4','State 5'})
yticks(1:1:5); yticklabels({'State 1','State 2','State 3','State 4','State 5'})

set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
set(gca,'FontSize',18);
set(findall(gcf,'-property','FontSize'),'FontSize',24)

saveas(k1,[Directory '/Figures/Fig6/MatrixPlotsPTFigure6Revision' Extension],'jpg')
saveas(k1,[Directory '/Figures/Fig6/MatrixPlotsPTFigure6Revision' Extension],'fig')

%% Spider plots

AxesLimits = [0, 0, 0, 0, 0; 1, 1, 1, 1, 1]; % Axes limits [min axes limits; max axes limits]
AxesLabels = {'State 1', 'State 2', 'State 3', 'State 4', 'State 5'}; % Axes properties

% SPIDER PLOT FRACTIONAL OCCUPANCY

spider_plot(P_ICC,'AxesLimits', AxesLimits,'AxesLabels', AxesLabels,...
    'LabelFontSize', 18, 'LineWidth', 4,'Marker', 'none',...
    'FillOption', 'on','Color', [139, 0, 0]/255,...
    'AxesInterval', 5,'AxesPrecision', 2,...
    'AxesDisplay', 'one','FillTransparency', 0.1);
title('ICC: Fractional Occupancy')
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
saveas(gcf,[Directory '/Figures/Fig6/SpiderPlotsFOFigure6Revision' Extension],'jpg')
saveas(gcf,[Directory '/Figures/Fig6/SpiderPlotsFOFigure6Revision' Extension],'fig')

% SPIDER PLOT LIFE TIMES
spider_plot(LT_ICC,'AxesLimits', AxesLimits,'AxesLabels', AxesLabels,...
    'LabelFontSize', 18, 'LineWidth', 4,'Marker', 'none',...
    'FillOption', 'on','Color', [139, 0, 0]/255,...
    'AxesInterval', 4,'AxesPrecision', 2,...
    'AxesDisplay', 'one','FillTransparency', 0.1);
title('ICC: Dwell Time')
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

saveas(gcf,[Directory '/Figures/Fig6/SpiderPlotsLTFigure6Revision' Extension],'jpg')
saveas(gcf,[Directory '/Figures/Fig6/SpiderPlotsLTFigure6Revision' Extension],'fig')
%% plotting scatter plots
figureMarker = {'ks', 'k*', 'ko', 'kx', 'kd'}
figureStates = {'State 1', 'State 2', 'State 3', 'State 4', 'State 5'}
fSI10 = figure
for i=1:5
    subplot(2,5,i)
    plot(squeeze(P_K5_LR(:,i))./max(squeeze(P_K5_LR(:,i))),squeeze(P_K5_RL(:,i))./max(squeeze(P_K5_RL(:,i))),'ks')
    hold on
    plot(0:0.1:1,0:0.1:1,'r','LineWidth',2)
    title(figureStates{i});
    if i==1
        ylabel('RL session: Fractional Occupancy');xlabel('LR session: Fractional Occupancy')
    end
    axis square
    
    subplot(2,5,i+5)
    plot(squeeze(LT_K5_LR(:,i))./max(squeeze(LT_K5_LR(:,i))),squeeze(LT_K5_RL(:,i))./max(squeeze(LT_K5_RL(:,i))),'ks')
    hold on
    plot(0:0.1:1,0:0.1:1,'r','LineWidth',2)
    title(figureStates{i});
    if i==1
        ylabel('RL session: Dwell Times');xlabel('LR session: Dwell Times');
    end
    axis square
end
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

saveas(fSI10,[Directory '/Figures/Fig6/LRRLCorrelationLTPFigure6SI' Extension],'jpg')
saveas(fSI10,[Directory '/Figures/Fig6/LRRLCorrelationLTPFigure6SI' Extension],'fig')

%% Plotting
 % perhaps to play in the future with plotting the ICC MSEb and MSEw
% % % covSbj = zeros(5,n_Subjects,n_Subjects);
% % % for p=1:5
% % %     for i =1:n_Subjects
% % %         for j=1:n_Subjects
% % %             covSbj(p,i,j) = [P_K5_LR(i,p) - P_K5_RL(j,p)];
% % %         end
% % %     end
% % % end
% % % 
% % % figure,
% % % for i = 1:5
% % %     subplot(1,5,i)
% % %     histogram(nonzeros(triu(squeeze(covSbj(i,:,:)),1)),20,'Normalization','Probability','BinWidth',0.05,'BinLimits',[-0.5,0.5],'EdgeColor','none','FaceColor',[1 0 0])
% % %     hold on
% % %     histogram(diag(squeeze(covSbj(i,:,:))),20,'Normalization','Probability','BinWidth',0.05,'BinLimits',[-0.5,0.5],'EdgeColor','none','FaceColor',[0 0 1])
% % % end