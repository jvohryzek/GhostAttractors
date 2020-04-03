% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)
% %
% % Function to plot Figure 5
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
%% load data
Directory='/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/';
addpath(genpath(Directory))

% Define the dataset
N_areas=90;
Extension='_unfiltered_100unrelated'; %'_filtered_100unrelated';
Parcellation='AAL116';

% Load the Cluster Centroids for all K
load([Directory 'LEiDA_HCP/Centroids/LEiDA' num2str(N_areas) '_Centroids_V1' Extension])

% just to get the vornoi graph right
Extension_voronoi= '_unfiltered';
load([Directory 'LEiDA_HCP/Centroids/LEiDA' num2str(N_areas) '_Centroids_V1' Extension_voronoi],'Centroids','rangeK')

Var=cov(struct2array(Centroids{rangeK==5}));
[pc2, ~]=eigs(Var,2);
clear Var
cmap=[ .7 .7 .7; 0 0 1 ; 1 0 0 ; 1 0.5 0;  0 1 1; 1 0 1 ; 1 1 0];

%%
for i=1:99
    transMat(i,:,:) = PT{i,4};
    transMatnorm(i,:,:) = PTnorm{i,4};

    probMat(i,:) = nonzeros(P(i,4,:))';

end
transMatAll = squeeze(mean(transMat));
transMatAllnorm = squeeze(mean(transMatnorm));

probMatAll = mean(probMat,1);
%%
nodeNames =  {'State 1', 'State 2', 'State 3', 'State 4', 'State 5'};

G= digraph(round(transMatAllnorm,2))
LWidths = 15*G.Edges.Weight/max(G.Edges.Weight);
k1 = figure,

subplot(6,9,[6,7,8,9,15,16,17,18,24,25,26,27,33,34,35,36,42,43,44,45,51,52,53,54])
% p = plot(G,'Layout','force','EdgeLabel',G.Edges.Weight,'LineWidth',LWidths,...
%   'NodeFontWeight','bold','EdgeFontWeight','bold')
p = plot(G,'Layout','force','EdgeLabel',G.Edges.Weight,'LineWidth',LWidths)
%view([0 90])
p.Marker = 'o';
p.NodeColor = [0.6350 0.0780 0.1840];%'b';
p.EdgeColor = 'k';%[0.6350 0.0780 0.1840];
p.MarkerSize = 12;
p.ArrowSize = 10;
p.NodeLabel = nodeNames;
p.NodeFontSize = 18;
p.EdgeFontSize = 18;
p.EdgeLabelColor = [0.6350 0.0780 0.1840];
p.NodeLabelColor = 'k';%[0.6350 0.0780 0.1840];
%axis square
axis off
title('Transition Graph')
set(gca,'FontSize',20,'FontWeight','bold')

subplot(8,9,[13,14,22,23,31,32,40,41])
imagesc(transMatAllnorm)
for j1 = 1:5
    for j2 = 1:5
        caption = sprintf('%.2f', transMatAllnorm(j1,j2));
        text(j2-0.35, j1, caption, 'FontSize', 15, 'Color', [1, 1, 1],'FontWeight','bold');
    end
end
xticks(1:1:5); xticklabels({'State 1','State 2','State 3','State 4','State 5'})
yticks(1:1:5); yticklabels({'State 1','State 2','State 3','State 4','State 5'})
xtickangle(45);ytickangle(45)
set(gca,'FontSize',20,'FontWeight','bold')
axis square
ylabel('From')
xlabel('To')
subplot(6,9,[4,5])
axis off
title('Transition Matrix')

set(gca,'FontSize',20,'FontWeight','bold')

subplot(8,9,[10,11,19,20,28,29,37,38])

k=5;
Centers=struct2array(Centroids{rangeK==k});
dist_Centroids=zeros(2,k);
for c=1:k
    V=Centers(c,:);
    dist_Centroids(1,c)=dot(V,pc2(:,1))/norm(V);%/norm(pc2(:,1));
    dist_Centroids(2,c)=dot(V,pc2(:,2))/norm(V);%/norm(pc2(:,2));
    % Note that the norm of PC1 and PC2 is 1
end

voronoi(dist_Centroids(1,:),dist_Centroids(2,:));

% To fill Voronoi cells
caption = {'State 1', 'State 2', 'State 3', 'State 4', 'State 5'};

if k==5
    cells_K5={[3,5,6,4],[2,4,3,1],[1,3,5,7],[7,5,6,8],[8,6,4,2]};
    [xv,yv]=voronoi(dist_Centroids(1,:),dist_Centroids(2,:));
    [vert(:,1), IA]=unique(xv); 
    vert(:,2)= yv(IA);
    for i = 1:k
        fill(vert(cells_K5{i},1),vert(cells_K5{i},2),cmap(i,:))
        hold on
        plot(dist_Centroids(1,i),dist_Centroids(2,i),'.k','Markersize',10)
        text(dist_Centroids(1,i)-0.15, dist_Centroids(2,i)-0.15, caption{i}, 'FontSize', 16, 'Color', [0, 0, 0],'FontWeight','bold');

    end
end


xlim([-1 1])
ylim([-1 1])
set(gca,'XTick',[])
set(gca,'YTick',[])
axis square
%view(23,17)
%axis off
set(gca,'FontSize',20,'FontWeight','bold')

subplot(6,9,[1,2])
axis off
title('Decomposition')


set(gca,'FontSize',20,'FontWeight','bold')
% set(gca,'DataAspectRatio',[1 2 1])

set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

saveas(k1,[Directory 'Figures/Fig5/graphRepresentationFigure5' Extension],'jpg')
saveas(k1,[Directory 'Figures/Fig5/graphRepresentationFigure5' Extension],'eps')
%% histogram for all the transition 

load([Directory '/LEiDA_HCP/K5_results_LR_RL' Extension '.mat'])
n_Subjects = 99;
for i=1:n_Subjects
    tmpLR = squeeze(PTnorm_K5_LR(i,:,:))';
    tmp2LR = tmpLR(:);
    tmpRL = squeeze(PTnorm_K5_RL(i,:,:))';
    tmp2RL = tmpRL(:);
    for j=1:25
        
        histPTnorm_K5_LR(i,j)=tmp2LR(j);
        histPTnorm_K5_RL(i,j)=tmp2RL(j);
    end
end
%%
fSI3 = figure
m=0;n=0;
for k=1:25
    subplot(5,5,k)
    max_PTnorm=max(max([histPTnorm_K5_LR(:,k) histPTnorm_K5_RL(:,k)]));

    histogram(histPTnorm_K5_LR(:,k),0:round(max_PTnorm,1)/20:max_PTnorm,'Normalization','Probability','FaceColor',[0 .5 .2],'EdgeColor',[0 .5 .2]);
    hold on
    histogram(histPTnorm_K5_RL(:,k),0:round(max_PTnorm,1)/20:max_PTnorm,'Normalization','Probability','FaceColor',[.5 .8 .4],'EdgeColor',[.5 .8 .4]);

%     title(sprintf('Edge %d',k))

    if k==5
        legend({'LR fMRI session','RL fMRI session'},'location','northeast')
        ylabel('Number of Subjects')
        xlabel('Probability of Transition')
    end
    
    listStates2 = [21 22 23 24 25];
    if mod(k,5) == 1
        m = m+1;
        ylabel(sprintf('State %d',listStates2(m)-20));
    end
    if sum(listStates2 == k)
        n = n+1;
        xlabel(sprintf('State %d',listStates2(n)-20));
    end
    
    xlim([0 1])
    set(gca,'FontSize',14,'FontWeight','bold')
end
[ax1,h1]=suplabel('Switch From','y');set(h1,'FontSize',30)
[ax2,h2]=suplabel('Switch To');set(h2,'FontSize',30)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
saveas(fSI3,[Directory 'Figures/Fig5/IndSbjFigure5' Extension],'jpg')
saveas(fSI3,[Directory 'Figures/Fig5/IndSbjFigure5' Extension],'eps')

%%
k1sbj = figure;%('Renderer', 'painters', 'Position', [10 10 9000 6000])
addSbj = 10
for i=1:10
    G= digraph(round(squeeze(PTnorm_K5_LR(i+addSbj,:,:)),2))
    LWidths = 15*G.Edges.Weight/max(G.Edges.Weight);
    subplot(2,5,i)
    p = plot(G,'Layout','force','EdgeLabel',G.Edges.Weight,'LineWidth',LWidths,...
        'NodeFontWeight','bold','EdgeFontWeight','bold')
    %view([0 90])
    p.Marker = 'o';
    p.NodeColor = [0.6350 0.0780 0.1840];%'b';
    p.EdgeColor = 'k';%[0.6350 0.0780 0.1840];
    p.MarkerSize = 12;
    p.ArrowSize = 10;
    p.NodeLabel = nodeNames;
    p.NodeFontSize = 8;
    p.EdgeFontSize = 8;
    p.EdgeLabelColor = [0.6350 0.0780 0.1840];
    p.NodeLabelColor = 'k';%[0.6350 0.0780 0.1840];
   % axis square
    axis off
    title(sprintf('Subject %d',i+addSbj))
    set(gca,'FontSize',16,'FontWeight','bold')
end

set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
saveas(k1sbj,[Directory 'Figures/Fig5/Figure5TransitionGraphIndSbj11_20' Extension],'jpg')
saveas(k1sbj,[Directory 'Figures/Fig5/Figure5TransitionGraphIndSbj11_20' Extension],'fig')

%%
% % % for i=1:n_Subjects
% % %     tmp = PTnorm{i,4}';
% % %     tmp2 = tmp(:);
% % %     for j=1:25
% % %         histPT(i,j)=tmp2(j);
% % %     end
% % % end
% % % fSI3 = figure
% % % for k=1:25
% % %     subplot(5,5,k)
% % %     hist(histPT(:,k),10)
% % %     title(sprintf('Edge %d',k))
% % %     ylabel('# of Subjects')
% % %     xlabel('# of tps')
% % %     xlim([0 1])
% % % 
% % % end