function LEiDA_Clouds_HCP
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)
% %
% % Code to visualize all LEiDA observations colored according to the k-means clustering solution
% % - plots each observation as a dot in 2 Principal Dimensions
% % - plots the cluster centroids as a +
% % - plots the corresponding matrices below
% % - plots the probability of each
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


Directory='/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/LEiDA_HCP/';
addpath(genpath(Directory))

% Define the dataset
N_areas=90;
Extension='_filtered_100unrelated'; %'_filtered_100unrelated'
typeRendering ='_filtered_100unrelated_5k'
Parcellation='AAL116';
plotDimension = 'plot3D'
% Load the Cluster Centroids for all K
load([Directory 'Centroids/LEiDA' num2str(N_areas) '_Centroids_V1' Extension],'Centroids','rangeK')

%% Define a 2-Dimensional Phase Space
% Ideally, we should obtain the 2 first principal components of the 
% covariance matrix of all the observations. However, there are too 
% many observations and there is not enough memory.
% So we find an alternative solution by taking the 2 PCs of the
% centroids for max K
%Var=cov(struct2array(Centroids{length(rangeK)}));
Var=cov(struct2array(Centroids{rangeK==5}));
[pc2, ~]=eigs(Var,3);
clear Var
% pc2 are 2 vectors of 1xN_areas, representing the x and y axis of the 
% 2D low-resolution landscape

cmap=[ .7 .7 .7; 0 0 1 ; 1 0 0 ; 1 0.5 0;  0 1 1; 1 0 1 ; 1 1 0];


% Make plots of the Voronoi Diagrams for 3 different K
K_selection=[5 10 20];

figure
for Partition=1:length(K_selection)
    k=K_selection(Partition);
    Centers=struct2array(Centroids{rangeK==k});
    dist_Centroids=zeros(2,k);
    for c=1:k
        V=Centers(c,:);
        dist_Centroids(1,c)=dot(V,pc2(:,1))/norm(V);%/norm(pc2(:,1));
        dist_Centroids(2,c)=dot(V,pc2(:,2))/norm(V);%/norm(pc2(:,2));
        % Note that the norm of PC1 and PC2 is 1
    end

    subplot(length(K_selection),3,2+3*(Partition-1))
    hold on
    voronoi(dist_Centroids(1,:),dist_Centroids(2,:));
    
    % To fill Voronoi cells
    if k==5
        switch typeRendering
        case '_unfiltered_100unrelated_7k'
            cells_K5={[2,4,5,7,8],[8,5,6],[3,6,5,4],[1,2,7],[1,2,4,3]};
        case '_unfiltered_100unrelated_5k'
            cells_K5={[4,6,3,2],[4,6,8,7],[3,5,8,6],[1,7,4,2],[1,5,3,2]};
        case '_filtered_100unrelated_7k'
            cells_K5={[2,4,5,7,8],[3,6,5,4],[1,2,7],[8,5,6],[1,2,4,3]};
        case '_filtered_100unrelated_5k'
            cells_K5={[4,6,3,2],[4,6,8,7],[3,5,8,6],[1,7,4,2],[1,5,3,2]};

        end
        [xv,yv]=voronoi(dist_Centroids(1,:),dist_Centroids(2,:));
        [vert(:,1), IA]=unique(xv); 
        vert(:,2)= yv(IA);
        for i = 1:k
            fill(vert(cells_K5{i},1),vert(cells_K5{i},2),cmap(i,:))
            plot(dist_Centroids(1,i),dist_Centroids(2,i),'.k','Markersize',5)
        end
    end
    switch typeRendering
        case '_unfiltered_100unrelated'   
        xlim([-0.4 1])
        ylim([0 1])
        case '_filtered_100unrelated'
    end
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    %view(23,17)
    %axis off
    title(['K = ' num2str(k)])
    set(gca,'DataAspectRatio',[1 2 1])
end




% NOW LOAD ALL THE EIGENVECTORS FROM ALL SESSIONS
load(['LEiDAaal' num2str(N_areas) '_LR' Extension],'V1_all')

% IF USING THE CONTINUOUS V1 USE THIS CODE
%load(['LEiDAaal' num2str(N_areas) '_V1' Extension(2:end)],'V1previous_all')
% V1_all=V1previous_all;
% clear V1previous_all;

%% Plot the trajectory in phase space from single subject colored as a function of time
subplot(2,3,1)
T_subject=1200-2;
cmap_time=winter(T_subject);
sub=1;
hold on
for t=(sub-1)*(1200-2)+1:sub*(1200-2)   %1:length(IDX)
    V=V1_all(t,:);
    dist1=dot(V,pc2(:,1))/norm(V);
    dist2=dot(V,pc2(:,2))/norm(V);
    plot(dist1 ,dist2,'.','Markersize',10,'Color',cmap_time(t,:))
end
xlabel('1st PC')
ylabel('2nd PC')
    xlim([-1 1])
    ylim([-1 1])
title(['Trajectory from one fMRI experiment (subject ' num2str(sub) ')'])
colormap(winter)
colorbar
set(gca,'DataAspectRatio',[1 2 1])
clear cmap_time

%% Plot all observations in phase space from all subjects
subplot(2,3,4)
hold on
switch plotDimension
    case 'plot2D'
        for t=1:10000 %size(V1_all,1)
            V=V1_all(t,:);
            dist1=dot(V,pc2(:,1))/norm(V);
            dist2=dot(V,pc2(:,2))/norm(V);
            plot(dist1 ,dist2,'.k','Markersize',2)
        end
        xlabel('1st PC')
        ylabel('2nd PC')
            xlim([-1 1])
            ylim([-1 1])
        title(['Observations from all fMRI experiments (99 subjects, total ' num2str(size(V1_all,1))  ' timepoints)'])
        set(gca,'DataAspectRatio',[1 2 1])
    case 'plot3D'
        for t=1:10000 %size(V1_all,1)
            V=V1_all(t,:);
            dist1=dot(V,pc2(:,1))/norm(V);
            dist2=dot(V,pc2(:,2))/norm(V);
            dist3=dot(V,pc2(:,3))/norm(V);
            plot3(dist1 ,dist2,dist3,'.k','Markersize',2)
            view(3)
        end
        xlabel('1st PC')
        ylabel('2nd PC')
        zlabel('3rd PC')
        xlim([-1 1])
        ylim([-1 1])
        zlim([-1 1])
        title(['Observations from all fMRI experiments (99 subjects, total ' num2str(size(V1_all,1))  ' timepoints)'])
        set(gca,'DataAspectRatio',[1 2 1])
end

% Load the cluster time courses for k=5
k=5;
load([Directory 'Centroids/LEiDA' num2str(N_areas) '_Centroids_V1' Extension],'Kmeans_results')
IDX=Kmeans_results{4}.IDX;

%% Plot the trajectory in phase space from single subject colored as a function of cluster
subplot(2,3,3)
T_subject=1200-2;
sub=1;
hold on
switch plotDimension
    case 'plot2D'
        for t=(sub-1)*(1200-2)+1:sub*(1200-2)   %1:length(IDX)
            V=V1_all(t,:);
            dist1=dot(V,pc2(:,1))/norm(V);
            dist2=dot(V,pc2(:,2))/norm(V);
            plot(dist1 ,dist2,'.','Markersize',10,'Color',cmap(IDX(t),:))
        end
        xlabel('1st PC')
        ylabel('2nd PC')
            xlim([-1 1])
            ylim([-1 1])
        title(['Trajectory from one fMRI experiment (subject ' num2str(sub) ')'])
        set(gca,'DataAspectRatio',[1 2 1])
    case 'plot3D'
        for t=(sub-1)*(1200-2)+1:sub*(1200-2)   %1:length(IDX)
            V=V1_all(t,:);
            dist1=dot(V,pc2(:,1))/norm(V);
            dist2=dot(V,pc2(:,2))/norm(V);
            dist3=dot(V,pc2(:,3))/norm(V);
            plot3(dist1, dist2, dist3,'.','Markersize',10,'Color',cmap(IDX(t),:))
        end
        xlabel('1st PC')
        ylabel('2nd PC')
        zlabel('3rd PC')
        xlim([-1 1])
        ylim([-1 1])
        zlim([-1 1])
        title(['Trajectory from one fMRI experiment (subject ' num2str(sub) ')'])
        set(gca,'DataAspectRatio',[1 2 1])
end

%% Plot all observations in phase space from all subjects
subplot(2,3,6)
hold on
switch plotDimension
    case 'plot2D'
        for t=1:10000 %size(V1_all,1)
            V=V1_all(t,:);
            dist1=dot(V,pc2(:,1))/norm(V);
            dist2=dot(V,pc2(:,2))/norm(V);
            plot(dist1 ,dist2,'.k','Markersize',2,'Color',cmap(IDX(t),:))
        end
        xlabel('1st PC')
        ylabel('2nd PC')
            xlim([-1 1])
            ylim([-1 1])
        title(['Observations from all fMRI experiments (99 subjects, total ' num2str(size(V1_all,1))  ' timepoints)'])
        set(gca,'DataAspectRatio',[1 2 1])
        
    case 'plot3D'
        for t=1:10000 %size(V1_all,1)
            V=V1_all(t,:);
            dist1=dot(V,pc2(:,1))/norm(V);
            dist2=dot(V,pc2(:,2))/norm(V);
            dist3=dot(V,pc2(:,3))/norm(V);
            plot3(dist1 ,dist2,dist3,'.k','Markersize',2,'Color',cmap(IDX(t),:))
            view(3)

        end
        xlabel('1st PC')
        ylabel('2nd PC')
        zlabel('3nrd PC')
        xlim([-1 1])
        ylim([-1 1])
        zlim([-1 1])
        title(['Observations from all fMRI experiments (99 subjects, total ' num2str(size(V1_all,1))  ' timepoints)'])
        set(gca,'DataAspectRatio',[1 2 1])
end

%%
% 

% Supplementary Figure, plot the trajectories for 20 subjects
figure('name',['Cluster Clouds for AAL' num2str(N_areas) ' K' num2str(k) Extension],'Color','w')

Subjects=49:96;

for s=1:length(Subjects)
    sub=Subjects(s);
    subplot_tight(8,6,s,0.01)
    
    hold on
    for t=(sub-1)*(1200-2)+1:sub*(1200-2)   %1:length(IDX)
        V=V1_all(t,:);
        dist1=dot(V,pc2(:,1))/norm(V);
        dist2=dot(V,pc2(:,2))/norm(V);
        plot(dist1 ,dist2,'.','Markersize',3,'Color',cmap(IDX(t),:))
    end
    
    % Then plot the Centroids in the 2-dimensional space
%     for c=1:k
%         % Note that the norms of PC1 and PC2 is 1
% 
%         plot(dist_Centroids(1,:),dist_Centroids(2,:),'.','Markersize',50,'Color','k');
%         plot(dist_Centroids(1,:),dist_Centroids(2,:),'*','Markersize',10,'Color',cmap(c,:),'Linewidth',2);
%     end

    title(['Subject ' num2str(sub)])
    box on
    xlim([-1 1])
    ylim([-1 1])
    set(gca,'DataAspectRatio',[1 2 1])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    
%     if s==41 
%         xlabel('1st PC')
%         ylabel('2nd PC')
%     end
end

%kho;

% 
% 
% 
% 
% YeoColor = [120 18 134; 70 30 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78]./256;
% 

% 
% % CODE TO DEFINE COLOR OF VORONOI CELLS
% k=5;
% Centers=struct2array(Centroids{rangeK==k});
% dist_Centroids=zeros(2,k);
% for c=1:k
%     V=Centers(c,:);
%     dist_Centroids(1,c)=dot(V,pc2(:,1))/norm(V);%/norm(pc2(:,1));
%     dist_Centroids(2,c)=dot(V,pc2(:,2))/norm(V);%/norm(pc2(:,2));
%     % Note that the norm of PC1 and PC2 is 1
% end
%     
% clear vert
% [xv,yv]=voronoi(dist_Centroids(1,:),dist_Centroids(2,:));
% [vert(:,1), IA]=unique(xv); 
% vert(:,2)= yv(IA);
% figure
% hold on
% for n=1:length(vert)
%     plot(vert(n,1),vert(n,2),'*k')
%   text(vert(n,1)+.1,vert(n,2),num2str(n))
% end
% 
% clear c
% for n=1:k
%     c{n}=input(['Boundaries of cell ' num2str(n) ':']);
%     cmap(n,:)=input(['Color code for cell ' num2str(n) ':']);
% end
% for i = 1:length(c)
%     fill(vert(c{i},1),vert(c{i},2),cmap(i,:))
%     plot(dist_Centroids(1,i),dist_Centroids(2,i),'*k','Markersize',10)
% end