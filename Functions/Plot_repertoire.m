function Plot_repertoire(Parcellation,N_areas,Extensions)

%%%%%%%%%
%
%  Function to plot the Repertoire obtained from a given clustering solution
%  of the LEiDA results.
%
%  Inputs are:
%  Parcellation = 'AAL116' or 'glasser378' or 'dbs80'
%  N_areas      = 90 or 116 or 80 or 378
%  Extensions   = [] or '_filtered' or '_contCORRECTED' or
%  '_contCORRECTEDfiltered' or _RL_filtered'
%
%  Ex.: Plot_repertoire('AAL116',90,[])
% 
%  It creates a new Figure with all K states, their probabilities, 
%  lifetimes and list of brain areas.
% 
%  Joana Cabral October 2019
%
%%%%%%%%%%%

Directory='/Users/joana/Documents/Work/LEiDA general/LEiDA_HCP/';
cd(Directory)
addpath(genpath(Directory))

load([Directory 'Centroids/LEiDA' num2str(N_areas) '_Centroids_V1' Extensions],'Centroids','P','LT','rangeK')

disp(' ')
disp('%%% PLOTS %%%%')
disp(['Choose number of clusters between ' num2str(rangeK(1)) ' and ' num2str(rangeK(end)) ])

K = input('Number of clusters: ');

% Get the K patterns
V=Centroids{rangeK==K}.C;

[cc_V_yeo7,p_V_yeo7] = Overlap_LEiDA_Yeo (Parcellation,N_areas,Extensions);
Yeo_names={'V','SM','DA','VA','L','FP','DM'};
Order=[1:2:N_areas N_areas:-2:2];

figure('Name',['LEiDA ' num2str(N_areas) ' K=' num2str(K) ' Centroids V1 ' Extensions])
colormap(jet) 
% Pannel A - Plot the BOLD PL patterns on the cortex (top and side) 
% Pannel B - Plot the FC patterns in matrix format
% Pannel C - Plot the probability of each state in each condition
  
for c=1:K
    subplot(6,K,c)
    % This needs function plot_nodes_in_cortex.m and aal_cog.txt
    plot_nodes_in_cortex_general(Parcellation,V(c,:))
    title({['State #' num2str(c)]})
    
    subplot(6,K,K+c)
    % This needs function plot_nodes_in_cortex.m and aal_cog.txt
    plot_nodes_in_cortex_general(Parcellation,V(c,:))
    view(0,0)
    
    subplot(6,K,2*K+c)
    FC_V=V(c,Order)'*V(c,Order);  
    li=max(abs(FC_V(:)));
    imagesc(FC_V,[-li li])   
    axis square
    title('FC pattern') 
    ylabel('Brain area #')
    xlabel('Brain area #')   
    
    subplot(6,K,3*K+c)  
    histogram(P(:,rangeK==K,c),0:0.05:1,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
    hold on
    plot([mean(P(:,rangeK==K,c)) mean(P(:,rangeK==K,c))],[0 100],':k')
    xlabel('Probability')
    ylabel('Nr of Subjects')
    xlim([0 1])
    ylim([0 60])
    box off
  
    subplot(6,K,4*K+c)  
    histogram(LT(:,rangeK==K,c),0:round(max(max(LT(:,rangeK==K,:)))/2)/20:max(max(LT(:,rangeK==K,:)))+1,'FaceColor',[.5 .5 1],'EdgeColor',[.5 .5 1]);
    hold on
    plot([mean(LT(:,rangeK==K,c)) mean(LT(:,rangeK==K,c))],[0 100],':k')
    xlabel('Lifetime (seconds)')
    ylabel('Nr of Subjects')
    xlim([0 max(max(LT(:,rangeK==K,:)))+1])
    ylim([0 70])
    box 
    
    subplot(6,K,5*K+c)  
    bar(1:7,squeeze(cc_V_yeo7(rangeK==K,c,:)),'FaceColor',[.5 0.8 .5],'EdgeColor',[.5 0.8 .5]);
    hold on
    set(gca,'XTickLabel',Yeo_names,'Fontsize',8)
    xlabel('Yeo Network')
    ylabel('Correlation')
    ylim([-1 1])
    box off 
    for net=1:7
        if p_V_yeo7(rangeK==K,c,net)<(0.05/K)
            plot(net,cc_V_yeo7(rangeK==K,c,net)*1.1,'*k')
        end
    end
            
    
end
    
    
