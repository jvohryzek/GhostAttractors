function LEiDA_HCP_parsing(parse) %(mode,parcellation)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)
% %
% % This function clusterises and runs in parallel the LeadingEigs
% %
% %  Choose mode:
% % 'kmeans'                => runs kmeans for individual k-solutions for 
% %                            subsequent analysis in 
% %                            LEiDA_HCP_standalone.m (mode: 'parse')
% % 'kmeans_filtrationBand' => runs kmeans for different fequency bands
% %                            with K=5 for subsequent analysis in 
% %                            LEiDA_HCP_standalone.m (mode: 'parse_filterBand')
% % 'louvain'               => runs louvain detection algortihm in batches
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
% Orginal version
% Joana Cabral March 2017
% joanacabral@med.uminho.pt
% Cabral, et al. 2017 Scientific reports 7, no. 1 (2017): 5135.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mode = 'kmeans_filtrationBand'; % 'kmeans' 'louvain'

%% 2 - Cluster the Leading Eigenvectors
Directory  = '/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/';
addpath('/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/LEiDA_Codes/BCT')
load('/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/LEiDA_HCP/Eigenvectors/LEiDAaal90_LR_unfiltered_100unrelated.mat')

disp('Clustering the eigenvectors into')
% Leading_Eig is a matrix containing all the eigenvectors:
% Collumns: N_areas are brain areas (variables)
% Rows: (Tmax-2)*n_Subjects are all time points (independent observations)

% Set maximum/minimum number of clusters
% There is no fixed number of FC states that the brain can display
% Keep the range small for the first trials
% Extend the range depending on the hypothesis of each work
if strcmp(mode,'kmeans') % compute LEiDA   
    % Set the parameters for Kmeans clustering
    Kmeans_results_single=cell(1);

    disp(['- ' num2str(parse) ' FC states'])
    distance = 'Cosine';
    %%%%%%
    %%
    [IDX, C, SUMD, D] = kmeans(V1_all, parse, 'Distance',distance, 'Replicates',100,'MaxIter',400,'Display','final'); %,'Options',statset('UseParallel',1));
    [~, ind_sort]=sort(hist(IDX,1:parse),'descend');
    [~,idx_sort]=sort(ind_sort,'ascend');
    Kmeans_results_single{1}.IDX=idx_sort(IDX);   % Cluster time course - numeric collumn vectors
    Kmeans_results_single{1}.C=C(ind_sort,:);       % Cluster centroids (FC patterns)
    Kmeans_results_single{1}.SUMD=SUMD(ind_sort); % Within-cluster sums of point-to-centroid distances
    Kmeans_results_single{1}.D=D(:,ind_sort);       % Distance from each point to every centroid   end
    clust = idx_sort(IDX);

    %% Evaluate Clustering performance

    % 2. Silhouette

    if strcmp(distance,'euclidean')
        eva = evalclusters(V1_all',clust','Silhouette','Distance','Euclidean');
        evaSilh = eva.CriterionValues;
    elseif strcmp(distance,'cityblock')
       eva = evalclusters(V1_all',clust','Silhouette','Distance','cityblock');
       evaSilh = eva.CriterionValues;
    elseif strcmp(distance,'Cosine')
       eva = evalclusters(V1_all,clust','Silhouette','Distance','cosine');
       evaSilh = eva.CriterionValues;
    elseif strcmp(distance,'correlation')
       eva = evalclusters(V1_all',clust','Silhouette','Distance','correlation');
       evaSilh = eva.CriterionValues;
    end

    % 3. Davies-Bouldin => for euclidean distance only
    % [1] Davies, D. L., and D. W. Bouldin. ?A Cluster Separation Measure.? IEEE Transactions on Pattern Analysis and Machine Intelligence. Vol. PAMI-1, No. 2, 1979, pp. 224?227.

    if strcmp(distance,'euclidean')

        eva = evalclusters(Kmeans_results_single{1}.IDX',clust','DaviesBouldin');
        evaDB(:,it) = eva.CriterionValues;
    else
        evaDB = [];
    end

    % save Kmeans results
    clear C D s t seed t_all idx_sort SUMD IDX
    save([Directory '/LEiDA_HCP/V1_LR_aal90_unfiltered_100unrelated/LEiDA' num2str(N_areas) '_K' num2str(parse) 'results.mat'], 'Kmeans_results_single', 'evaSilh', 'N_areas', 'n_Subjects','TR') 
 elseif strcmp(mode,'kmeans_filtrationBand')
     filterBand = {'0.01-0.07', '0.01-0.1', '0.01-0.2', '0.01-0.6'};

     Directory  = '/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/LEiDA_HCP/';
     addpath('/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/LEiDA_Codes/BCT')
     load(strcat(Directory,'Eigenvectors/LEiDAaal90_LR_filtered_',filterBand{parse},'_100unrelated.mat'))

      % Set the parameters for Kmeans clustering
      Kmeans_results_single=cell(1);

      disp(['filtering ' filterBand{parse} ' Hz'])
      distance = 'Cosine';
      %%
      [IDX, C, SUMD, D] = kmeans(V1_all, 5, 'Distance',distance, 'Replicates',100,'MaxIter',400,'Display','final'); %,'Options',statset('UseParallel',1));
      [~, ind_sort]=sort(hist(IDX,1:5),'descend');
      [~,idx_sort]=sort(ind_sort,'ascend');
      Kmeans_results_single{1}.IDX=idx_sort(IDX);   % Cluster time course - numeric collumn vectors
      Kmeans_results_single{1}.C=C(ind_sort,:);       % Cluster centroids (FC patterns)
      Kmeans_results_single{1}.SUMD=SUMD(ind_sort); % Within-cluster sums of point-to-centroid distances
      Kmeans_results_single{1}.D=D(:,ind_sort);       % Distance from each point to every centroid   end
      clust = idx_sort(IDX);

      %% Evaluate Clustering performance

      % 2. Silhouette

        if strcmp(distance,'euclidean')
            eva = evalclusters(V1_all',clust','Silhouette','Distance','Euclidean');
            evaSilh = eva.CriterionValues;
        elseif strcmp(distance,'cityblock')
           eva = evalclusters(V1_all',clust','Silhouette','Distance','cityblock');
           evaSilh = eva.CriterionValues;
        elseif strcmp(distance,'Cosine')
           eva = evalclusters(V1_all,clust','Silhouette','Distance','cosine');
           evaSilh = eva.CriterionValues;
        elseif strcmp(distance,'correlation')
           eva = evalclusters(V1_all',clust','Silhouette','Distance','correlation');
           evaSilh = eva.CriterionValues;
        end

        % 3. Davies-Bouldin => for euclidean distance only
        % [1] Davies, D. L., and D. W. Bouldin. ?A Cluster Separation Measure.? IEEE Transactions on Pattern Analysis and Machine Intelligence. Vol. PAMI-1, No. 2, 1979, pp. 224?227.

        if strcmp(distance,'euclidean')

            eva = evalclusters(Kmeans_results_single{1}.IDX',clust','DaviesBouldin');
            evaDB(:,it) = eva.CriterionValues;
        else
            evaDB = [];
        end  
        
         % save Kmeans results
        clear C D s t seed t_all idx_sort SUMD IDX
        save([Directory '/V1_LR_aal90_filtered_100unrelated_filterBand/LEiDA' num2str(N_areas) '_' filterBand{parse} '_results.mat'], 'Kmeans_results_single', 'evaSilh', 'N_areas', 'n_Subjects','TR') 

 elseif strcmp(mode,'louvain') 
     
        btchVec = 1:12;
        btch = btchVec(parse);
        tmpOn = ((btch-1)*9600+1);
        tmpOff = (btch*9600);
        V1_batch = V1_all(tmpOn:tmpOff,:);
        sim_cosine_matrix = pdist(V1_batch,'cosine');
        sim_mat = squareform(sim_cosine_matrix);
        clear sim_cosine_matrix V1_all
        sim_mat = ones(size(sim_mat)) - sim_mat;
        sim_mat = sim_mat - eye(size(sim_mat));
        % figure,imagesc(sim_mat)
        %% STEP1: LOUVAIN COMMUNITY DETECTION ALGORITHM

        %  gamma>1 = detects smaller modules
        %  0<=gamma<1 = detects larger modules
        %  gamma=1 = classic modularity (default)
        %  gamma = 1; % resolution parameter needs to be chosen appropriately

        numInit = 1;
        p = 1;
        mat_CIU = zeros(length(sim_mat),numInit);
        mat_Q = zeros(4, numInit);
        partitionMatrix = zeros(length(sim_mat), numInit, 4); % numofdim x numofinit x numofgamma
        gammaVec = [1:0.02:1.06];
        for gamma = gammaVec

            % creating community

            for i=1:numInit
                i
                [indComm,Q] = community_louvain(sim_mat,gamma,[],'negative_asym'); % get community assignments
                mat_Q(p,i) = Q;
                partitionMatrix(:,i,p) = indComm;
            end
            p = p + 1
        end

    %     for j=1:4
    %         D = agreement(partitionMatrix(:,1:3,j));
    %         CIU = consensus_und(D,0,3);
    %         mat_CIU(:,j) = CIU;
    %     end

        % figure; plot(mat_Q(1:4,1)); hold on; plot(mat_Q(1:4,2)); hold on; plot(mat_Q(1:4,3));

        % figure
        % [X,Y,INDSORT] = grid_communities(partitionMatrix(:,3,4)); % call function
        % imagesc(sim_mat(INDSORT,INDSORT));           % plot ordered adjacency matrix
        % hold on;                                 % hold on to overlay community visualization
        % plot(X,Y,'r','linewidth',2);
        %%
        for m = 1:4
            for l = 1:length(unique(partitionMatrix(:,1,m)))
                centroidsLouvain{m}(:,l) = sum(V1_batch(partitionMatrix(:,1,m) == l,:),1)./(sum(sum(V1_batch(partitionMatrix(:,1,m) == l,:),1)));
            end
        disp('Next Gamma')

        end  
        disp('Saving')
        save([Directory '/LEiDA_HCP/V1_aal90_Louvain/LEiDA' num2str(N_areas) '_batch' num2str(parse) 'results.mat'], 'centroidsLouvain', 'N_areas', 'V1_batch','TR','partitionMatrix','gammaVec') 


 end