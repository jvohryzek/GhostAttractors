function LEiDA_HCP_standalone_Pub

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)
%
% This function processes, clusters and analyses BOLD data using LEiDA
%
% Adpated for the analysis of Healthy data from Human Connectome Project
% on different parcellations (aal116, aal90, Glasser378, Glasser360, dbs80)
%
%  Choose mode:
%  'run'              = > calculates and save the V1_all that is
%                         subsequently run on a cluster by
%                         LEiDA_HCP_parsing.m (mode 'kmeans' and 'kmeans_filtrationBand')
%  'parse'            = > computes and saves all the LEiDA metric such
%                         as P, LT and PTR
%  'parse_filterBand' = > computes and saves all the LEiDA metric for
%                         different frequency bands
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
mode='parse'; % 'run' 'parse_filterBand''
parcellation='aal90'; %'glasser378','dbs80','glasser360','aal116','lausanne2008_scale1','lausanne2008_scale2','lausanne2008_scale3','lausanne2008_scale4','lausanne2008_scale5'
filtering = 'unfiltered'; % 'filtered'
encoding = 'LR'; % 'RL'
filterBand = '0.01-0.6'; % '0.01-0.07'; '0.01-0.1';'0.01-0.2'
Directory  = '/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/';

% Select folder where data is stored

if strcmp(mode,'run') % compute LEiDA
    
    %% 1 - Compute the Leading Eigenvectors from the BOLD datasets
    
    disp('Processing the eigenvectors from BOLD data')
    
    % Define here data parameters:
    
    %N_areas    = 90;
    TR         = 0.72; % based on Shine et al. 2016 and HCP reference manual
    Tmax       = 1200;  % Set here the total number of frames in each scan
    % FILTER SETTINGS => to be checked
    fnq = 1/(2 * TR);               % Nyquist frequency
    switch filterBand
        case '0.01-0.07'
            flp = 0.01; %0.04;       % lowpass frequency of filter
            fhi = 0.07;            % highpass
        case '0.01-0.1'
            flp = 0.01; %0.04;       % lowpass frequency of filter
            fhi = 0.1;             % highpass
        case '0.01-0.2'
            flp = 0.01; %0.04;       % lowpass frequency of filter
            fhi = 0.2              % highpass
        case '0.01-0.6'
            flp = 0.01; %0.04;       % lowpass frequency of filter
            fhi = fnq-0.01              % highpass
    end
    Wn = [flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
    k = 2;                          % 2nd order butterworth filter
    [bfilt,afilt] = butter(k,Wn);   % construct the filter

    % Note: if Tmax is different at each scan let me know to adapt.
    
    if strcmp(parcellation,'glasser378')
        % GLASSER
        file='hcpunrelated100_rest_glasser378.mat';
        load([Directory file],'s_glasser378');
        N_areas=378;
        n_Subjects = length(s_glasser378); % Total number of scans
        % I removed fields to reduce file size to half
        for s=1:n_Subjects
            s_glasser378{1,s}=rmfield(s_glasser378{1,s},'glasserts_stderr');
            s_glasser378{1,s}=rmfield(s_glasser378{1,s},'fc');
        end
    elseif strcmp(parcellation,'dbs80')
            file='hcpunrelated100_REST_dbs80.mat';
            load([Directory file],'subject');
            N_areas=80;
            n_Subjects = length(subject); % Total number of scans
            % I removed fields to reduce file size to half
            for s=1:n_Subjects
                subject{1,s}=rmfield(subject{1,s},'dbs80ts_stderr');
                subject{1,s}=rmfield(subject{1,s},'fc_symm');
            end
        elseif strcmp(parcellation,'glasser360')
            file='hcpunrelated100_resting_glasser.mat';
            load([Directory file],'subject');
            N_areas=360;
            n_Subjects = length(subject); % Total number of scans
            % I removed fields to reduce file size to half
            for s=1:n_Subjects
                subject{1,s}=rmfield(subject{1,s},'glasserts_stderr');
                subject{1,s}=rmfield(subject{1,s},'fc');
            end 
        elseif strcmp(parcellation,'aal116')
            file='hcpunrelated100_rest_aal116_LR.mat';
            load([Directory file],'aal_LRts');
            N_areas=116;
            n_Subjects = length(aal_LRts); % Total number of scans
        elseif strcmp([parcellation encoding],'aal90LR')
            file='hcpunrelated100_rest_aal116_LR.mat';
            load([Directory file],'aal_LRts');
            N_areas=90;
            n_Subjects = length(aal_LRts); % Total number of scans
        elseif strcmp([parcellation encoding],'aal90RL')
            file='hcpunrelated100_rest_aal116_RL.mat';
            load([Directory file],'aal_RLts');
            N_areas=90;
            n_Subjects = length(aal_RLts); % Total number of scans
        elseif strcmp(parcellation,'lausanne2008_scale1')
            file='hcpunrelated100_rest_lausanne2008_scale1_LR.mat';
            load([Directory file],'lausanne2008_scale1_LRts');
            N_areas=83;
            n_Subjects = length(lausanne2008_scale1_LRts); % Total number of scans
        elseif strcmp(parcellation,'lausanne2008_scale2')
            file='hcpunrelated100_rest_lausanne2008_scale2_LR.mat';
            load([Directory file],'lausanne2008_scale2_LRts');
            N_areas=129;
            n_Subjects = length(lausanne2008_scale2_LRts); % Total number of scans
        elseif strcmp(parcellation,'lausanne2008_scale3')
            file='hcpunrelated100_rest_lausanne2008_scale3_LR.mat';
            load([Directory file],'lausanne2008_scale3_LRts');
            N_areas=234;
            n_Subjects = length(lausanne2008_scale3_LRts); % Total number of scans
        elseif strcmp(parcellation,'lausanne2008_scale4')
            file='hcpunrelated100_rest_lausanne2008_scale4_LR.mat';
            load([Directory file],'lausanne2008_scale4_LRts');
            N_areas=463;
            n_Subjects = length(lausanne2008_scale4_LRts); % Total number of scans
        elseif strcmp(parcellation,'lausanne2008_scale5')
            file='hcpunrelated100_rest_lausanne2008_scale5_LR.mat';
            load([Directory file],'lausanne2008_scale5_LRts');
            N_areas=1015;
            n_Subjects = length(lausanne2008_scale5_LRts); % Total number of scans
    end
    
    % Allocate space for Eigenvectors and associated information
    V1_all = zeros((Tmax-2)*n_Subjects,N_areas); % All leading eigenvectors
    t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*(Tmax-2))
    Time_all= zeros((Tmax-2)*n_Subjects,1); % Vector that links each frame to a subject
    
    for s=1:n_Subjects
        
        if strcmp(parcellation,'glasser378')
            BOLD=s_glasser378{1,s}.glasserts;
        elseif strcmp(parcellation,'glasser360')
            BOLD=subject{1,s}.glasserts;
        elseif strcmp(parcellation,'dbs80')
            BOLD=subject{1,s}.dbs80ts;
        elseif strcmp(parcellation,'aal116')
            BOLD=aal_LRts{1,s}.tps;
        elseif strcmp([parcellation encoding],'aal90LR')
            BOLD=aal_LRts{1,s}.tps;
        elseif strcmp([parcellation encoding],'aal90RL')
            BOLD=aal_RLts{1,s}.tps;
        elseif strcmp(parcellation,'lausanne2008_scale1')
            BOLD=lausanne2008_scale1_LRts{1,s}.tps;
        elseif strcmp(parcellation,'lausanne2008_scale2')
            BOLD=lausanne2008_scale2_LRts{1,s}.tps;
        elseif strcmp(parcellation,'lausanne2008_scale3')
            BOLD=lausanne2008_scale3_LRts{1,s}.tps;
        elseif strcmp(parcellation,'lausanne2008_scale4')
            BOLD=lausanne2008_scale4_LRts{1,s}.tps;
        elseif strcmp(parcellation,'lausanne2008_scale5')
            BOLD=lausanne2008_scale5_LRts{1,s}.tps;
        end
                
        Phase_BOLD=zeros(N_areas,Tmax);

        switch filtering
        % OPTION: filtered - Get the BOLD phase using the Hilbert transform
            case 'filtered'
                for seed=1:N_areas
                    BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
                    BOLD_processed = filtfilt(bfilt,afilt,BOLD(seed,:));            
                    Phase_BOLD(seed,:) = angle(hilbert(BOLD_processed));
                end
        % OPTION: unfiltered - Get the BOLD phase using the Hilbert transform
            case 'unfiltered'
                for seed=1:N_areas
                    BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
                    Phase_BOLD(seed,:) = angle(hilbert(BOLD(seed,:)));
                end
        end
        
       [V1previous,~]=eigs(cos(Phase_BOLD(:,1)-Phase_BOLD(:,1)'),1);

        % Slide over time discarding the first and last epochs
        for t=2:Tmax-1
            
            % Get the leading eigenvector of Instantaneous BOLD Phases
            [V1,~]=eigs(cos(Phase_BOLD(:,t)-Phase_BOLD(:,t)'),1);
            
            % OPTION 1: half-switching
            if sum(V1)>0
                V1=-V1;
            end
            
            % OPTION 2: continous-switching
            
            distV1vsV1previous = 1-pdist([V1,V1previous]','cosine');
            distmirrorV1vsV1previous = 1-pdist([-V1,V1previous]','cosine');
            if distV1vsV1previous > distmirrorV1vsV1previous
                V1previous = V1;
            elseif distmirrorV1vsV1previous > distV1vsV1previous
                V1previous = -V1;
            end
                
            % Save V1 from all frames in all fMRI sessions
            t_all = t_all+1; % Update time
            Time_all(t_all) = s;

            % OPTION 1: half-switching
            V1_all(t_all,:) = V1;
            
            % OPTION 2: continous-switching
            V1cont_all(t_all,:) = V1previous;
        end
    end
    clear BOLD iFC V1 Phase_BOLD V1previous
    save([Directory 'LEiDA_HCP/Eigenvectors/LEiDA' parcellation '_' encoding '_' filtering '_' filterBand '_100unrelated.mat'],'TR', 'N_areas','n_Subjects', 'V1_all', 'V1cont_all','Time_all')

    disp('%%%%% LEiDA --run-- SUCCESSFULLY COMPLETED %%%%%%%')
    disp('LEiDA results are now saved.')

elseif strcmp(mode,'parse') 
    % Load kmeans cluster results and compute clustering results analysis

    mink=2;
    maxk=20;
    rangeK=mink:maxk;
    Kmeans_results={};
    
    if strcmp([parcellation encoding],'aal90LR')
        % loading the 'Time_all' from /Eigenvectors/LEiDAaal90_LR_filtered_100unrelated.mat'
        load([Directory 'LEiDA_HCP/Eigenvectors/LEiDAaal90_LR_' filtering '_100unrelated.mat'],'Time_all');
        N_areas=90;
        for K=mink:maxk
            load([Directory 'LEiDA_HCP/V1_LR_aal90_' filtering '_100unrelated/LEiDA' num2str(N_areas) '_K' num2str(K) 'results.mat']) %V1_all Time_all Kmeans_results P LT P_pval LT_pval rangeK n_Controls n_Patients
            Kmeans_results = [Kmeans_results; Kmeans_results_single];
            evaSilh_clusters(K) = evaSilh;
        end
    elseif strcmp([parcellation encoding],'aal90RL') % usually calculate this from centorids of LR
        load([Directory 'LEiDA_HCP/Eigenvectors/LEiDAaal90_RL_V1.mat'],'Time_all');
        N_areas=90;
        for K=mink:maxk
            load([Directory 'LEiDA_HCP/V1_RL_aal90/LEiDA' num2str(N_areas) '_K' num2str(K) 'results.mat']) %V1_all Time_all Kmeans_results P LT P_pval LT_pval rangeK n_Controls n_Patients
            Kmeans_results = [Kmeans_results; Kmeans_results_single];
            evaSilh_clusters(K) = evaSilh;
        end
    end
    %% Plotting the quality measure for individiual clustering solutions
    fSI1 = figure,plot(2:20,evaSilh_clusters(2:20),'k','LineWidth',4)
    ylabel('Silhouette Value'); xlabel('# of clusters');
    set(gca, 'FontSize',18)
    
    %%
    % For every fMRI scan calculate probability and lifetimes of each state c.
    P=zeros(n_Subjects,maxk-mink+1,maxk);
    LT=zeros(n_Subjects,maxk-mink+1,maxk);
    PT = cell(n_Subjects,size(rangeK,2));
    PTnorm = cell(n_Subjects,size(rangeK,2));

    for K=1:length(rangeK)

        for s=1:n_Subjects
            
            % Select the time points representing this subject
            T = (Time_all==s);
            Ctime=Kmeans_results{K}.IDX(T);
            
            for c=1:rangeK(K)
                % Probability
                P(s,K,c)=mean(Ctime==c) + eps; % normalised by T
                %Punnorm(s,K,c)=sum(Ctime==c)+eps; % for the PT normalisation
                % Mean Lifetime
                Ctime_bin=Ctime==c;
                
                % Detect switches in and out of this state
                a=find(diff(Ctime_bin)==1);
                b=find(diff(Ctime_bin)==-1);
                
                % We discard the cases where state sarts or ends ON
                if length(b)>length(a)
                    b(1)=[];
                elseif length(a)>length(b)
                    a(end)=[];
                elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
                    b(1)=[];
                    a(end)=[];
                end
                if ~isempty(a) && ~isempty(b)
                    C_Durations=b-a;
                else
                    C_Durations=0;
                end
                LT(s,K,c)=mean(C_Durations)*TR;
            end
            % Probability Transition Matrix
            transferMatrix = zeros(rangeK(K),rangeK(K));
            for tp = 2:size(Ctime,2)%-5
                transferMatrix(Ctime(tp-1),Ctime(tp)) = transferMatrix(Ctime(tp-1),Ctime(tp)) + 1;
            end
            PT{s,K} = transferMatrix./(size(Ctime,2)-1); % normalised by T-1
            PTnorm{s,K} = PT{s,K}./squeeze(nonzeros(P(s,K,:)));

        end
    end
    %% SAVING FOR JOANA's FIGURES 1,2,3,4
    Centroids = Kmeans_results;

    for i =1:19
        Centroids{i,1} = rmfield(Centroids{i,1},'IDX');
        Centroids{i,1} = rmfield(Centroids{i,1},'SUMD');
        Centroids{i,1} = rmfield(Centroids{i,1},'D');

    end
    Centroids = Centroids';
    save([Directory 'LEiDA_HCP/Centroids/LEiDA90_Centroids_V1_' filtering '_100unrelated.mat'], 'rangeK','P','LT','PTnorm','PT','Centroids','Kmeans_results')
    %% PLOT first 20 sbujects' Probability Matrices
    fSI2 = figure,
    for i=1:20
        subplot(4,5,i)
        imagesc(PTnormOrder6{i+79,4})
        colorbar
        title(sprintf('Sbj%d',i))
        axis square
        hold on
    end
    
    %% histogram for all the transition 
    for i=1:n_Subjects
        tmp = PTnorm{i,4}';
        tmp2 = tmp(:);
        for j=1:25
            histPT(i,j)=tmp2(j);
        end
    end
    fSI3 = figure
    for k=1:25
        subplot(5,5,k)
        histogram(histPT(:,k),10)
        hold on
        title(sprintf('Edge %d',k))
        ylabel('# of Subjects')
        xlabel('# of tps')
        xlim([0 1])

    end
    clear Ctime_bin Ctime C D a b n p s t k c seed t_all C_Durations ind_sort idx_sort SUMD IDX
    saveas(fSI1,[Directory 'Figures/SI_Figs/Silh_quality_measure' filtering],'png')
    saveas(fSI1,[Directory 'Figures/SI_Figs/Silh_quality_measure' filtering],'fig')
    saveas(fSI1,[Directory 'Figures/SI_Figs/Silh_quality_measure' filtering],'pdf')
    saveas(fSI1,[Directory 'Figures/SI_Figs/Silh_quality_measure' filtering],'eps')
    
    saveas(fSI2,[Directory 'Figures/SI_Figs/ProbMat_first20sbj' filtering],'png')
    saveas(fSI2,[Directory 'Figures/SI_Figs/ProbMat_first20sbj' filtering],'fig')
    saveas(fSI2,[Directory 'Figures/SI_Figs/ProbMat_first20sbj' filtering],'pdf')
    saveas(fSI2,[Directory 'Figures/SI_Figs/ProbMat_first20sbj' filtering],'eps')
    
    saveas(fSI3,[Directory 'Figures/SI_Figs/ProbMat_across_sbj' filtering],'png')
    saveas(fSI3,[Directory 'Figures/SI_Figs/ProbMat_across_sbj' filtering],'fig')
    saveas(fSI3,[Directory 'Figures/SI_Figs/ProbMat_across_sbj' filtering],'pdf')
    saveas(fSI3,[Directory 'Figures/SI_Figs/ProbMat_across_sbj' filtering],'eps')

    disp('%%%%% LEiDA "PARSE" SUCCESSFULLY COMPLETED %%%%%%%')
    disp('LEiDA results are now saved.')
    
 %%
elseif strcmp(mode,'parse_filterBand')
    % Load kmeans cluster results and compute clustering results analysis
    % for filtration bands

	filterBand = {'0.01-0.07', '0.01-0.1', '0.01-0.2', '0.01-0.6'};
    mink=2;
    maxk=20;
    rangeK=mink:maxk;
    Kmeans_results={};
    
    if strcmp([parcellation encoding],'aal90LR')
        % loading the 'Time_all' from /Eigenvectors/LEiDAaal90_LR_filtered_100unrelated.mat'
        load([Directory 'LEiDA_HCP/Eigenvectors/LEiDAaal90_LR_' filtering '_100unrelated.mat'],'Time_all');
        N_areas=90;
        for b=1:4
            load([Directory 'LEiDA_HCP/V1_LR_aal90_' filtering '_100unrelated_filterBand/LEiDA90_' filterBand{b} '_results.mat']) %V1_all Time_all Kmeans_results P LT P_pval LT_pval rangeK n_Controls n_Patients
            Kmeans_results = [Kmeans_results; Kmeans_results_single];
        end
    end
    %%
    % For every fMRI scan calculate probability and lifetimes of each state c.
    P=zeros(n_Subjects,4,5);
    LT=zeros(n_Subjects,4,5);
    PT = cell(n_Subjects,4);
    PTnorm = cell(n_Subjects,4);

    for band=1:4

        for s=1:n_Subjects
            
            % Select the time points representing this subject and task
            T = (Time_all==s);
            Ctime=Kmeans_results{band}.IDX(T);
            
            for c=1:5
                % Probability
                P(s,band,c)=mean(Ctime==c) + eps; % normalised by T
                %Punnorm(s,K,c)=sum(Ctime==c)+eps; % for the PT normalisation
                % Mean Lifetime
                Ctime_bin=Ctime==c;
                
                % Detect switches in and out of this state
                a=find(diff(Ctime_bin)==1);
                b=find(diff(Ctime_bin)==-1);
                
                % We discard the cases where state sarts or ends ON
                if length(b)>length(a)
                    b(1)=[];
                elseif length(a)>length(b)
                    a(end)=[];
                elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
                    b(1)=[];
                    a(end)=[];
                end
                if ~isempty(a) && ~isempty(b)
                    C_Durations=b-a;
                else
                    C_Durations=0;
                end
                LT(s,band,c)=mean(C_Durations)*TR;
            end
            % Probability Transition Matrix
            transferMatrix = zeros(5,5);

            for tp = 2:size(Ctime,2)%-5
                transferMatrix(Ctime(tp-1),Ctime(tp)) = transferMatrix(Ctime(tp-1),Ctime(tp)) + 1;
            end
            PT{s,band} = transferMatrix./(size(Ctime,2)-1); % normalised by T-1
            PTnorm{s,band} = PT{s,band}./squeeze(nonzeros(P(s,band,:)));

        end
    end
    %% SAVING FOR JOANA's FIGURES 1,2,3,4
    Centroids = Kmeans_results;

    for i =1:4
        Centroids{i,1} = rmfield(Centroids{i,1},'IDX');
        Centroids{i,1} = rmfield(Centroids{i,1},'SUMD');
        Centroids{i,1} = rmfield(Centroids{i,1},'D');

    end
    Centroids = Centroids';
    save([Directory 'LEiDA_HCP/Centroids/LEiDA90_Centroids_V1_' filtering '_' mode '_100unrelated.mat'], 'rangeK','P','LT','PTnorm','PT','Centroids','Kmeans_results')
     clear Ctime_bin Ctime C D a b n p s t k c seed t_all C_Durations ind_sort idx_sort SUMD IDX

    disp('%%%%% LEiDA "PARSE" SUCCESSFULLY COMPLETED %%%%%%%')
    disp('LEiDA results are now saved.')
    
end
