function Cluster_RL_into_LR_Centroids_Pub %(mode,parcellation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)
%
% This function projects the RL V1 timepoints onto the LR centrodis
%
%  Choose mode:
%  'parse'            = > computes the RL centroids projected onto the LR session
%  'parse_filterBand' = > computes a sequence of filtrations of RL and
%                         their appropriate LR centroids (so far works for 4)
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

mode='parse_filterBand'; % 'parse'
parcellation='aal90'; 
filtering = 'filtered'; % 'unfiltered'
encoding = 'RL'; % 'LR'
filterBand = '0.01-0.07'; % '0.01-0.1'; '0.01-0.2';'0.01-0.6'

% Code to cluster the eigenvectors from the Right-Left session using the
% centroids from the Left-Right session.
if strcmp(mode,'parse') % compute LEiDA


%% Directory
% First load the Centroids from the Left-Right session
Directory='/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/LEiDA_HCP/';

addpath(genpath(Directory))
Parcellation='AAL116';
N_areas=90; %max(Volume(:));
Extension='_unfiltered_100unrelated'; % '_filtered_100unrelated
K=5;
load([Directory 'Centroids/LEiDA' num2str(N_areas) '_Centroids_V1' Extension])

Centroids_K5_LR = struct2array(Centroids{rangeK==K});
P_K5_LR=squeeze(P(:,rangeK==K,1:K));
LT_K5_LR=squeeze(LT(:,rangeK==K,1:K));

% kmeans results
Kmeans_results_LR = cell(1);
Kmeans_results_LR{1}=Kmeans_results{rangeK==K};

for sub=1:99
    for k = 1:5
        for l = 1:5
            PTnorm_K5_LR(sub,k,l) = PTnorm{sub,K-1}(k,l);
            PT_K5_LR(sub,k,l) = PT{sub,K-1}(k,l);
        end
    end
end
% Then load the eigenvectors from the Right-Left
load([Directory 'Eigenvectors/LEiDAaal90_RL' Extension '.mat'])

%%
% Now run only 1 iteration of kmeans setting the centroids as start vectors
Kmeans_results_RL=cell(1);

disp(['- ' num2str(K) ' FC states'])
[IDX, C, SUMD, D]=kmeans(V1_all,K,'Distance','Cosine','MaxIter',1,'Display','off','Start',Centroids_K5_LR);

Kmeans_results_RL{1}.IDX=IDX';  % Cluster time course - numeric collumn vectors
Kmeans_results_RL{1}.C=C;       % Cluster centroids (FC patterns)
Kmeans_results_RL{1}.SUMD=SUMD; % Within-cluster sums of point-to-centroid distances
Kmeans_results_RL{1}.D=D;

% Calculate Probabilities and Lifetimes
P_K5_RL=zeros(n_Subjects,K);
LT_K5_RL=zeros(n_Subjects,K);
PT_K5_RL=zeros(n_Subjects,K,K);
PTnorm_K5_RL=zeros(n_Subjects,K,K);

for s=1:n_Subjects
    
    % Select the time points representing this subject and task
    Ctime=Kmeans_results_RL{1}.IDX(Time_all==s);
    
    for c=1:K
        % Probability
        P_K5_RL(s,c)=mean(Ctime==c);
        
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
        LT_K5_RL(s,c)=mean(C_Durations)*TR;
    end
    % Probability Transition Matrix
    transferMatrix_K5_RL = zeros(K,K);
    for tp = 2:size(Ctime,2)
        transferMatrix_K5_RL(Ctime(tp-1),Ctime(tp)) = transferMatrix_K5_RL(Ctime(tp-1),Ctime(tp)) + 1;
    end
    PT_K5_RL(s,:,:) = transferMatrix_K5_RL./(size(Ctime,2)-1); % normalised by T-1
    PTnorm_K5_RL(s,:,:) = squeeze(PT_K5_RL(s,:,:))./P_K5_RL(s,:)'; % normalising the rows by the probabiliy of occurence of each state

end
%%
% Finally to save all in one file

save([Directory '/K5_results_LR_RL' Extension '.mat'],...
   'Centroids_K5_LR', 'Kmeans_results_LR', 'Kmeans_results_RL', ...
   'LT_K5_LR', 'LT_K5_RL', 'P_K5_LR', 'P_K5_RL', 'PT_K5_LR', 'PT_K5_RL', 'PTnorm_K5_LR', 'PTnorm_K5_RL')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmp(mode,'parse_filterBand') % compute LEiDA
%% Directory
% First load the Centroids from the Left-Right session
Directory='/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/LEiDA_HCP/';
addpath(genpath(Directory))
N_areas=90;
Extension='_filterBand_100unrelated'; % '_filtered_100unrelated
K=5;

% loading the filterBands for all the filtering solutions
load(strcat([Directory,'Centroids/LEiDA90_Centroids_V1_filtered_parse_filterBand_100unrelated.mat']))

for i=1:4
    Centroids_K5_LR{i} = struct2array(Centroids{i});
    Kmeans_results_LR{i} = Kmeans_results{i};
end
% Centroids_K5_LR=Centroids{rangeK==K};
P_K5_LR=P;
LT_K5_LR=LT;
PTnorm_K5_LR = PTnorm;
PT_K5_LR = PT;
Kmeans_results_RL=cell(4);
% Calculate Probabilities and Lifetimes
n_Subjects = 99;
P_K5_RL=zeros(n_Subjects,4,5);
LT_K5_RL=zeros(n_Subjects,4,5);
PT_K5_RL=cell(n_Subjects,4);
PTnorm_K5_RL=cell(n_Subjects,4);

filterBand = {'0.01-0.07', '0.01-0.1', '0.01-0.2', '0.01-0.6'};
for filt=1:4
    load(strcat(Directory, 'Eigenvectors/LEiDAaal90_RL_filtered_', filterBand{filt}, '_100unrelated.mat'))

    % Now run only 1 iteration of kmeans setting the centroids as start vectors

    disp(['- ' filterBand{filt} ' Hz'])
    [IDX, C, SUMD, D]=kmeans(V1_all,5,'Distance','Cosine','MaxIter',1,'Display','off','Start',Centroids_K5_LR{filt});

    Kmeans_results_RL{filt}.IDX=IDX';   % Cluster time course - numeric column vectors
    Kmeans_results_RL{filt}.C=C;       % Cluster centroids (FC patterns)
    Kmeans_results_RL{filt}.SUMD=SUMD; % Within-cluster sums of point-to-centroid distances
    Kmeans_results_RL{filt}.D=D;

    for s=1:n_Subjects

        % Select the time points representing this subject and task
        Ctime=Kmeans_results_RL{filt}.IDX(Time_all==s);

        for c=1:K
            % Probability
            P_K5_RL(s,filt,c)=mean(Ctime==c);

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
            LT_K5_RL(s,filt,c)=mean(C_Durations)*TR;
        end
        % Probability Transition Matrix
        transferMatrix_K5_RL = zeros(K,K);
        for tp = 2:size(Ctime,2)
            transferMatrix_K5_RL(Ctime(tp-1),Ctime(tp)) = transferMatrix_K5_RL(Ctime(tp-1),Ctime(tp)) + 1;
        end
        PT_K5_RL{s,filt} = transferMatrix_K5_RL./(size(Ctime,2)-1); % normalised by T-1
        PTnorm_K5_RL{s,filt} = squeeze(PT_K5_RL{s,filt})./squeeze(P_K5_RL(s,filt,:)); % normalising the rows by the probabiliy of occurence of each state

    end
end
%%
% Finally to save all in one file
save([Directory '/K5_results_LR_RL' Extension '.mat'],...
    'Centroids_K5_LR', 'Kmeans_results_LR', 'Kmeans_results_RL', ...
    'LT_K5_LR', 'LT_K5_RL', 'P_K5_LR', 'P_K5_RL', 'PT_K5_LR', 'PT_K5_RL', 'PTnorm_K5_LR', 'PTnorm_K5_RL')

end
clear all