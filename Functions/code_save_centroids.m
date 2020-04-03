function code_save_centroids

% Code to save Centroids
% and calculate Probabilities and Lifetimes from clustering results
% Saves the necessary info to plot repertoire and PMS in local computer

% Define folder name
Parcellation='aal90'; % 'aal116' % '378' 
Options='_filtered'; %[]; % %[]; % '_contCORRECTED'; % '_contCORRECTEDfiltered'; % 
%folder=['V1_' Parcellation Options];
folder=['V1_RL_' Parcellation Options];
Directory='/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/LEiDA_HCP/';
N_areas=90;
rangeK=2:20;
mink=rangeK(1);
maxk=rangeK(end);

load([Directory 'Eigenvectors/LEiDA' Parcellation '_RL_V1' Options(2:end) '.mat'],'n_Subjects','Time_all','TR')
     
% For every fMRI scan calculate probability and lifetimes of each state c.
P=zeros(n_Subjects,maxk-mink+1,maxk);
LT=zeros(n_Subjects,maxk-mink+1,maxk);

for K=1:length(rangeK)
      load([Directory folder '/LEiDA' num2str(N_areas) '_K' num2str(rangeK(K)) 'results.mat'],'Kmeans_results_single')
     
      Centroids{K}.C=Kmeans_results_single{1}.C;
      
      % Calculate Probabilities  
      for s=1:n_Subjects
          
          % Select the time points representing this subject and task
          Ctime=Kmeans_results_single{1}.IDX(Time_all==s);
          
          for c=1:rangeK(K)
              % Probability
              P(s,K,c)=mean(Ctime==c);
              
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
      end
end

save([Directory 'Centroids/LEiDA' num2str(N_areas) '_RL_Centroids_V1' Options],'Centroids','rangeK','P','LT')