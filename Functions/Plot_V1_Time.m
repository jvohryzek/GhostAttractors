function Plot_V1_Time
%%%%%%%%%%%%%
%
%  Function to plot V1 over time on the cortex
%  The assigned cluster is reported on top 
%
%%%%%%%%%%%%%

% Joana local directory
Directory='/Users/joana/Documents/Work/LEiDA general/LEiDA_HCP/';
addpath(genpath(Directory))

% Aarhus Server Directory and folder
% Directory='/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/LEiDA_HCP/'
% Folder='V1_aal90_contCORRECTEDfiltered';
% Directory = [Directory Folder];

K=6;
N_areas=90;
Extensions = '_contCORRECTED'; % [],'_filtered','_contCORRECTED','_contCORRECTEDfiltered','_RL_filtered'
Parcellation='aal90';
ParcelMNIname='AAL116';

% LOAD THE EIGENVECTORS 
V1_file=['LEiDA' Parcellation '_V1' Extensions(2:end) '.mat'];
if numel(Extensions) && contains(Extensions,'cont')
    V1_all=struct2array(load(V1_file,'V1previous_all'));
else
    V1_all=struct2array(load(V1_file,'V1_all'));  
end
disp('Loading the eigenvectors from:')
disp(V1_file)


kmeans_file=[Directory 'LEiDA' num2str(N_areas) '_K' num2str(K) Extensions '.mat'];
load(kmeans_file,'Kmeans_results_single')
disp(['Loading kmeans results for K='  num2str(K) ' from:'])
disp(kmeans_file)

State = input('Choose State of Interest: ');
% For example, look for the time points where the selected state
% is maximally expressed (minimal distance to the centroid):
[~, Best_Time]=min(Kmeans_results_single{1}.D(:,State));
Tsubject=1200-2;
Time_interval=Best_Time-60:Best_Time+60;

disp(['Closest observation of centroid ' num2str(State) ' detected in subject ' num2str(ceil(Time_interval(1)/Tsubject))])

% Create warning if the time points capture different subjects
if ceil(Time_interval(1)/Tsubject)~=ceil(Time_interval(1)/Tsubject)
    disp('Attention: This time interval contains frames from different 2 subjects')
end

Vecs_to_plot=V1_all(Time_interval,:);
IDX_Vecs=Kmeans_results_single{1}.IDX(Time_interval);
clear V1_all

% Create a colormap with one color per state
color_vecs=jet(K);

%%

figure('Name',[Parcellation Extensions],'color','white')
Volume=struct2array(load('ParcelsMNI2mm',['V_' ParcelMNIname])); 

for t=1:size(Vecs_to_plot,1)
        
   V=Vecs_to_plot(t,:);

   subplot_tight(6,20,t,0.01)
   hold on
        
        Centroid_Vol=zeros(size(Volume));
        
        % To color all areas above zero
        for n=1:N_areas
            if V(n)>=0
                Centroid_Vol(Volume==n)=1;
            elseif V(n)<0
                Centroid_Vol(Volume==n)=-1;
            end
        end     
           
        sregion=smooth3(Centroid_Vol>0);
        %patch(isosurface(sregion,0.3), 'FaceColor', cmap(find(ind_color-V(n)>0,1),:), 'EdgeColor', 'none')%'FaceAlpha', abs(V(n)))
        %patch(isosurface(sregion,0.3), 'FaceColor', color_vecs(IDX_Vecs(t),:), 'EdgeColor', 'none') %'FaceAlpha', V(n))
        patch(isosurface(sregion,0.3), 'FaceColor', 'r', 'EdgeColor', 'none') %,'FaceAlpha', 0.1)      
        sregion=smooth3(Centroid_Vol<0);
        patch(isosurface(sregion,0.3), 'FaceColor', 'b', 'EdgeColor', 'none') %,'FaceAlpha', 0.1)
        
        material dull
        lighting gouraud
        % Top view
        view(-90,90)
        % Side view
        % view(0,0)
        daspect([1 1 1])
        camlight;
        xlim([5 105])
        ylim([5 85])
        zlim([10 80])
        axis off
        
        title(IDX_Vecs(t),'Color',color_vecs(IDX_Vecs(t),:),'Fontsize',16,'FontWeight','bold')
end
