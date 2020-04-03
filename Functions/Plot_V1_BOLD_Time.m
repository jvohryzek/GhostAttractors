function Plot_V1_BOLD_Time
%%%%%%%%%%%%%
%
%  Function to plot V1 over time on the cortex
%  The assigned cluster is reported on top 
%
%%%%%%%%%%%%%

% Aarhus Server Directory and folder
% Directory='/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/LEiDA_HCP/'
% Folder='V1_aal90_contCORRECTEDfiltered';
% Directory = [Directory Folder];

% Joana local directory
Directory='/Users/joana/Documents/Work/LEiDA general/LEiDA_HCP/';
addpath(genpath(Directory))
K=6;
N_areas=90;
Extensions = '_contCORRECTEDfiltered'; % [],'_filtered','_contCORRECTED','_contCORRECTEDfiltered','_RL_filtered'
Parcellation='aal90';
ParcelMNIname='AAL116';
BOLD_file='hcpunrelated95_rest_aal116_LR.mat';

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
Time_interval=Best_Time-10:Best_Time+10;
Vecs_to_plot=V1_all(Time_interval,:);
IDX_Vecs=Kmeans_results_single{1}.IDX(Time_interval);
clear V1_all
clear Kmeans_results_single

disp(['Closest observation of centroid ' num2str(State) ' detected in subject ' num2str(ceil(Time_interval(1)/Tsubject))])

% Create warning if the time points capture different subjects
if ceil(Time_interval(1)/Tsubject)~=ceil(Time_interval(end)/Tsubject)
    disp('Attention: This time interval contains frames from different 2 subjects')
end

% Load corresponding BOLD signal for this subject in this Time interval
BOLD=struct2array(load(BOLD_file));
BOLD_subject=BOLD{ceil(Time_interval(1)/Tsubject)}.tps(1:N_areas,:);
clear BOLD

% Remove mean from BOLD, filter and compute the Hilbert Phase

% FILTER SETTINGS
TR=0.78;
fnq = 1/(2*TR);                 % Nyquist frequency
flp = 0.04;                     % lowpass frequency of filter
fhi = 0.07;                     % highpass
Wn = [flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k = 2;                          % 2nd order butterworth filter
[bfilt,afilt] = butter(k,Wn);   % construct the filter


BOLD_Phase=zeros(size(BOLD_subject));

for seed=1:N_areas
    BOLD_subject(seed,:)=BOLD_subject(seed,:)-mean(BOLD_subject(seed,:));
    BOLD_subject(seed,:) = filtfilt(bfilt,afilt,BOLD_subject(seed,:));
    BOLD_Phase(seed,:) = angle(hilbert(BOLD_subject(seed,:)));
end

% Save the relevant interval only and clear the rest 
Time_BOLDSubject=mod(Time_interval,Tsubject)+1; % We need to add 1 because we remove first frame for LEiDA
BOLD_interval=BOLD_subject(:,Time_BOLDSubject);
Phase_interval=BOLD_Phase(:,Time_BOLDSubject);

disp('Generating Figure')

figure('Name',[Parcellation Extensions ' centred in centroid ' num2str(State) ' in subject ' num2str(ceil(Time_interval(1)/Tsubject))],'color','white')

subplot_tight(5,1,1)
hold on
lim=max(abs(BOLD_subject(:)));
% plot([Time_BOLDSubject(1) Time_BOLDSubject(1)],[-lim lim],':k');
% plot([Time_BOLDSubject(end) Time_BOLDSubject(end)],[-lim lim],':k');
x = [Time_BOLDSubject(1) Time_BOLDSubject(end) Time_BOLDSubject(end) Time_BOLDSubject(1)];
y = [-lim -lim lim lim];
patch(x,y,'k','LineStyle','none','FaceAlpha',0.1);
plot(BOLD_subject')
xlabel ('Time (TR)')
ylabel('BOLD change')
ylim([-lim lim])
title(['A - Mean BOLD signals in 90 brain areas band-pass filtered [' num2str(flp) '-' num2str(fhi) 'Hz] in subject ' num2str(ceil(Time_interval(1)/Tsubject))])
clear BOLD_subject BOLD_Phase

%%

subplot_tight(5,1,2)
hold on
plot(Time_BOLDSubject,BOLD_interval)
xlabel ('Time (TR)')
ylabel('BOLD change')
lim=max(abs(BOLD_interval(:)));
ylim([-lim lim])
xlim([Time_BOLDSubject(1) Time_BOLDSubject(end)])
for t=1:length(Time_interval)
    V=Vecs_to_plot(t,:);
    hold on
    if sum(V<0)
        plot(Time_BOLDSubject(t),BOLD_interval(V<0,t),'ob')
    end
    if sum(V>=0)
        plot(Time_BOLDSubject(t),BOLD_interval(V>=0,t),'.r')
    end
end
box off
set(gca,'XTick',Time_BOLDSubject)
set(gca,'XGrid','on')
title(['B - Zoom of A into an interval of ' num2str(numel(Time_BOLDSubject)) ' frames'])


subplot_tight(5,1,3)
hold on
R0=zeros(N_areas,1);
for t=1:length(Time_interval)
    plot3([Time_BOLDSubject(t) Time_BOLDSubject(t)],[R0(Vecs_to_plot(t,:)>=0) sin(Phase_interval(Vecs_to_plot(t,:)>=0,t))]',[R0(Vecs_to_plot(t,:)>=0) cos(Phase_interval(Vecs_to_plot(t,:)>=0,t))]','r')
    plot3([Time_BOLDSubject(t) Time_BOLDSubject(t)],[R0(Vecs_to_plot(t,:)<0) sin(Phase_interval(Vecs_to_plot(t,:)<0,t))]',[R0(Vecs_to_plot(t,:)<0) cos(Phase_interval(Vecs_to_plot(t,:)<0,t))]','b')
end
box off
ylabel ('Imaginary')
zlabel('Real axis')
ylim([-1 1])
zlim([-1 1])
xlim([Time_BOLDSubject(1) Time_BOLDSubject(end)])
view(2,17)
grid on
set(gca,'XTick',Time_BOLDSubject)
%     [Vec,Val]=eigs(cos(Phase_interval(:,t)-Phase_interval(:,t)'),2);
%     Phase_V1_mean=mean(acos(cos(angle(Val(1,1)*Vec(:,1)+1i*(Val(2,2)*Vec(:,2)))-Phase_interval(:,t))));    
%     plot([-2*sin(Phase_V1_mean) 2*sin(Phase_V1_mean)]',[-2*cos(Phase_V1_mean) 2*cos(Phase_V1_mean)]','k','LineWidth',2)
%     plot([-2*sin(Phase_V1_mean+pi/2) 2*sin(Phase_V1_mean+pi/2)]',[-2*cos(Phase_V1_mean+pi/2) 2*cos(Phase_V1_mean+pi/2)]','k','LineWidth',2)
title('C - BOLD Phase portrait of the 90 Brain areas over time')

Volume=struct2array(load('ParcelsMNI2mm',['V_' ParcelMNIname]));
% Create a colormap with one color per state
Order=[1:2:N_areas N_areas:-2:2];

for t=1:length(Time_interval)
    
    V=Vecs_to_plot(t,:);
    
    Centroid_Vol=zeros(size(Volume));  
    % To color all areas above zero
    for n=1:N_areas
        if V(n)>=0
            Centroid_Vol(Volume==n)=1;
        elseif V(n)<0
            Centroid_Vol(Volume==n)=-1;
        end
    end
    
    subplot_tight(5,length(Time_interval)+2,t+3*(length(Time_interval)+2),0.01)
    hold on
    sregion=smooth3(Centroid_Vol>0);
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
    
    if t==1
        title('D - Brain areas colored according to projection into BOLD Phase Leading Eigenvector')
        axis on
    end
    
    %title(IDX_Vecs(t),'Color',color_vecs(IDX_Vecs(t),:),'Fontsize',16,'FontWeight','bold')
        
    subplot_tight(5,length(Time_interval)+2,t+4*(length(Time_interval)+2),0.01)
    colormap(jet)
    iFC=cos(Phase_interval(:,t)-Phase_interval(:,t)');
    imagesc(iFC(Order,Order),[-1 1])
    axis square
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    
    title(num2str(Time_BOLDSubject(t)))
    
   if t==1
        title('E - Instantaneous matrices of BOLD Phase synchronization')
        xlabel('90 Brain areas')
        ylabel('90 Brain areas')
    end
    
end

