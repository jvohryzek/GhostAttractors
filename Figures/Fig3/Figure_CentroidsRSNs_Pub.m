% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)
% %
% % Script to see Pyramid of LEiDA centroids colored according to best fit to YEO RSNs
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
%% LOAD DIRECTORY

Directory='/scratch1/MINDLAB2012_21-Olfaction-MEG/HCP/';

addpath(genpath(Directory))
Parcellation='AAL116';
N_areas=90; %max(Volume(:));
Extension='_filtered_100unrelated'; %'_unfiltered_100unrelated';
CortexDirection = 'TopView' %SideView

[cc_V_yeo7,p_V_yeo7] = Overlap_LEiDA_Yeo (Parcellation,N_areas,Extension,0);

%% Load the Cluster Centroids for all K
load([Directory 'LEiDA_HCP/Centroids/LEiDA' num2str(N_areas) '_Centroids_V1' Extension],'Centroids','rangeK')

YeoColor = [120 18 134; 70 30 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78]./256;
%color_vecs=[0.8 0.8 0.8; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0; 1 .5 0];

% Plot P-values
figure
for net=1:7    
        semilogy(rangeK,squeeze(p_V_yeo7(:,:,net)),'*','Color',YeoColor(net,:));
        hold on
end     
semilogy(rangeK(1)-1:rangeK(end)+1,0.05/numel(p_V_yeo7)*ones(1,length(rangeK)+2),'k','LineWidth',1)       
title('Overlap with Yeo RSNs')
ylabel('p-value')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
ylim([0 1])
box off

% LEiDA networks colored according to closest RSN
Volume=struct2array(load('ParcelsMNI2mm',['V_' Parcellation])); 
Brain_Mask_tmp=load_nifti('MNI152_T1_2mm_brain_mask.nii');
Brain_Mask=Brain_Mask_tmp.vol;
scortex=smooth3(Brain_Mask>0);

figure
for k=1:length(rangeK)
    
    VLeida=Centroids{k}.C;
    
    disp(['K= ' num2str(rangeK(k))])
    
    for Centroid=1:rangeK(k)
        
        V=VLeida(Centroid,:);
%         V=V-min(V);
%         V=2*(V/max(V));
%         V=V-1;
%        V=round(V*1e3)/1e3;  

        [cc_net, net]= max(cc_V_yeo7(k,Centroid,:));    
        %[p_net, net]= min(p_V_yeo7(k,Centroid,:)); 
        
        subplot_tight(rangeK(end),rangeK(end),Centroid+(k-1)*rangeK(end),0.01)

        hold on
        % First plot a transparent cortex
        cortexpatch=patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
        reducepatch(cortexpatch,0.01);
        isonormals(scortex,cortexpatch);
        
        Centroid_Vol=zeros(size(Volume));
        
        % To color all areas above zero
        for n=find(V>0)
            Centroid_Vol(Volume==n)=1;
        end
        
        % To color only the strongest areas
%         [~, ind_V]=sort(V,'descend');
%         for n=ind_V(1:40)
%             Centroid_Vol(Volume==n)=1;
%         end        
           
        sregion=smooth3(Centroid_Vol>0);
        %patch(isosurface(sregion,0.3), 'FaceColor', cmap(find(ind_color-V(n)>0,1),:), 'EdgeColor', 'none')%'FaceAlpha', abs(V(n)))
        if p_V_yeo7(k,Centroid,net)<0.05/k
            patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %'FaceAlpha', V(n))
        else
            % Color in Black if no significant overlap with any Yeo RSN
            patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %'FaceAlpha', V(n))
        end
        
        material dull
        lighting gouraud
        switch CortexDirection
            case 'TopView'
                % Top view
                view(-90,90)
            case 'SideView'
                % Side view
                view(0,0)
        end
        daspect([1 1 1])
        camlight;
        xlim([5 105])
        ylim([5 85])
        zlim([10 80])
        axis off
        
        if p_V_yeo7(k,Centroid,net)<0.05/k
            p_label=num2str(p_V_yeo7(k,Centroid,net)*10,'%1.0e');
            title(['r=' num2str(cc_V_yeo7(k,Centroid,net),2) ' p<1' p_label(2:end) ''],'Fontsize',6)
        end

    end
end
%% Saving
saveas(gcf,[Directory 'Figures/Fig3/Figure3_Repertoire_FNs' Extension CortexDirection],'jpg')
saveas(gcf,[Directory 'Figures/Fig3/Figure3_Repertoire_FNs' Extension CortexDirection],'fig')


