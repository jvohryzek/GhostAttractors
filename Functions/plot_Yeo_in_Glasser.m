function plot_Yeo_in_Glasser

%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function set the color and/or transparency of brain areas.
%
%  Brain areas are defined according to a given parcellation scheme.
%  One patch is created for each brain area. 
%  Color/Transparency of each brain region are scaled by input vector V.
% 
%  Vector V has the size Nx1, where N is the number of brain areas, sorted
%  in the same order as in the parcellation scheme.
%  
%  3D volumes of a number of parcellations are saved in ParcelsMNI2mm, but
%  new volumes can be easily created from the NIFTI files
%  e.g. V_dbs80 = niftiread('dbs80symm_2mm.nii.gz');
%
%  Joana Cabral
%  October 2019
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define here the name of the Parcellation scheme to load:
% glasser378 AAL116 dbs80 Yeo7 glasser360
Parcellation='glasser378';

Volume=struct2array(load('ParcelsMNI2mm',['V_' Parcellation])); 

% Then load the mask of the Yeo parcellation in MNI space 2mm
V_Yeo=struct2array(load('ParcelsMNI2mm','V_Yeo7'));

N_areas=max(Volume(:));
N_Yeo=max(V_Yeo(:));
Yeo_in_new_Parcel=zeros(N_Yeo,N_areas);

% Create 7 vectors representing the 7 Yeo RSNs in new Parcellation scheme
for n=1:N_areas
    indn=Volume==n;
    for Net=1:7
        Yeo_in_new_Parcel(Net,n)=numel(find(V_Yeo(indn)==Net))/sum(indn(:));
    end
end
clear V_Yeo V_Parcels indn

color_vecs=[0.8 0.8 0.8; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0; 1 .5 0];

figure

for net=1:size(Yeo_in_new_Parcel,1)

    subplot(2,size(Yeo_in_new_Parcel,1),net)
    V=Yeo_in_new_Parcel(net,:);
    %V=V/max(abs(V));
    %V=round(V*10)/10;
    hold on

    scortex=smooth3(Volume>0);
    cortexpatch=patch(isosurface(scortex,0.3), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.05);
    reducepatch(cortexpatch,0.1);
    isonormals(scortex,cortexpatch);

    for n=find(V>0)
            sregion=smooth3(Volume==n);
            %patch(isosurface(sregion,0.3), 'FaceColor', cmap(find(ind_color-V(n)>0,1),:), 'EdgeColor', 'none')%'FaceAlpha', abs(V(n)))       
            patch(isosurface(sregion,0.3), 'FaceColor', color_vecs(net,:), 'EdgeColor', 'none','FaceAlpha', V(n))       
    end

    material dull
    lighting gouraud
    % Top view
    view(-90,90)
    daspect([1 1 1])
    axis on
    camlight;
    xlim([5 105])
    ylim([5 85])
    zlim([10 80])
    axis off

    subplot(2,size(Yeo_in_new_Parcel,1),net+size(Yeo_in_new_Parcel,1))
    hold on

    scortex=smooth3(Volume>0);
    cortexpatch=patch(isosurface(scortex,0.3), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.05);
    reducepatch(cortexpatch,0.1);
    isonormals(scortex,cortexpatch);

    for n=find(V>0)
            sregion=smooth3(Volume==n);
            %patch(isosurface(sregion,0.3), 'FaceColor', cmap(find(ind_color-V(n)>0,1),:), 'EdgeColor', 'none')%'FaceAlpha', abs(V(n)))       
            patch(isosurface(sregion,0.3), 'FaceColor', color_vecs(net,:), 'EdgeColor', 'none','FaceAlpha', V(n))       
    end

    material dull
    lighting gouraud
    % Top view
    view(0,0)
    daspect([1 1 1])
    axis on
    camlight;
    xlim([5 105])
    ylim([5 85])
    zlim([10 80])
    axis off
    
end




% 
% % Image of the left hemisphere
% subplot_tight(2,2,3,0.0001)
% hold on
% for n=find(V>=0)
%     if find(181:369==n)
%         sregion=smooth3(V_glasser378==n);
%         patch(isosurface(sregion,0.3), 'FaceColor', 'r', 'EdgeColor', 'none')%'FaceAlpha', abs(V(n)))      
%     end
% end
% for n=find(V<0)
%     if find(181:369==n)
%         sregion=smooth3(V_glasser378==n);
%         patch(isosurface(sregion,0.3), 'FaceColor', 'b', 'EdgeColor', 'none')%'FaceAlpha', abs(V(n))) 
%     end
% end
% 
% rotate3d
% material dull
% lighting gouraud
% % Left side view
% view(0,0)
% axis equal
% camlight;
% xlim([5 105])
% ylim([5 85])
% zlim([10 80])
% axis off
% 
% % Image of the right hemisphere
% subplot_tight(2,2,4,0.0001)
% hold on
% for n=find(V>=0)
%     if find([1:180 370:378]==n)
%         sregion=smooth3(V_glasser378==n);
%         patch(isosurface(sregion,0.3), 'FaceColor', 'r', 'EdgeColor', 'none')%'FaceAlpha', abs(V(n)))     
%     end
% end
% for n=find(V<0)
%     if find([1:180 370:378]==n)
%         sregion=smooth3(V_glasser378==n);
%         patch(isosurface(sregion,0.3), 'FaceColor', 'b', 'EdgeColor', 'none')%'FaceAlpha', abs(V(n))) 
%     end
% end
% 
% rotate3d
% material dull
% lighting gouraud
% % Left side view
% view(-180,0)
% axis equal
% camlight;xlim([5 105])
% ylim([5 85])
% zlim([10 80])
% axis off
% 
% 

