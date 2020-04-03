function plot_arrows_in_cortex_general (Parcellation,V)

hold on
Directory='/Users/joana/Documents/Work/LEiDA general/';
addpath(genpath(Directory))


% PLOT TRANSPARENT CORTEX
Brain_Mask_tmp=load_nifti('MNI152_T1_2mm_brain_mask.nii');
Brain_Mask = Brain_Mask_tmp.vol;
scortex=smooth3(Brain_Mask>0);
cortexpatch=patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
reducepatch(cortexpatch,0.01);
isonormals(scortex,cortexpatch);

% PLOT NODES IN THE CENTER OF GRAVITY OF EACH ROI
N_areas=length(V);
V_Parcels=struct2array(load('ParcelsMNI2mm',['V_' Parcellation]));
V_Parcels((V_Parcels>N_areas))=0;
sz=size(V_Parcels);
for n=1:max(V_Parcels(:))
    ind_N=find(V_Parcels==n);
    [I1,I2,I3] = ind2sub(sz,ind_N);
    Node_cog(n,:)=mean([I1,I2,I3]);
end
V=V/max(abs(V));
V=round(V*100)/100;

scale=3;
P1=[Node_cog(:,2) Node_cog(:,1) Node_cog(:,3)];
P2=P1+[4*sign(V)'+scale*V' zeros(N_areas,1)  zeros(N_areas,1)];
daspect([1 1 1])
blueharrow3=arrow3(P1(V<=0,:),P2(V<=0,:),'1.2b',1,.7,.3);
if find(V>0)
    redharrow3=arrow3(P1(V>0,:),P2(V>0,:),'1.2r',1,.7,.3);
end
% axis off;
% view(-90,90)
% camlight;
% %rotate3d;
% %daspect([1 1 1])
