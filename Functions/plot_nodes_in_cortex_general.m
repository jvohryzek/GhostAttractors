function plot_nodes_in_cortex_general (Parcellation,V)

hold on
Directory='/Users/joana/Documents/Work/LEiDA general/LEiDA_HCP/';

% PLOT TRANSPARENT CORTEX
cortex=niftiread([Directory 'MNI152_T1_2mm_brain_mask.nii']);
sregion=smooth3(cortex);
psregion=patch(isosurface(sregion,0.3), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
reducepatch(psregion,0.1);
isonormals(sregion,psregion);
set(psregion,'FaceAlpha', 0.1); %transparency

% PLOT NODES IN THE CENTER OF GRAVITY OF EACH ROI
N_areas=length(V);
V_Parcels=struct2array(load([Directory 'ParcelsMNI2mm'],['V_' Parcellation]));
V_Parcels((V_Parcels>N_areas))=0;
sz=size(V_Parcels);
for n=1:max(V_Parcels(:))
    ind_N=find(V_Parcels==n);
    [I1,I2,I3] = ind2sub(sz,ind_N);
    Node_cog(n,:)=mean([I1,I2,I3]);
end
V=V/max(abs(V));
V=round(V*100)/100; 

a=2.5;
[x,y,z] = sphere;
x=a*x;
y=a*y;
z=a*z;

colormap(jet)
ind_color=-1:0.01:1;
cmap=jet(numel(ind_color));

for n=1:length(V)
        surf(x+Node_cog(n,2), y+Node_cog(n,1),z+Node_cog(n,3),'FaceColor',cmap(find(abs(ind_color-V(n))<0.01,1),:),'EdgeColor','none','FaceAlpha',0.7);
end
  
% PLOT LINKS BETWEEN NODES 
n_strong=find(V>0);
if numel(n_strong)>1
    u=1;
    
    for a=1:numel(n_strong)
        n=n_strong(a);
        for b=1:a
            p=n_strong(b);
            c1=[Node_cog(n,2) Node_cog(n,1) Node_cog(n,3)];
            c2=[Node_cog(p,2) Node_cog(p,1) Node_cog(p,3)];
            
            plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)],'Color',[1 0.3 0.3]); %cmap(find(ind_colormap-(V(n)+V(p))/2>0,1),:)); 
            u=u+1;
        end
    end
end

% Setting image properties - light, material, angle
 
material dull
lighting gouraud
% Top view
view(-90,90)
daspect([1 1 1])
camlight;
xlim([5 105])
ylim([5 85])
zlim([10 80])
axis off
end