function plot_in_Glasser378 (V)

load Glasser360_2mm.mat

V=V/max(abs(V));
V=round(V*100)/100; 

sz=size(V_glasser378);
%figure('Color','k')

cmap=jet(length(-1:0.01:1));
ind_color=-1:0.01:1;

hold on
for n=[1:180 190:369]
    ind_N=find(V_glasser378==n);
    [I1,I2,I3] = ind2sub(sz,ind_N);
    plot3(I1,I2,I3,'.','color',cmap(find(ind_color-V(n)>=0,1),:))
end
for n=[181:189 370:378]
    ind_N=find(V_glasser378==n);
    [I1,I2,I3] = ind2sub(sz,ind_N);
    plot3(I1,I2,I3-50,'.','markersize',10,'color',cmap(find(ind_color-V(n)>=0,1),:))
end

axis equal
xlim([0 100])
ylim([0 200])
zlim([-50 100])
axis off
