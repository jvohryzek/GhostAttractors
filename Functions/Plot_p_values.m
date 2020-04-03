function Plot_p_values(N_areas)

load([Directory '/LEiDA_' num2str(N_areas) '_results/LEiDA' num2str(N_areas) '_K' num2str(mink) '-' num2str(maxk) 'results.mat','P_pval','LT_pval','rangeK','Kmeans_results','Index_Controls','Index_Patients','LT','P')
mink=rangeK(1);
maxk=rangeK(end);

colormap(jet)
ind_colormap=-1:0.01:1;
cmap=jet(numel(ind_colormap));


%% PLOT P-VALUES OVER THE RANGE OF K

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK) 
    for c=1:rangeK(k)        
        if P_pval(k,c)>0.05
            semilogy(rangeK(k),P_pval(k,c),'*k');
            text(rangeK(k),P_pval(k,c),[' ' num2str(c)])
        end
        if P_pval(k,c)<0.05 && P_pval(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),P_pval(k,c),'*r');
            text(rangeK(k),P_pval(k,c),[' ' num2str(c)])
        end
        if P_pval(k,c)<(0.05/rangeK(k)) && P_pval(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),P_pval(k,c),'*g');
            text(rangeK(k),P_pval(k,c),[' ' num2str(c)])
        end
        if P_pval(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),P_pval(k,c),'*b');
            text(rangeK(k),P_pval(k,c),[' ' num2str(c)])
        end
    end
end

title('Probability Patients vs Controls')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off

subplot(1,2,2)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK) 
    for c=1:rangeK(k)        
        if LT_pval(k,c)>0.05
            semilogy(rangeK(k),LT_pval(k,c),'*k');
            text(rangeK(k),LT_pval(k,c),[' ' num2str(c)])
        end
        if LT_pval(k,c)<0.05 && LT_pval(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),LT_pval(k,c),'*r');
            text(rangeK(k),LT_pval(k,c),[' ' num2str(c)])
        end
        if LT_pval(k,c)<(0.05/rangeK(k)) && LT_pval(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),LT_pval(k,c),'*g');
            text(rangeK(k),LT_pval(k,c),[' ' num2str(c)])
        end
        if LT_pval(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),LT_pval(k,c),'*b');
            text(rangeK(k),LT_pval(k,c),[' ' num2str(c)])
        end
    end
end

title('Dwell Times Controls vs Patients')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off


%%


figure('Name','Probabilities')
load AAL_labels.mat label116
Order=[1:2:N_areas N_areas:-2:2];
labels=label116(Order,:);

sigC=zeros(1,length(rangeK));
for k=1:length(rangeK)   
    [Min_p_value(k), sigC(k)]=min(P_pval(k,P_pval(k,:)>0));    
end

V_all=zeros(N_areas,length(rangeK));

for k=1:length(rangeK)
    c=sigC(k);
    V=Kmeans_results{k}.C(c,:);
    V=V/max(abs(V));
    V_all(:,k)=V;
    
    V=V(Order);
    V=V/max(abs(V));
    V=round(V*100)/100;   
    
    subplot(1,length(rangeK),k)
    hold on
    for n=1:length(V)
        barh(V.*(1:length(V)==n),'FaceColor',cmap(abs(ind_colormap-V(n))<0.00001,:),'EdgeColor','none','Barwidth',.5)
    end
    ylim([0 N_areas+1])
    xlim([-1 1])
    
    grid on
    set(gca,'YTick',1:N_areas,'Fontsize',9)    
    if k==1
        set(gca,'YTickLabel',labels(end:-1:1,:))
    else
        set(gca,'YTickLabel',[])
    end
    if Min_p_value(k)<0.05 && Min_p_value(k)>(0.05/rangeK(k))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k),'%.1e')]},'Fontsize',10,'color','r')
    elseif Min_p_value(k)>0.05
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k),'%.1e')]},'Fontsize',10,'color','k')
    elseif Min_p_value(k)<(0.05/k) && Min_p_value(k)>(0.05/sum(rangeK))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k),'%.1e')]},'Fontsize',10,'color','g')
    elseif Min_p_value(k)<(0.05/sum(rangeK))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k),'%.1e')]},'Fontsize',10,'color','b')        
    end
end

[~, kbest]=min(Min_p_value);
V_best=Kmeans_results{kbest}.C(sigC(kbest),:);
cbest=sigC(kbest);

figure

subplot(3,2,1)
plot_nodes_in_cortex(V_best)
rotate3d
view(0,0)
title({'Network with most different','Probability','between groups',['K=' num2str(rangeK(kbest)) ', c=' num2str(sigC(kbest)) ]})

subplot(4,2,3)
plot_nodes_in_cortex(V_best)
view(-90,90)
rotate3d

V_best=V_best(Order);
V_best=V_best./max(abs(V_best));
V_best=round(V_best*100)/100; 

subplot(4,2,[2 4 6 8])
hold on
for n=1:length(V_best)
    barh(V_best.*(1:length(V_best)==n),'FaceColor',cmap(abs(ind_colormap-V_best(n))<0.00001,:),'EdgeColor','none','Barwidth',.5)
end
set(gca,'YTick',1:N_areas,'Fontsize',9)
set(gca,'YTickLabel',labels(end:-1:1,:))
ylim([0 N_areas+1])
xlim([-1 1])
title('Relative BOLD Phase','Fontsize',12)
grid on

subplot(4,2,5)
colormap(jet)
FC_V=V_best'*V_best;  
li=max(abs(FC_V(:)));
imagesc(FC_V,[-li li])
axis square
title('FC pattern') 
ylabel('Brain area #')
xlabel('Brain area #')  

subplot(4,2,7)
% Plot Lifetimes
Controls=squeeze(P(Index_Controls,kbest,cbest));
Patients=squeeze(P(Index_Patients,kbest,cbest));
bar([mean(Controls) mean(Patients)],'EdgeColor','w','FaceColor',[.5 .5 .5])
hold on
% Error bar containing the standard error of the mean
errorbar([mean(Controls) mean(Patients)],[std(Controls)/sqrt(numel(Controls)) std(Patients)/sqrt(numel(Patients))],'LineStyle','none','Color','k')
set(gca,'XTickLabel',{'HEALTHY', 'OCD'})
ylabel('Probability')
box off
title(['p= ' num2str(P_pval(kbest,cbest))])
grid on


%% LIFETIMES

sigC=zeros(length(rangeK));

for k=1:length(rangeK)   
    [Min_p_value(k), sigC(k)]=min(LT_pval(k,LT_pval(k,:)>0));    
end

figure
colormap(jet)
ind_colormap=-1:0.01:1;
cmap=jet(numel(ind_colormap));

V_all=zeros(N_areas,length(rangeK));

for k=1:length(rangeK)
    c=sigC(k);
    V=Kmeans_results{k}.C(c,:);
    V_all(:,k)=V;
    
    V=V(Order);
    V=V/max(abs(V));
    V=round(V*100)/100;   
    
    subplot(1,length(rangeK),k)
    hold on
    for n=1:length(V)
        barh(V.*(1:length(V)==n),'FaceColor',cmap(abs(ind_colormap-V(n))<0.00001,:),'EdgeColor','none','Barwidth',.5)
    end
    ylim([0 N_areas+1])
    xlim([-1 1])
    
    grid on
    set(gca,'YTick',1:N_areas,'Fontsize',9)    
    if k==1
        set(gca,'YTickLabel',labels(end:-1:1,:))
    else
        set(gca,'YTickLabel',[])
    end
    if Min_p_value(k)<0.05 && Min_p_value(k)>(0.05/rangeK(k))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k),'%.1e')]},'Fontsize',10,'color','r')
    elseif Min_p_value(k)>0.05
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k),'%.1e')]},'Fontsize',10,'color','k')
    elseif Min_p_value(k)<(0.05/k) && Min_p_value(k)>(0.05/sum(rangeK))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k),'%.1e')]},'Fontsize',10,'color','g')
    elseif Min_p_value(k)<(0.05/sum(rangeK))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k),'%.1e')]},'Fontsize',10,'color','b')        
    end
    title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k),'%.1e')]},'Fontsize',10)
end

[~, kbest]=min(Min_p_value);
cbest=sigC(kbest);
V_best=Kmeans_results{kbest}.C(cbest,:);

figure

subplot(4,2,1)
plot_nodes_in_cortex(V_best)
rotate3d
view(0,0)
title({'Network with most different','Dwell Time','between groups',['K=' num2str(rangeK(kbest)) ', c=' num2str(cbest) ]})

subplot(4,2,3)
plot_nodes_in_cortex(V_best)
view(-90,90)
rotate3d

V_best=V_best(Order);
V_best=V_best./max(abs(V_best));
V_best=round(V_best*100)/100; 

subplot(4,2,[2 4 6 8])
hold on
for n=1:length(V_best)
    barh(V_best.*(1:length(V_best)==n),'FaceColor',cmap(abs(ind_colormap-V_best(n))<0.00001,:),'EdgeColor','none','Barwidth',.5)
end
set(gca,'YTick',1:N_areas,'Fontsize',9)
set(gca,'YTickLabel',labels(end:-1:1,:))
ylim([0 N_areas+1])
xlim([-1 1])
title('Relative BOLD Phase','Fontsize',12)
grid on

subplot(4,2,5)
colormap(jet)
FC_V=V_best'*V_best;  
li=max(abs(FC_V(:)));
imagesc(FC_V,[-li li])
axis square
title('FC pattern') 
ylabel('Brain area #')
xlabel('Brain area #')  

subplot(4,2,7)
% Plot Lifetimes
Controls=squeeze(LT(Index_Controls,kbest,cbest));
Patients=squeeze(LT(Index_Patients,kbest,cbest));
bar([mean(Controls) mean(Patients)],'EdgeColor','w','FaceColor',[.5 .5 .5])
hold on
% Error bar containing the standard error of the mean
errorbar([mean(Controls) mean(Patients)],[std(Controls)/sqrt(numel(Controls)) std(Patients)/sqrt(numel(Patients))],'LineStyle','none','Color','k')
set(gca,'XTickLabel',{'HEALTHY', 'OCD'})
ylabel('Dwell Time (seconds)')
box off
title(['p= ' num2str(LT_pval(kbest,cbest))])
grid on