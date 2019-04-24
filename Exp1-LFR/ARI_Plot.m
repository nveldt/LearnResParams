addpath('../src/')

mus = 20:5:50;

LearnAri = zeros(5,numel(mus));
ModAri = LearnAri;
next = 1;
for mu = mus

    load(strcat('Output/Scores_ari_nmi_SmallComm_mu_',num2str(mu),'_Gurobi.mat'))
    LearnAri(:,next) = Learned(:,2);
    ModAri(:,next) = Modularity(:,2);
    next = next + 1;

end

y = median(LearnAri,1);
y2 = median(ModAri,1);

%% Display the median
figure(1); clf;
hax=axes; 
lw = 1.25;
design = '.-';
ms = 12;
plot(mus,y,design,'LineWidth',lw,'markersize',ms,'color','b');
hold on
plot(mus,y2,design,'LineWidth',lw,'markersize',ms,'color','r');


%% Show all

mus = .01*(20:5:50);
ms = 12;

LearnAri = zeros(5,numel(mus));
ModAri = LearnAri;
next = 1;
for mu = 20:5:50
    
    load(strcat('Output/Scores_ari_nmi_SmallComm_mu_',num2str(mu),'_Gurobi.mat'))
    LearnAri(:,next) = Learned(:,2);
    ModAri(:,next) = Modularity(:,2);
    next = next + 1;

end

figure(1); clf;
y = mean(LearnAri,1);
y2 = mean(ModAri,1);
X = LearnAri';
l = prctile(X',0);
u = prctile(X',100);
xvals = mus;
p1 = patch([xvals(:); flipud(xvals(:))], [l(:); flipud(u(:))],0.5*[0,0,1]+0.5*[1,1,1],'EdgeColor','none','HandleVisibility','off');
hold on
X = ModAri';
l = prctile(X',0);
u = prctile(X',100);
p2 = patch([xvals(:); flipud(xvals(:))], [l(:); flipud(u(:))],0.5*[1,0,0]+0.5*[1,1,1],'EdgeColor','none','HandleVisibility','off');
p2.FaceAlpha = .5;
p1.FaceAlpha = .5;


plot(mus,y,design,'LineWidth',lw,'markersize',ms,'color','b');
plot(mus,y2,design,'LineWidth',lw,'markersize',ms,'color','r');

legend('Location','Southwest');
legend boxoff
xlabel('\mu','fontsize',15);
ylabel({'ARI scores with ground truth'});
legend('Learning \lambda','Modularity')
xlim([.2,mus(end)])
ylim([.5,1.02])
hold off

%% Save plot if desired
set_figure_size([2.25*2,1.75*2]);
saveas(gcf,'Figures/LRFsmallcom_ARI.png')
saveas(gcf,'Figures/LRFsmallcom_ARI','jpeg')