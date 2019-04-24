datasets = {'DBLP','Amazon','LiveJournal','Orkut'};
names = {'dblp','amazon','livejournal','orkut'};

  % Stats format:
  % mincond_stats[:,commID] = [t1; length(S1); c1; pr1; re1; f11; f21]
  % learnalpha_stats[:,commID] = [t2; length(Sbest); c2; pr2; re2; f12; f22; alphaBest; BestP]
  
sz = 40;

for graph = datasets
   
  load(strcat('Output/',char(graph),'_5X.mat'))
  learnalpha_stats(9,:)
  learnalpha_stats(6,:)
    sz = 40;
    figure(1)
    scatter(learnalpha_stats(9,:),learnalpha_stats(6,:),sz,'filled')
    hold on
end

figure(1)
legend(names)
legend('Location','Northeast');
% legend boxoff
 ylabel('F1: learning \alpha');
% xlabel({'P_X(\alpha)'});
%ylabel('F1: learning $\alpha$','Interpreter','latex');
xlabel('$\mathcal{P}_X(\alpha)$','Interpreter','latex');
hold off

%%
scale = 1.25;
addpath('../src')
set_figure_size([2.25*scale,1.75*scale]);
print(gcf,'Figures/plot1_leaningalpha.eps','-depsc2');

%% Need to have the datasets properly formatted and stored in the graphs folder

for graph = datasets
   
  load(strcat('Output/',char(graph),'_5X.mat'))
    load(strcat('graphs/',char(graph),'-top10.mat'))
    sz = 40;
    figure(2)
    scatter(Conds,mincond_stats(6,:),sz,'filled')
    hold on
end

figure(2)
% legend(names)
% legend('Location','North');
% legend boxoff
ylabel('F1: min cond subset');
%xlabel('$\phi(X)$','Interpreter','latex');
xlabel('cond(X)','Interpreter','latex');
hold off
set_figure_size([2.25*scale,1.75*scale]);
print(gcf,'Figures/plot1_mincond.eps','-depsc2');
