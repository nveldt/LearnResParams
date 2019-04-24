addpath('Output')
datasets = {'DBLP','Amazon','LiveJournal','Orkut'};

clc

fprintf(' Graph | F1 Scores | Conductance | Runtimes | Sizes| \n')
for graph = datasets
   
  load(strcat('Output/',char(graph),'_5X.mat'))
  
  % Stats format:
  % mincond_stats[:,commID] = [t1; length(S1); c1; pr1; re1; f11; f21]
  % learnalpha_stats[:,commID] = [t2; length(Sbest); c2; pr2; re2; f12; f22; alphaBest; BestP]

  
  % Conductance averages
  mc_cond = mean(mincond_stats(3,:));
  la_cond = mean(learnalpha_stats(3,:));

  % F1 averages
  mc_f1 = mean(mincond_stats(6,:));
  la_f1 = mean(learnalpha_stats(6,:));
  
  % Time averages
  mc_time = mean(mincond_stats(1,:));
  la_time = mean(learnalpha_stats(1,:));
  
  % size averages
  mc_size = round(mean(mincond_stats(2,:)));
  la_size = round(mean(learnalpha_stats(2,:)));
  

  fprintf('%s & %.2f & %.2f & %.2f & %.2f & %.1f & %.1f & %d & %d \\\\\n',char(graph), mc_f1, la_f1, mc_cond, la_cond,mc_time, la_time, mc_size, la_size)

   
  
end
  