fid = fopen('Facebook_Sets.txt','r');
addpath('../include')
% Must add path to Facebook100 dataset
addpath('~/data/Facebook100/')
f = fgets(fid);
names = split(f);
minSize = 20;
YearSets = zeros(700,7);
next = 1;

% For each network
for i = 1:numel(names)
    name = names{i};
    load(strcat(name))
    volA = sum(nonzeros(A));
    
    % Get the graduation year labels, exclude nodes with a 0 label
    yearlabels = unique(local_info(:,6));
    yearlabels = setdiff(yearlabels,0);
    
    % For each graduation year...
    numyears = numel(yearlabels);
    for j = 1:numel(yearlabels)
        y = yearlabels(j);
        
        % ...get stats for the set of nodes with that label in the network
        year = find(local_info(:,6) == y);
        if numel(year) > minSize
            [cut,vol,edges,cond] = set_stats(A,year,volA);
            YearSets(next,:) = [y,i,numel(year),cut, vol, edges, cond];
            next = next+1;
        end
    end
    
end

%% Bar plots
fprintf('\n\n')
years = unique(YearSets(:,1));
sz = 10;
mainyears = intersect(find(YearSets(:,1) > 2005),find(YearSets(:,1) < 2010));

vols = YearSets(mainyears,5);
cuts = YearSets(mainyears,4);
conds = YearSets(mainyears,7);

ind6 = find(YearSets(mainyears,1) == 2006);
ind7 = find(YearSets(mainyears,1) == 2007);
ind8 = find(YearSets(mainyears,1) == 2008);
ind9 = find(YearSets(mainyears,1) == 2009);


figure(1)
boxplot([vols(ind6),vols(ind7),vols(ind8), vols(ind9)])
xlabel('Graduation year')
ylabel('volume(X)')
hold off
scale = 1.25
set_figure_size([2.25*scale,1.75*scale]);
print(gcf,sprintf('Figures/fb_year_vol.eps'),'-depsc2');

figure(2)
boxplot([cuts(ind6),cuts(ind7),cuts(ind8), cuts(ind9)])
xlabel('Graduation year')
ylabel('cut(X)')
hold off
set_figure_size([2.25*scale,1.75*scale]);
print(gcf,sprintf('Figures/fb_year_cut.eps'),'-depsc2');

figure(3)
boxplot([conds(ind6),conds(ind7),conds(ind8), conds(ind9)])
xlabel('Graduation year')
ylabel('\phi(X)')
hold off
set_figure_size([2.25*scale,1.75*scale]);
print(gcf,sprintf('Figures/fb_year_maincond.eps'),'-depsc2');

Vols = [vols(ind6),vols(ind7),vols(ind8), vols(ind9)];
Cuts = [cuts(ind6),cuts(ind7),cuts(ind8), cuts(ind9)];
Conds = [conds(ind6),conds(ind7),conds(ind8), conds(ind9)];


save('FBvolcutcond_data','Vols','Cuts','Conds')

%%  Violin plots

scale = 1.5;
addpath('../include/Violinplot-Matlab-master/')
figure(1)
vs = violinplot(Vols, {'2006','2007','2008','2009'},'ShowData',false); %,'ViolinColor',[0 1 0; 0 0 1; 1 0 0; 1 0 0])
xlabel('Graduation year')
ylabel('volume(X)')
hold off
set_figure_size([2.25*scale,1.75*scale]);
print(gcf,sprintf('Figures/fb_year_vol.eps'),'-depsc2');
saveas(gcf,'Figures/fb_year_vol.pdf')

figure(2)
vs = violinplot(Cuts, {'2006','2007','2008','2009'},'ShowData',false); %,'ViolinColor',[0 1 0; 0 0 1; 1 0 0; 1 0 0])
xlabel('Graduation year')
ylabel('cut(X)')
hold off
set_figure_size([2.25*scale,1.75*scale]);
print(gcf,sprintf('Figures/fb_year_cut.eps'),'-depsc2');
saveas(gcf,'Figures/fb_year_cut.pdf')

figure(3)
vs = violinplot(Conds, {'2006','2007','2008','2009'},'ShowData',false); %,'ViolinColor',[0 1 0; 0 0 1; 1 0 0; 1 0 0])
xlabel('Graduation year')
ylabel('\phi(X)')
hold off
set_figure_size([2.25*scale,1.75*scale]);
print(gcf,sprintf('Figures/fb_year_maincond.eps'),'-depsc2');
saveas(gcf,'Figures/fb_year_maincond.pdf')