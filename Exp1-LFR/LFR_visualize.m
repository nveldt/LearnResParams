clear
mu = 30;
    
addpath('../src')
load(strcat('Output/Scores_ari_nmi_SmallComm_mu_',num2str(mu),'_Gurobi.mat'))

example = 2;
[c,p] = sort(Truths(:,example));

clam = Clams(p,example);
cmod = Cmods(p,example);
A = As(p,p,example);

[C,p] = sortrows([c, cmod],2);

clam = clam(p);
c = C(:,1);
cmod = C(:,2);

c = rename_clustering(c);
A = A(p,p);


xy = ClusterLayout(c);
 
%% Plots
colormap = lines
d = sum(A,2);
lw = .1;
ms = 10;
sms = 35; % set marker size
next = 1;
figure(1); clf;
[px,py] = gplot(A,xy);
grayscale = .7;
plot(px,py,'.-','LineWidth',lw,'MarkerSize',ms,'Color',[grayscale grayscale grayscale]);
axis off;
axis tight;
axis square;
hold on
plot(xy(:,1),xy(:,2),'.','MarkerSize',ms,'Color','k');
nmodc = max(cmod);
for k = 1:nmodc
    set = find(cmod == k);
    plot(xy(set,1),xy(set,2),'.','MarkerSize',ms,'color',colormap(next,:));
    %next = mod((next+7),64);
    next = next+1;
end
set_figure_size([2.25*2,1.75*2]);
saveas(gcf,'Figures/LFR_smallcom_visualize_mod.png')
saveas(gcf,'Figures/LFR_smallcom_visualize_mod.pdf')
print(gcf,'Figures/LFR_smallcom_visualize_mod.eps','-depsc2');

%%

next = 1;
figure(2); clf;
[px,py] = gplot(A,xy);
grayscale = .7;
plot(px,py,'.-','LineWidth',lw,'MarkerSize',ms,'Color',[grayscale grayscale grayscale]);
axis tight;
axis square;
axis off;
hold on
plot(xy(:,1),xy(:,2),'.','MarkerSize',ms,'Color','k');
nlamc = max(clam);
clc
for k = 1:nlamc
    set = find(clam == k);
    plot(xy(set,1),xy(set,2),'.','MarkerSize',ms,'color',colormap(next,:));
    next = mod((next+2),64);
end
set_figure_size([2.25*2,1.75*2]);
saveas(gcf,'Figures/LFR_smallcom_visualize_lam.png')
print(gcf,sprintf(strcat('Figures/LFR_smallcom_visualize_lam.eps')),'-depsc2');



