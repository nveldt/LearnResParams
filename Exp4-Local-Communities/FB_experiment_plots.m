
% # Cat1 = local_info(:,1); % student-faculty flag
% # Cat2 = local_info(:,2); % gender
% # Cat3 = local_info(:,3); % major
% # Cat4 = local_info(:,4); % 2nd major
% # Cat5 = local_info(:,5); % dorm/residence
% # Cat6 = local_info(:,6); % Year
% # Cat7 = local_info(:,7); % high school

addpath('../include')
fid = fopen('Facebook_Sets.txt','r');
f = fgets(fid);
names = split(f);
fclose(fid);
All = [];
uni = 0;

%% Collect and filter data

for i = 1:numel(names)
    name = char(names(i));
    uni = uni+1;
    try
        load(strcat('all_output_saved/FB_metacom_',name,'_20_3.mat'))
        ss = size(outputdata,1) -1;
        outputdata = [outputdata(2:end,:) uni*ones(ss,1)];
        set = find(outputdata(:,4) < .7);
        All = [All; outputdata(set,:)];
        % [category, label of group, group size, conductance of group, learned alphaBest, min of parameter fitness, precision, recall, f1_minP,f1_mincond]'

    catch
    end
end

%% Plot 1
All = sortrows(All,4,'ascend');
category = 2;
cat = find(All(:,1) > 0);

p = 1;
sz = 20;
color =  ['k', 'r', 'g', 'm','y'];
cs = [2 5 1];
for ic = 1:3
    c = cs(ic);
    cat = find(All(:,1) == c);
    numplot = round(numel(cat)*p);
    range = cat(1:numplot);
    figure(10)
    scatter(All(range,6),All(range,9),sz,'filled') %,'MarkerFaceColor',color(ic))
    hold on
 
end


figure(10)
legend({'Gender','Residence','S/F'})
legend('Location','Southeast');
ylabel('F1: learning $\alpha$','Interpreter','latex');
xlabel('$\mathcal{P}_X(\alpha)$','Interpreter','latex');

%%
set_figure_size([2.25*1.5,1.75*1.5]);
print(gcf,sprintf('Figures/fb_gender_res_sf_P.eps'),'-depsc2');
%Process_AtendHeader('Figures/fb_gender_res_sf_P.eps','');
hold off

%% Plot 2

for ic = 1:3
    c = cs(ic);
    cat = find(All(:,1) == c);
    numplot = round(numel(cat)*p);
    range = cat(1:numplot);
    
    figure(4)
    scatter(All(range,4),All(range,10),sz,'filled') %,'MarkerFaceColor',color(ic))

    hold on
end

figure(4)
legend({'Gender','Residence','S/F'})
legend('Location','Southwest');
ylabel('F1: learning $\alpha$','Interpreter','latex');
xlabel('cond(X)','Interpreter','latex');
hold off

set_figure_size([2.25*1.5,1.75*1.5]);
print(gcf,sprintf('Figures/fb_gender_res_sf_cond.eps'),'-depsc2');


%% Just 2006-2009 for graduation years
year = find(All(:,1) == 6);


y6 = find(All(:,2) == 2006);
y7 = find(All(:,2) == 2007);
y8 = find(All(:,2) == 2008);
y9 = find(All(:,2) == 2009);
not9 = [y6; y7; y8];
figure(2)
scatter(All(not9,6),All(not9,9),sz,'filled','MarkerFaceColor',[0 .75 0])
hold on
% scatter(All(y7,6),All(y7,9),sz,'filled','MarkerFaceColor','b')
% scatter(All(y8,6),All(y8,9),sz,'filled','MarkerFaceColor','b')
scatter(All(y9,6),All(y9,9),sz,'filled','MarkerFaceColor',[0 0 .75])
legend({'2006-2008','2009'})
legend('Location','Southeast');
% ylabel('F1: learning \alpha');
% xlabel({'P_X(\alpha)'});
ylabel('F1: learning $\alpha$','Interpreter','latex');
xlabel('$\mathcal{P}_X(\alpha)$','Interpreter','latex');
hold off

set_figure_size([2.25*1.5,1.75*1.5]);
print(gcf,sprintf('Figures/fb_year_P.eps'),'-depsc2');

%% Plot 4

figure(3)
scatter(All(not9,4),All(not9,9),sz,'filled','MarkerFaceColor',[0 .75 0])
hold on
% scatter(All(y7,4),All(y7,9),sz,'filled','MarkerFaceColor','b')
% scatter(All(y8,4),All(y8,9),sz,'filled','MarkerFaceColor','b')
scatter(All(y9,4),All(y9,9),sz,'filled','MarkerFaceColor',[0 0 .75])
legend({'2006-2008','2009'})
legend('Location','Southeast');
ylabel('F1: learning $\alpha$','Interpreter','latex');
xlabel('cond(X)','Interpreter','latex');
hold off

set_figure_size([2.25*1.5,1.75*1.5]);
print(gcf,sprintf('Figures/fb_year_cond.eps'),'-depsc2');
%Process_AtendHeader('Figures/fb_year_cond.eps','');

