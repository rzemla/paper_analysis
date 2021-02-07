function [outputArg1,outputArg2] = master_plotter_fig4(perf_data_plotting, frac_tuned_across_days)


%% Figure 4 - short term learning vs. recall


%unload learning performance data
mean_learn_filt_day = perf_data_plotting.stl.mean_learn_filt_day;
sem_learn_filt_day = perf_data_plotting.stl.sem_learn_filt_day;
learning_stage_labels = perf_data_plotting.stl.learning_stage_labels;
%standalone performance plot for learning cohort (try integarting as part
%of the whole plot for Figure 4

fig = figure;
fig.Units = 'centimeters';
fig.Position(1) = 8;
fig.Position(2) = 1;
fig.Position(3) = 20;
fig.Position(4) = 30;

%master layout
gridSize = [4,3];
t1 = tiledlayout(fig,gridSize(1),gridSize(2),'TileSpacing','normal','Padding','compact','Units','centimeters');

%common subplot
%s1 = tiledlayout(t1,1,2,'TileSpacing','compact','Padding','compact','Units','centimeters');
%s1.Layout.Tile = 1;
%s1.Layout.TileSpan = [1,1];
%subplot title
%set_subplot_title(s1,'Common', 'Arial', 16,'bold')

nexttile(t1,1)
hold on
axis square
xlim([0.5 7.5])
ylim([0 1.1])
yticks(0:0.2:1)
ylabel('Fraction of correct trials')
xticks((1:7))
xticklabels(learning_stage_labels)
xtickangle(45)
errorbar(1:7,mean_learn_filt_day(1,1:7),sem_learn_filt_day(1,1:7),'Color', [139, 0, 139]/255, 'LineStyle', '-','LineWidth',1.5)

%unload fractional data
mean_learning_ts_day = frac_tuned_across_days.mean_learning_ts_day;
mean_recall_ts = frac_tuned_across_days.mean_recall_ts;

bar_colorset = [65,105,225; 220,20,60; 139, 0, 139; 128 128 128]./255;

nexttile(t1,2)
hold on
axis square
ylim([0 1.2])
xlabel('Day')
ylabel('Fraction tuned')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
title('Learning - day - T.S.')
sh1 = bar(1:7,mean_learning_ts_day(1:7,:),'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh1,2)
    sh1(ii).FaceColor = 'flat';
    sh1(ii).CData = bar_colorset(ii,:);
end

nexttile(t1,3)
hold on
axis square
ylim([0 1.2])
xlabel('Day')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
xticklabels({'1','2','3','6','7','8','9'})
title('Recall - day - T.S.')
sh2 = bar(mean_recall_ts,'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh2,2)
    sh2(ii).FaceColor = 'flat';
    sh2(ii).CData = bar_colorset(ii,:);
end

%set axis font/label and font size
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12, ...
'FontWeight','normal', 'LineWidth', 1.5,'layer','top')

%% Equivalent SI data plot (Ex. Data Fig. 9)


%% Long term recall supplementary figure (Ex. Data Fig. 10)


end

