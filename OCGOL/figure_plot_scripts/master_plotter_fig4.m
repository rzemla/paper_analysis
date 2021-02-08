function [outputArg1,outputArg2] = master_plotter_fig4(perf_data_plotting, frac_tuned_across_days,PV_TC_plot_data,...
    centroid_data)


%% Figure 4 - short term learning vs. recall

%unload learning performance data
mean_learn_filt_day = perf_data_plotting.stl.mean_learn_filt_day;
sem_learn_filt_day = perf_data_plotting.stl.sem_learn_filt_day;
learning_stage_labels = perf_data_plotting.stl.learning_stage_labels;
%standalone performance plot for learning cohort (try integarting as part
%of the whole plot for Figure 4

fig = figure;
fig.Units = 'centimeters';
fig.Position(1) = 7;
fig.Position(2) = 0;
fig.Position(3) = 30.8;
fig.Position(4) = 39.2;

%master layout
gridSize = [4,3];
t1 = tiledlayout(fig,gridSize(1),gridSize(2),'TileSpacing','none','Padding','compact','Units','centimeters');

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
%title('Learning - day - T.S.')
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
%title('Recall - day - T.S.')
sh2 = bar(mean_recall_ts,'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh2,2)
    sh2(ii).FaceColor = 'flat';
    sh2(ii).CData = bar_colorset(ii,:);
end

legend({'A','B','A&B','Neither'},'Location','northoutside','Orientation','horizontal')

%unload data
PV_mean_sem = PV_TC_plot_data.PV_mean_sem;
TC_mean_sem = PV_TC_plot_data.TC_mean_sem;

nexttile(t1,5)
hold on
axis square
xlabel('Relative Day to D1')
ylabel('Corr. score')
xticks([1:9])
xticklabels({'0','1','2','3','4','5','6','7','8'})
%title('PV correlation - raw')
%learn
lA = plot_error_line(PV_mean_sem.st_learn.A,'--',2,[65,105,225]/255);
lB = plot_error_line(PV_mean_sem.st_learn.B,'--',2,[220,20,60]/255);

%recall
rA = plot_error_line(PV_mean_sem.st_recall.d4_d5.A,'-',2,[65,105,225]/255);
rB = plot_error_line(PV_mean_sem.st_recall.d4_d5.B,'-',2,[220,20,60]/255);

legend([lA,lB,rA,rB],{'Learning A','Learning B','Recall A','Recall B'},'location','northeast')

nexttile(t1,6)
hold on
axis square
xlabel('Relative Day to D1')
ylabel('Corr. score')
xticks([1:9])
xticklabels({'0','1','2','3','4','5','6','7','8'})
%title('TC TS correlation - raw')
%learn
lA = plot_error_line(TC_mean_sem.ts.st_learn.A,'--',2,[65,105,225]/255);
lB = plot_error_line(TC_mean_sem.ts.st_learn.B,'--',2,[220,20,60]/255);

%recall
rA = plot_error_line(TC_mean_sem.ts.st_recall.d4_d5.A,'-',2,[65,105,225]/255);
rB = plot_error_line(TC_mean_sem.ts.st_recall.d4_d5.B,'-',2,[220,20,60]/255);

legend([lA,lB,rA,rB],{'Learning A','Learning B','Recall A','Recall B'},'location','northeast')

nexttile(t1,8)
hold on
axis square
ylim([0 1])
xlabel('Day comparison')
ylabel('Correlation score')
xticks([1:8])
xlim([0 9])
xticklabels({'1 vs. 2','2 vs. 3','3 vs. 4','4 vs. 5','5 vs. 6','6 vs. 7','7 vs. 8','8 vs. 9'})
xtickangle(45)
%title('PV correlation neighbor')
%learn
lA = plot_error_line(PV_mean_sem.neighbor.st_learn.A(:,1:6),'--',2,[65,105,225]/255);
lB = plot_error_line(PV_mean_sem.neighbor.st_learn.B(:,1:6),'--',2,[220,20,60]/255);
%recall
rA = plot_error_line(PV_mean_sem.neighbor.st_recall.A,'-',2,[65,105,225]/255);
rB = plot_error_line(PV_mean_sem.neighbor.st_recall.B,'-',2,[220,20,60]/255);

set(gca,'FontSize',16)
set(gca,'Linewidth',2)

legend([lA,lB,rA,rB],{'Learning A','Learning B','Recall A','Recall B'},'location','southeast')

nexttile(t1,9)
axis square
hold on
xlim([0 9])
ylim([0 1])
xlabel('Day comparison')
ylabel('Correlation score')
xticks([1:8])
xlim([0 9])
xticklabels({'1 vs. 2','2 vs. 3','3 vs. 4','4 vs. 5','5 vs. 6','6 vs. 7','7 vs. 8','8 vs. 9'})
xtickangle(45)
%title('TC correlation - TS neighbor')
%learn
lA = plot_error_line(TC_mean_sem.neighbor.ts.st_learn.A(:,1:6),'--',2,[65,105,225]/255);
lB = plot_error_line(TC_mean_sem.neighbor.ts.st_learn.B(:,1:6),'--',2,[220,20,60]/255);
%recall
rA = plot_error_line(TC_mean_sem.neighbor.ts.st_recall.A,'-',2,[65,105,225]/255);
rB = plot_error_line(TC_mean_sem.neighbor.ts.st_recall.B,'-',2,[220,20,60]/255);

set(gca,'FontSize',16)
set(gca,'Linewidth',2)

legend([lA,lB,rA,rB],{'Learning A','Learning B','Recall A','Recall B'},'location','southeast')

set(gca,'FontSize',16)
set(gca,'Linewidth',2)

%by animal centroid difference here
%unload centroid data
%A
st_learn.animal.ts.A = centroid_data.st_learn.ts.A;
st_recall.animal.ts.A = centroid_data.st_recall.ts.A;
%B
st_learn.animal.ts.B = centroid_data.st_learn.ts.B;
st_recall.animal.ts.B = centroid_data.st_recall.ts.B;


nexttile(t1,11)
cm2rad = (2*pi)./196;
%rad ticks that correspond to 0 25 50 cm
rad_ticks = [0 10 20 30 40].*cm2rad;

hold on
axis square
hold on
xlabel('Relative day')
ylabel('Centroid distance [cm]')
ylim([0 1.4])
yticks(rad_ticks)
yticklabels({'0','10','20','30','40'});

xticks([2:9])
xlim([1 10])
xticklabels({'1','2','3','4','5','6','7','8'})
%learn
lA = plot_error_line(st_learn.animal.ts.A(:,2:end),'--',2,[65,105,225]./255);
%recall
rA = plot_error_line(st_recall.animal.ts.A(:,2:end),'-',2,[65,105,225]./255);

legend([lA,rA],{'Learning A','Recall A'},'location','southeast')

nexttile(t1,12)
hold on
axis square
hold on
xlabel('Relative day')
ylabel('Centroid distance [cm]')
ylim([0 1.4])
yticks(rad_ticks)
yticklabels({'0','10','20','30','40'});

xticks([2:9])
xlim([1 10])
xticklabels({'1','2','3','4','5','6','7','8'})
%learn
lB = plot_error_line(st_learn.animal.ts.B(:,2:end),'--',2,[220,20,60]./255);
%recalll
rB = plot_error_line(st_recall.animal.ts.B(:,2:end),'-',2,[220,20,60]./255);

legend([lB,rB],{'Learning B','Recall B'},'location','southeast')

%set axis font/label and font size
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12, ...
    'FontWeight','normal', 'LineWidth', 1.5,'layer','top')

%% Equivalent SI data plot (Ex. Data Fig. 9)


%% Long term recall supplementary figure (Ex. Data Fig. 10)


end

