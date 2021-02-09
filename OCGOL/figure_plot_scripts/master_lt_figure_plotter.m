function [outputArg1,outputArg2] = master_lt_figure_plotter(frac_tuned_across_days,perf_data_plotting,PV_TC_plot_data,...
                                    performance_mean_sem)



%% Figure layout


fig = figure;
fig.Units = 'centimeters';
fig.Position(1) = 7;
fig.Position(2) = 0;
fig.Position(3) = 30;
fig.Position(4) = 38;

%master layout
gridSize = [4,3];
t1 = tiledlayout(fig,gridSize(1),gridSize(2),'TileSpacing','compact','Padding','compact','Units','centimeters');

frac_tuned_across_days.mean_long_recall_si

%fraction SI
%color theme
bar_colorset = [65,105,225; 220,20,60; 139, 0, 139; 128 128 128]./255;

%fraction tuned
nexttile(t1,1)
hold on
axis square
ylim([0 1.2])
xlim([0 7])
xlabel('Day')
ylabel('Fraction tuned')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
xticklabels({1,6,16,20,25,30})
%title('Learning - day - T.S.')
sh1 = bar(1:6,frac_tuned_across_days.mean_long_recall_si(1:6,:),'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh1,2)
    sh1(ii).FaceColor = 'flat';
    sh1(ii).CData = bar_colorset(ii,:);
end

%fraction TS
nexttile(t1,2)
hold on
axis square
ylim([0 1.2])
xlim([0 7])
xlabel('Day')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:6)
xticklabels({1,6,16,20,25,30})
%title('Recall - day - T.S.')
sh2 = bar(frac_tuned_across_days.mean_long_recall_ts,'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh2,2)
    sh2(ii).FaceColor = 'flat';
    sh2(ii).CData = bar_colorset(ii,:);
end

legend({'A','B','A&B','Neither'},'Location','northoutside','Orientation','horizontal','NumColumns',2)

%unload long term performance data
mean_lt_recall_perf = perf_data_plotting.ltr.mean_lt_recall_filt_day(:,[1,6,16,20,25,30]);
sem_lt_recall_perf = perf_data_plotting.ltr.sem_lt_recall_filt_day(:,[1,6,16,20,25,30]);
lt_perf_labels = perf_data_plotting.ltr.long_recall_labels;

%performance
nexttile(t1,3)
hold on
axis square
xlim([0.5 6.5])
ylim([0 1.1])
yticks(0:0.2:1)
ylabel('Fraction of correct trials')
xlabel('Day')
xticks((1:6))
xticklabels(lt_perf_labels)
%xtickangle(45)
errorbar(1:6,mean_lt_recall_perf(1,1:6),sem_lt_recall_perf (1,1:6),'Color', [139, 0, 139]/255, 'LineStyle', '-','LineWidth',1.5)

%PV

nexttile(t1,4)
hold on
axis square
xlabel('Day')
ylabel('Correlation score')
xticks([1,6,16,20,25,30])
xlim([0 32]) 
ylim([0 1])

%recall
rA = plot_error_line(PV_TC_plot_data.PV_mean_sem.lt_recall.A,'-',2,[65,105,225]/255);
rB = plot_error_line(PV_TC_plot_data.PV_mean_sem.lt_recall.B,'-',2,[220,20,60]/255);

legend([rA,rB],{'A','B'},'location','northeast')
%TC SI

nexttile(t1,5)
hold on
axis square
xlabel('Day')
ylabel('Correlation score')
xticks([1,6,16,20,25,30])
xlim([0 32]) 
ylim([0 1])

%recall
rA = plot_error_line(PV_TC_plot_data.TC_mean_sem.si.lt_recall.A,'-',2,[65,105,225]/255);
rB = plot_error_line(PV_TC_plot_data.TC_mean_sem.si.lt_recall.B,'-',2,[220,20,60]/255);

legend([rA,rB],{'A','B'},'location','northeast')

%TC TS
nexttile(t1,6)
hold on
axis square
xlabel('Day')
ylabel('Correlation score')
xticks([1,6,16,20,25,30])
xlim([0 32]) 
ylim([0 1])

%recall
rA = plot_error_line(PV_TC_plot_data.TC_mean_sem.ts.lt_recall.A,'-',2,[65,105,225]/255);
rB = plot_error_line(PV_TC_plot_data.TC_mean_sem.ts.lt_recall.B,'-',2,[220,20,60]/255);

legend([rA,rB],{'A','B'},'location','northeast')

%PV neighbor
nexttile(t1,7)
axis square
hold on
xlim([0 6])
ylim([0 1])
xlabel('Day comparison')
ylabel('Correlation score')
xticks([1:5])
xlim([0 6])
xticklabels({'1 vs. 6','6 vs. 16','16 vs. 20','20 vs. 25','25 vs. 30'})
xtickangle(45)
%recall
rA = plot_error_line(PV_TC_plot_data.PV_mean_sem.neighbor.lt_recall.A,'-',2,[65,105,225]/255);
rB = plot_error_line(PV_TC_plot_data.PV_mean_sem.neighbor.lt_recall.B,'-',2,[220,20,60]/255);

% set(gca,'FontSize',16)
% set(gca,'Linewidth',2)

legend([rA,rB],{'A','B'},'location','northeast')

%TC neighbor SI
nexttile(t1,8)
axis square
hold on
xlim([0 6])
ylim([0 1])
xlabel('Day comparison')
ylabel('Correlation score')
xticks([1:5])
xlim([0 6])
xticklabels({'1 vs. 6','6 vs. 16','16 vs. 20','20 vs. 25','25 vs. 30'})
xtickangle(45)
%recall
rA = plot_error_line(PV_TC_plot_data.TC_mean_sem.neighbor.si.lt_recall.A,'-',2,[65,105,225]/255);
rB = plot_error_line(PV_TC_plot_data.TC_mean_sem.neighbor.si.lt_recall.B,'-',2,[220,20,60]/255);

% set(gca,'FontSize',16)
% set(gca,'Linewidth',2)

legend([rA,rB],{'A','B'},'location','northeast')
%TC neighbor TS
nexttile(t1,9)
axis square
hold on
xlim([0 6])
ylim([0 1])
xlabel('Day comparison')
ylabel('Correlation score')
xticks([1:5])
xlim([0 6])
xticklabels({'1 vs. 6','6 vs. 16','16 vs. 20','20 vs. 25','25 vs. 30'})
xtickangle(45)
%recall
rA = plot_error_line(PV_TC_plot_data.TC_mean_sem.neighbor.ts.lt_recall.A,'-',2,[65,105,225]/255);
rB = plot_error_line(PV_TC_plot_data.TC_mean_sem.neighbor.ts.lt_recall.B,'-',2,[220,20,60]/255);

% set(gca,'FontSize',16)
% set(gca,'Linewidth',2)

legend([rA,rB],{'A','B'},'location','northeast')

%A&B correlation SI
nexttile(t1,10)
hold on
axis square
yyaxis left
ylim([0 1.2])
yticks([0:0.2:1])
xlabel('Relative day')
ylabel('Normalized correlation score')
        xticks([1:6])
        xlim([0 7])
        xticklabels({'1','6','16','20','25','30'})
%dashed 1 reference line
plot([0 10],[1 1],'--','Color',[0.5 0.5 0.5])

rA = plot_error_line(PV_TC_plot_data.ABcorr.si.lt_recall.AB  ,'-',2,[139,0,139]/255);

%plot correlation on right y axis
yyaxis right
ylabel('Performance')
ylim([0 1.2])
yticks([0.2 0.4 0.6 0.8 1])

plot_error_line(performance_mean_sem.lt_recall,'-',2,[34,139,34]/255);

set(gca,'FontSize',12)
set(gca,'Linewidth',2)

ax = gca;
%set left axis color
ax.YAxis(1).Color = [139,0,139]/255;
%set right axis color
ax.YAxis(2).Color = [34,139,34]/255;

%A&B correlation TS
nexttile(t1,11)
hold on
axis square
yyaxis left
yticks([0:0.2:1])
ylim([0 1.2])
xlabel('Relative day')
ylabel('Normalized correlation score')
        xticks([1:6])
        xlim([0 7])
        xticklabels({'1','6','16','20','25','30'})
%dashed 1 reference line
plot([0 10],[1 1],'--','Color',[0.5 0.5 0.5])

rA = plot_error_line(PV_TC_plot_data.ABcorr.ts.lt_recall.AB  ,'-',2,[139,0,139]/255);

%plot correlation on right y axis
yyaxis right
ylabel('Performance')
ylim([0 1.2])
yticks([0.2 0.4 0.6 0.8 1])

plot_error_line(performance_mean_sem.lt_recall,'-',2,[34,139,34]/255);

set(gca,'FontSize',12)
set(gca,'Linewidth',2)

ax = gca;
%set left axis color
ax.YAxis(1).Color = [139,0,139]/255;
%set right axis color
ax.YAxis(2).Color = [34,139,34]/255;

%set axis font/label and font size
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12, ...
    'FontWeight','normal', 'LineWidth', 1.5,'layer','top')


end

