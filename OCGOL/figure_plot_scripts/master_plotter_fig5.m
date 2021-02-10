function [outputArg1,outputArg2] = master_plotter_fig5(PV_TC_plot_data,performance_mean_sem)



%% Plot

fig = figure;
fig.Units = 'centimeters';
fig.Position(1) = 7;
fig.Position(2) = 0;
fig.Position(3) = 10.5;
fig.Position(4) = 24;

%master layout
gridSize = [2,1];
t1 = tiledlayout(fig,gridSize(1),gridSize(2),'TileSpacing','compact','Padding','compact','Units','centimeters');


nexttile(t1,1)
hold on
axis square
yyaxis left
ylim([0 1.2])
xlabel('Relative day')
ylabel('Normalized correlation score')
%if learn (left side)
xticks([1:7])
xlim([0 8])
xticklabels({'1','2','3','4','5','6','7'})
%dashed 1 reference line
plot([0 8],[1 1],'--','Color',[0.5 0.5 0.5])

lA = plot_error_line(PV_TC_plot_data.ABcorr.ts.st_learn.AB  ,'-',2,[139,0,139]/255);

%plot correlation on right y axis
yyaxis right
ylabel('Performance')
ylim([0 1.2])
yticks([0.2 0.4 0.6 0.8 1])

plot_error_line(performance_mean_sem.st_learn,'-',2,[34,139,34]/255);

set(gca,'FontSize',12)
set(gca,'Linewidth',2)

ax = gca;
%set left axis color
ax.YAxis(1).Color = [139,0,139]/255;
%set right axis color
ax.YAxis(2).Color = [34,139,34]/255;
 
%short term recall
nexttile(t1,2)
hold on
axis square
yyaxis left
ylim([0 1.2])
xlabel('Relative day')
ylabel('Normalized correlation score')
%if learn (left side)
        xticks([1:9])
        xlim([0 10])
        xticklabels({'1','2','3','4','5','6','7','8','9'})
%dashed 1 reference line
plot([0 10],[1 1],'--','Color',[0.5 0.5 0.5])

lA = plot_error_line(PV_TC_plot_data.ABcorr.ts.st_recall.AB  ,'-',2,[139,0,139]/255);

%plot correlation on right y axis
yyaxis right
ylabel('Performance')
ylim([0 1.2])
yticks([0.2 0.4 0.6 0.8 1])

plot_error_line(performance_mean_sem.st_recall,'-',2,[34,139,34]/255);

set(gca,'FontSize',12)
set(gca,'Linewidth',2)

ax = gca;
%set left axis color
ax.YAxis(1).Color = [139,0,139]/255;
%set right axis color
ax.YAxis(2).Color = [34,139,34]/255;

%dashed 1 reference line

plot_error_line(performance_mean_sem.st_recall,'-',2,[34,139,34]/255);

%set axis font/label and font size
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12, ...
    'FontWeight','normal', 'LineWidth', 1.5,'layer','top')

end

