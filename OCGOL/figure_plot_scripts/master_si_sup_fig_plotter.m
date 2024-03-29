function [outputArg1,outputArg2] = master_si_sup_fig_plotter(performance_mean_sem, frac_tuned_across_days,PV_TC_plot_data,centroid_data)


%% Plot master figure

fig = figure;
fig.Units = 'centimeters';
fig.Position(1) = 7;
fig.Position(2) = 0;
fig.Position(3) = 42;
fig.Position(4) = 24;

%master layout
gridSize = [2,4];
t1 = tiledlayout(fig,gridSize(1),gridSize(2),'TileSpacing','compact','Padding','compact','Units','centimeters');

%unload data for fraction tuned across time
mean_learning_si_day = frac_tuned_across_days.mean_learning_si_day;
mean_recall_si = frac_tuned_across_days.mean_recall_si;

%color theme
bar_colorset = [65,105,225; 220,20,60; 139, 0, 139; 128 128 128]./255;

%fraction tuned
nexttile(t1,1)
hold on
axis square
ylim([0 1.2])
xlabel('Day')
ylabel('Fraction tuned')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
%title('Learning - day - T.S.')
sh1 = bar(1:7,mean_learning_si_day(1:7,:),'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh1,2)
    sh1(ii).FaceColor = 'flat';
    sh1(ii).CData = bar_colorset(ii,:);
end

nexttile(t1,2)
hold on
axis square
ylim([0 1.2])
xlabel('Day')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
xticklabels({'1','2','3','6','7','8','9'})
%title('Recall - day - T.S.')
sh2 = bar(mean_recall_si,'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh2,2)
    sh2(ii).FaceColor = 'flat';
    sh2(ii).CData = bar_colorset(ii,:);
end

legend({'A','B','A&B','Neither'},'Location','northoutside','Orientation','horizontal','NumColumns',2)

%unload data
TC_mean_sem = PV_TC_plot_data.TC_mean_sem;

%TC correlation
nexttile(t1,3)
hold on
axis square
xlabel('Relative Day to D1')
ylabel('Corr. score')
xticks([1:9])
xticklabels({'0','1','2','3','4','5','6','7','8'})
%title('TC TS correlation - raw')
%learn
lA = plot_error_line(TC_mean_sem.si.st_learn.A,'--',2,[65,105,225]/255);
lB = plot_error_line(TC_mean_sem.si.st_learn.B,'--',2,[220,20,60]/255);

%recall
rA = plot_error_line(TC_mean_sem.si.st_recall.d4_d5.A,'-',2,[65,105,225]/255);
rB = plot_error_line(TC_mean_sem.si.st_recall.d4_d5.B,'-',2,[220,20,60]/255);

legend([lA,lB,rA,rB],{'Learning A','Learning B','Recall A','Recall B'},'location','northeast')

%neighboring day correlation
nexttile(t1,4)
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
%title('TC correlation - SI neighbor')
%learn
lA = plot_error_line(TC_mean_sem.neighbor.si.st_learn.A(:,1:6),'--',2,[65,105,225]/255);
lB = plot_error_line(TC_mean_sem.neighbor.si.st_learn.B(:,1:6),'--',2,[220,20,60]/255);
%recall
rA = plot_error_line(TC_mean_sem.neighbor.si.st_recall.A,'-',2,[65,105,225]/255);
rB = plot_error_line(TC_mean_sem.neighbor.si.st_recall.B,'-',2,[220,20,60]/255);

set(gca,'FontSize',16)
set(gca,'Linewidth',2)

legend([lA,lB,rA,rB],{'Learning A','Learning B','Recall A','Recall B'},'location','southeast')


%unload centroid data
%A
st_learn.animal.si.A = centroid_data.st_learn.si.A;
st_recall.animal.si.A = centroid_data.st_recall.si.A;
%B
st_learn.animal.si.B = centroid_data.st_learn.si.B;
st_recall.animal.si.B = centroid_data.st_recall.si.B;

%A and B lap centroid difference
cm2rad = (2*pi)./196;
%rad ticks that correspond to 0 25 50 cm
rad_ticks = [0 10 20 30 40].*cm2rad;

nexttile(t1,5)
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
lA = plot_error_line(st_learn.animal.si.A(:,2:end),'--',2,[65,105,225]./255);
%recall
rA = plot_error_line(st_recall.animal.si.A(:,2:end),'-',2,[65,105,225]./255);

legend([lA,rA],{'Learning A','Recall A'},'location','southeast')

nexttile(t1,6)
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
lB = plot_error_line(st_learn.animal.si.B(:,2:end),'--',2,[220,20,60]./255);
%recalll
rB = plot_error_line(st_recall.animal.si.B(:,2:end),'-',2,[220,20,60]./255);

legend([lB,rB],{'Learning B','Recall B'},'location','southeast')

%A&B tuned A lap vs. B lap TC correlation across time


PV_TC_plot_data.ABcorr.si.st_recall.AB  

nexttile(t1,7)
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

lA = plot_error_line(PV_TC_plot_data.ABcorr.si.st_learn.AB  ,'-',2,[139,0,139]/255);

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
nexttile(t1,8)
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

lA = plot_error_line(PV_TC_plot_data.ABcorr.si.st_recall.AB  ,'-',2,[139,0,139]/255);

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

