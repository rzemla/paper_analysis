function [outputArg1,outputArg2] = master_lt_figure_plotter(frac_tuned_across_days)



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

%fraction TS


%performance


%PV

%TC SI

%TC TS

%PV neighbor


%TC neighbor SI

%TC neighbor TS


%A&B correlation SI

%A&B correlation TS


end

