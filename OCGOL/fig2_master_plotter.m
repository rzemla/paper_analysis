function [outputArg1,outputArg2] = fig2_master_plotter(fraction_place_plot_data,si_ts_score_dist_data,....
                                    task_sel_STC_data,centroid_dist_data,AUC_data, tc_corr_sel_data)


%% Unload data from input structure

%fraction of task-specific place cells across animals
%means
frac_tuned_each_mean = fraction_place_plot_data.frac_tuned_each_mean;
%sems
frac_tuned_each_sem = fraction_place_plot_data.frac_tuned_each_sem;

%unload si/ts score distribution data
si = si_ts_score_dist_data.si;
ts = si_ts_score_dist_data.ts;
mean_si_each = si_ts_score_dist_data.mean_si_each;
mean_ts_each = si_ts_score_dist_data.mean_ts_each;

%unload rate map data
%reward A and B black marker lines
rewA_bin_pos = task_sel_STC_data.rewA_bin_pos;
rewB_bin_pos = task_sel_STC_data.rewB_bin_pos;
odor_bin_pos = task_sel_STC_data.odor_bin_pos;

%A sel on A laps
AselA = task_sel_STC_data.AselA;
%A sel on B laps
AselB = task_sel_STC_data.AselB;
%B sel on A laps
BselA = task_sel_STC_data.BselA; 
%B sel on B laps
BselB = task_sel_STC_data.BselB;

%for histogram distribution plot
pooled_A_counts_norm = centroid_dist_data.pooled_A_counts_norm;
pooled_B_counts_norm = centroid_dist_data.pooled_B_counts_norm;
A_zone_start = centroid_dist_data.A_zone_start;
B_zone_start = centroid_dist_data.B_zone_start;

%data for cdfplot
edge_center = centroid_dist_data.edge_center;
meanAsel = centroid_dist_data.meanAsel;
semAsel = centroid_dist_data.semAsel;

meanBsel = centroid_dist_data.meanBsel;
semBsel = centroid_dist_data.semBsel;

%AUC data
grouped_means_run = AUC_data.grouped_means_run;
grouped_sem_run = AUC_data.grouped_sem_run;
grouped_means_norun = AUC_data.grouped_means_norun;
grouped_sem_norun = AUC_data.grouped_sem_norun;
mean_AUC = AUC_data.mean_AUC;

%TC correlaton data boxplot
mean_TC.Asel = tc_corr_sel_data.mean_TC.Asel;
mean_TC.Bsel = tc_corr_sel_data.mean_TC.Bsel;
mean_TC.AB = tc_corr_sel_data.mean_TC.AB;

%% Plotter task-selective distribution

%paper color scheme for A,B, A&B, and neither
color_groups = [65,105,225; 220,20,60; 139, 0, 139; 128 128 128]./255;

%plot bar
%set size of the figure
fig = figure;
fig.Units = 'centimeters';
fig.Position(1) = 8;
fig.Position(2) = 1;
fig.Position(3) = 24;
fig.Position(4) = 18;

%tiledLayout (replaces subplot fxn)
gridSize = [2,3];
t1 = tiledlayout(fig,gridSize(1),gridSize(2),'TileSpacing','compact','Padding','normal','Units','centimeters');

% gridSize = [2,3];
% t2 = tiledlayout(t,gridSize(1),gridSize(2),'TileSpacing','compact','Padding','normal','Units','centimeters'); 
% t2.Layout.Tile = 1;
% t2.Layout.TileSpan = [1,1];

nexttile(t1)
hold on
axis square
title('Spatial Information');
%bar the mean for each group
b = bar(1:4,frac_tuned_each_mean.si,'FaceColor', 'flat');
pause(0.1)
%plot the sem for each mean for each group
for ib = 1:numel(b)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = b(ib).XData + b(ib).XOffset;
    
    errorbar(xData,frac_tuned_each_mean.si',frac_tuned_each_sem.si,'k.','LineWidth',1.5)
end

%set each bar to group color
b(1).CData(1:4,:) =  color_groups;
xticks([1 2 3 4]);
xticklabels({'A','B','A&B', 'Neither'});
ylabel('Fraction of neurons');
xlim([0 5])
ylim([0 0.6])
yticks([0:0.1:0.6])

% Scatterplots of mean scores
%figure('Position',[2295 322 1080 420]);
nexttile(2)
hold on
axis square
xlim([0 0.25])
ylim([0 0.25])
xticks([0 0.1 0.2])
yticks([0 0.1 0.2])
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
xlabel('A trials')
ylabel('B trials')
title('Mean S.I. score')

%plot connecting lines for same animal points
% for ee=1:size(path_dir,2)
%     plot([mean_AUC.A(ee,1),mean_AUC.B(ee,1)],[mean_AUC.A(ee,2),mean_AUC.B(ee,2)] ,'Color',[1 1 1]*0.5,'LineWidth',0.5)
% end

%plot center line
plot([0 10], [0 10],'Color',[1 1 1]*0,'LineStyle', '--','LineWidth',2)

s1 = scatter(mean_si_each.Aonly(:,1),mean_si_each.Aonly(:,2),'filled','MarkerFaceColor',color_groups(1,:));
s2 = scatter(mean_si_each.Bonly(:,1),mean_si_each.Bonly(:,2),'filled','MarkerFaceColor',color_groups(2,:));
s4 = scatter(mean_si_each.N(:,1),mean_si_each.N(:,2),'filled','MarkerFaceColor',color_groups(4,:));
s3 = scatter(mean_si_each.AB(:,1),mean_si_each.AB(:,2),'filled','MarkerFaceColor',color_groups(3,:));
legend([s1 s2 s3 s4],{'A','B','A&B','Neither'},'location','northeast')


nexttile(4)
hold on;
axis square
title('Tuning Specificity');
%bar the mean for each group
b2 = bar(1:4,frac_tuned_each_mean.ts,'FaceColor', 'flat');
pause(0.1)
%plot the sem for each mean for each group
for ib = 1:numel(b2)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = b2(ib).XData + b2(ib).XOffset;
    
    errorbar(xData,frac_tuned_each_mean.ts',frac_tuned_each_sem.ts,'k.','LineWidth',1.5)
end

%set each bar to group color
b2(1).CData(1:4,:) =  color_groups;
xticks([1 2 3 4]);
xticklabels({'A','B','A&B', 'Neither'});
ylabel('Fraction of neurons');
xlim([0 5])
ylim([0 0.6])
yticks([0:0.1:0.6])


nexttile(5) 
hold on
axis square
xlim([0 1])
ylim([0 1])
xticks([0 0.5 1])
yticks([0 0.5 1])
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
xlabel('A trials')
ylabel('B trials')
title('Mean T.S. score')

%plot connecting lines for same animal points
% for ee=1:size(path_dir,2)
%     plot([mean_AUC.A(ee,1),mean_AUC.B(ee,1)],[mean_AUC.A(ee,2),mean_AUC.B(ee,2)] ,'Color',[1 1 1]*0.5,'LineWidth',0.5)
% end

%plot center line
plot([0 10], [0 10],'Color',[1 1 1]*0,'LineStyle', '--','LineWidth',2)

s1 = scatter(mean_ts_each.Aonly(:,1),mean_ts_each.Aonly(:,2),'filled','MarkerFaceColor',color_groups(1,:),'MarkerEdgeColor','none');
s2 = scatter(mean_ts_each.Bonly(:,1),mean_ts_each.Bonly(:,2),'filled','MarkerFaceColor',color_groups(2,:),'MarkerEdgeColor','none');
s4 = scatter(mean_ts_each.N(:,1),mean_ts_each.N(:,2),'filled','MarkerFaceColor',color_groups(4,:),'MarkerEdgeColor','none');
s3 = scatter(mean_ts_each.AB(:,1),mean_ts_each.AB(:,2),'filled','MarkerFaceColor',color_groups(3,:),'MarkerEdgeColor','none');
legend([s1 s2 s3 s4],{'A','B','A&B','Neither'},'location','southeast')


% Plot scatter plots of all neurons for each category
marker_size = 5;

%try scatterplot
%figure('Position',[2208 244 1198 512])
%si
nexttile(3)
hold on
title('S.I. score')
axis square
xlim([-0.01 0.25])
ylim([-0.01 0.25])
xticks([0 0.1 0.2])
yticks([0 0.1 0.2])
xlabel('A laps')
ylabel('B laps')
scatter(si.N_cumul(:,1),si.N_cumul(:,2),marker_size,'MarkerFaceColor',color_groups(4,:),'MarkerEdgeColor',color_groups(4,:))
scatter(si.Aonly_cumul(:,1),si.Aonly_cumul(:,2),marker_size,'MarkerFaceColor',color_groups(1,:),'MarkerEdgeColor',color_groups(1,:))
scatter(si.Bonly_cumul(:,1),si.Bonly_cumul(:,2),marker_size,'MarkerFaceColor',color_groups(2,:),'MarkerEdgeColor',color_groups(2,:))
scatter(si.AB_cumul(:,1),si.AB_cumul(:,2),marker_size,'MarkerFaceColor',color_groups(3,:),'MarkerEdgeColor',color_groups(3,:))
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%center line
plot([-0.01 0.25],[-0.01 0.25],'k--')

%ts
nexttile(6)
hold on
title('T.S. score')
axis square
xlim([0 1.1])
ylim([0 1.1])
xticks([0 0.5 1])
yticks([0 0.5 1])
xlabel('A laps')
ylabel('B laps')
scatter(ts.N_cumul(:,1),ts.N_cumul(:,2),marker_size,'MarkerFaceColor',color_groups(4,:),'MarkerEdgeColor',color_groups(4,:))
scatter(ts.Aonly_cumul(:,1),ts.Aonly_cumul(:,2),marker_size,'MarkerFaceColor',color_groups(1,:),'MarkerEdgeColor',color_groups(1,:))
scatter(ts.Bonly_cumul(:,1),ts.Bonly_cumul(:,2),marker_size,'MarkerFaceColor',color_groups(2,:),'MarkerEdgeColor',color_groups(2,:))
scatter(ts.AB_cumul(:,1),ts.AB_cumul(:,2),marker_size,'MarkerFaceColor',color_groups(3,:),'MarkerEdgeColor',color_groups(3,:))
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%center line
plot([0 1.1],[0 1.1],'k--')


%apply formatting settings across all subplots
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10, ...
'FontWeight','normal', 'LineWidth', 1.5,'layer','top')

%remove whitespace around figure - not sure how this works so far
%set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))


%% Plot color-coded STCs

%define colors maps for rasters
cmap_blue=cbrewer('seq', 'Blues', 32);
cmap_red=cbrewer('seq', 'Reds', 32);

%%%%% A selective neurons - STC on A laps and STC on B laps %%%%%

%set 0 value to white for both blue and red
cmap_blue(1,:) = [1 1 1];
cmap_red(1,:) = [1 1 1];

% gridSize = [5,2];
% t3 = tiledlayout(t,gridSize(1),gridSize(2),'TileSpacing','compact','Padding','normal','Units','centimeters');
% %positioning of rasters with overall master layout
% t3.Layout.Tile = 2;
% %span of raster layout within overall master layout
% t3.Layout.TileSpan = [2,1];

fig2 = figure;
fig2.Units = 'centimeters';
fig2.Position(1) = 8;
fig2.Position(2) = 1;
fig2.Position(3) = 24;
fig2.Position(4) = 36;

%tiledLayout (replaces subplot fxn)
gridSize = [5,2];
t2 = tiledlayout(fig2,gridSize(1),gridSize(2),'TileSpacing','compact','Padding','normal','Units','centimeters');
%nest tiles within master tile

%how many cells for each raster to occupy
raster_size = [2,1];
%input raster for display (A and B neighboring)
%for A selective neurons

%A
nexttile(t2,raster_size);
%A laps
input_matrix = AselA;
%create blank alpha shading matrix where 
imAlpha=ones(size(input_matrix));
imAlpha(isnan(input_matrix))=0;
imagesc(input_matrix,'AlphaData',imAlpha);

hold on
%pbaspect([1 1.5 1])
title('Asel - A laps')

%set background axis color to black
set(gca,'color',0*[1 1 1]);
%set colormap to 
caxis([0 1])
colormap(gca,cmap_blue);
cbar= colorbar;
cbar.Label.String = 'Normalized activity';
cbar.Ticks = [0 0.5 1];
ax1 = gca;
ylabel('Neuron #');
%xlabel('Normalized position');
ax1.XTick = [1 100];
ax1.XTickLabel = {'0','1'};

%plot reward A and reward B zones
%A
plot([rewA_bin_pos rewA_bin_pos],[1 size(input_matrix,1)],'k--')
%B
plot([rewB_bin_pos rewB_bin_pos],[1 size(input_matrix,1)],'k--')
%odor pos
plot([odor_bin_pos odor_bin_pos],[1 size(input_matrix,1)],'k--')

%make ticks invisible
set(ax1, 'TickLength', [0 0]);

set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%B laps as input
input_matrix = AselB;

%B 
nexttile(raster_size);
%create blank alpha shading matrix where 
imAlpha=ones(size(input_matrix));
imAlpha(isnan(input_matrix))=0;
imagesc(input_matrix,'AlphaData',imAlpha);
hold on
%pbaspect([1 1.5 1])
title('Asel - B laps')
%set background axis color to black
set(gca,'color',0*[1 1 1]);
%set colormap to 
caxis([0 1])
colormap(gca,cmap_red);
cbar2 = colorbar;
cbar2.Label.String = 'Normalized activity';
cbar2.Ticks = [0 0.5 1];
ax2 = gca;
%ylabel('Neuron #');
%xlabel('Normalized position');
ax2.XTick = [1 100];
ax2.XTickLabel = {'0','1'};
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%plot reward A and reward B zones
%A
plot([rewA_bin_pos rewA_bin_pos],[1 size(input_matrix,1)],'k--')
%B
plot([rewB_bin_pos rewB_bin_pos],[1 size(input_matrix,1)],'k--')
%odor pos
plot([odor_bin_pos odor_bin_pos],[1 size(input_matrix,1)],'k--')

%make ticks invisible
set(ax2, 'TickLength', [0 0]);


%%%%% B selective neurons - STC on A laps and STC on B laps %%%%%

% %set 0 value to white for both blue and red
% cmap_blue(1,:) = [1 1 1];
% cmap_red(1,:) = [1 1 1];

%input raster for display (A and B neighboring)
%for A selective neurons
%figure('Position',[2182 356 934 559])
%A
nexttile(raster_size);
%A laps
input_matrix = BselA;

%create blank alpha shading matrix where 
imAlpha=ones(size(input_matrix));
imAlpha(isnan(input_matrix))=0;
imagesc(input_matrix,'AlphaData',imAlpha);
hold on
%pbaspect([1 1.5 1])
title('Bsel - A laps')
%set background axis color to black
set(gca,'color',0*[1 1 1]);
%set colormap to 
caxis([0 1])
colormap(gca,cmap_blue);
if 0
    cbar= colorbar;
    cbar.Label.String = 'Normalized activity';
    cbar.Ticks = [0 0.5 1];
end
ax1 = gca;
ylabel('Neuron #');
xlabel('Normalized position');
ax1.XTick = [1 100];
ax1.XTickLabel = {'0','1'};

%plot reward A and reward B zones
%A
plot([rewA_bin_pos rewA_bin_pos],[1 size(input_matrix,1)],'k--')
%B
plot([rewB_bin_pos rewB_bin_pos],[1 size(input_matrix,1)],'k--')
%odor pos
plot([odor_bin_pos odor_bin_pos],[1 size(input_matrix,1)],'k--')

%make ticks invisible
set(ax1, 'TickLength', [0 0]);

set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%B laps as input
input_matrix = BselB;

%B 
nexttile(raster_size);
%create blank alpha shading matrix where 
imAlpha=ones(size(input_matrix));
imAlpha(isnan(input_matrix))=0;
imagesc(input_matrix,'AlphaData',imAlpha);
hold on
%pbaspect([1 1.5 1])
title('Bsel - B laps')
%set background axis color to black
set(gca,'color',0*[1 1 1]);
%set colormap to 
caxis([0 1])
colormap(gca,cmap_red);

if 0
    cbar2 = colorbar;
    cbar2.Label.String = 'Normalized activity';
    cbar2.Ticks = [0 0.5 1];
end

ax2 = gca;
ylabel('Neuron #');
xlabel('Normalized position');
ax2.XTick = [1 100];
ax2.XTickLabel = {'0','1'};
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%plot reward A and reward B zones
%A
plot([rewA_bin_pos rewA_bin_pos],[1 size(input_matrix,1)],'k--')
%B
plot([rewB_bin_pos rewB_bin_pos],[1 size(input_matrix,1)],'k--')
%odor pos
plot([odor_bin_pos odor_bin_pos],[1 size(input_matrix,1)],'k--')

%make ticks invisible
set(ax2, 'TickLength', [0 0]);

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',14, ...
'FontWeight','normal', 'LineWidth', 1.5,'layer','top')

%% Plot Pooled Normalized histogram for centroid distributions
%return the colormaps used in paper and apply to bars
paper_cmap = return_paper_colormap;

% gridSize = [2,2];
% t4 = tiledlayout(t,gridSize(1),gridSize(2),'TileSpacing','normal','Padding','normal','Units','centimeters');
% %positioning of rasters with overall master layout
% t4.Layout.Tile = 26;
% %span of raster layout within overall master layout
% t4.Layout.TileSpan = [2,4];

raster_size = [1,1];

nexttile(t2,raster_size)
hold on
%axis square
xlim([0,25.5])
ylim([0 0.11])
yticks([0 0.05 0.1])
%mean by bin
b1 = bar(1:25,pooled_A_counts_norm,1);
%sem by bin
%errorbar(1:25,mean_norm_A,sem_norm_A,'LineStyle','none','Color','k')

xticks([0 25])
xticklabels({'0','1'})
xlabel('Normalized position')
ylabel('Normalized density')
b1.FaceColor = paper_cmap(1,:);
b1.FaceAlpha = 1;

%start zones
plot([B_zone_start B_zone_start]./4, [0 0.15],'LineStyle','--', 'LineWidth',2,'Color',[220,20,60]./255)
plot([A_zone_start A_zone_start]./4, [0 0.15],'LineStyle','--', 'LineWidth',2,'Color',[65,105,225]./255)

%add odor marker
plot([10 10]./4, [0 0.15],'LineStyle','--', 'LineWidth',2,'Color',paper_cmap(end,:))

set(gca,'FontSize',16)
set(gca,'LineWidth',2)
%end zones
% plot([B_zone_end B_zone_end]./4, [0 0.15],'LineStyle','--', 'LineWidth',2,'Color',[220,20,60]./255)
% plot([A_zone_end-1 A_zone_end-1]./4, [0 0.15],'LineStyle','--', 'LineWidth',2,'Color',[65,105,225]./255)

nexttile(t2, raster_size);
hold on
%axis square
xlim([0,25.5])
ylim([0 0.11])
yticks([0 0.05 0.1])
%mean by bin
b1 = bar(1:25,pooled_B_counts_norm,1);
%sem by bin
%errorbar(1:25,mean_norm_B,sem_norm_B,'LineStyle','none','Color','k')

xticks([0 25])
xticklabels({'0','1'})
xlabel('Normalized position')
ylabel('Normalized density')
b1.FaceColor = paper_cmap(2,:);
b1.FaceAlpha = 1;

%start zones
plot([B_zone_start B_zone_start]./4, [0 0.15],'LineStyle','--', 'LineWidth',2,'Color',[220,20,60]./255)
plot([A_zone_start A_zone_start]./4, [0 0.15],'LineStyle','--', 'LineWidth',2,'Color',[65,105,225]./255)

%add odor marker
plot([10 10]./4, [0 0.15],'LineStyle','--', 'LineWidth',2,'Color',paper_cmap(end,:))

set(gca,'FontSize',16)
set(gca,'LineWidth',2)
%end zones
% plot([B_zone_end B_zone_end]./4, [0 0.15],'LineStyle','--', 'LineWidth',2,'Color',[220,20,60]./255)
% plot([A_zone_end-1 A_zone_end-1]./4, [0 0.15],'LineStyle','--', 'LineWidth',2,'Color',[65,105,225]./255)
%[N, edges] = histcounts(1:100,25);
%bin centers for 25 bins
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',14, ...
'FontWeight','normal', 'LineWidth', 1.5,'layer','top')


%% Plot cdf for field difference
%blue, red, light blue, light red      
%color_mat = [65,105,225; 220,20,60; 135,206,250; 240,128,128]./255;
fig3 = figure;
fig3.Units = 'centimeters';
fig3.Position(1) = 8;
fig3.Position(2) = 1;
fig3.Position(3) = 24;
fig3.Position(4) = 9;

%tiledLayout (replaces subplot fxn)
gridSize = [1,2];
t3 = tiledlayout(fig3,gridSize(1),gridSize(2),'TileSpacing','compact','Padding','normal','Units','centimeters');

%cdf plot of place field distribution for A vs. B selective place cells
nexttile(t3,1);
lw = 2;
hold on
axis square
xlabel('Normalized position')
ylabel('Fraction of place cells')
xticks([0,50,100]);
xticklabels({'0','0.5','1'});
yticks([0 0.2 0.4 0.6 0.8 1])

e1 = errorbar(edge_center,meanAsel,semAsel,'LineWidth',lw,'Color',paper_cmap(1,:));
e2 = errorbar(edge_center,meanBsel,semBsel,'LineWidth',lw,'Color',paper_cmap(2,:));

legend([e1 e2],{'A','B'},'Location','southeast')

% TC correlation among A-sel, B-sel, and A&B place cells
%xlab={'1','2','3'};
col=[220,20,60, 255; 65,105,225, 255; 139, 0, 139, 255];  
col=col/255;

%create merge cell for boxplot fxn (using animal means vs. pooled neurons)
merge_new = cell(1,3);
merge_new{1} =  mean_TC.Asel;
merge_new{2} = mean_TC.Bsel;
merge_new{3} = mean_TC.AB;
%convert into linear vector with grouping variable
mean_TC_vec = [merge_new{1}'; merge_new{2}'; merge_new{3}'];
mean_TC_grp = [repmat(["A"],numel(merge_new{1}),1);...
                repmat(["B"],numel(merge_new{2}),1);...
                repmat(["A&B"],numel(merge_new{3}),1)];
       
nexttile(t3,2)
hold on
%title('Change colors manually in illustrator')
%order of categorical variables to plots
box_order = {'A','B','A&B'};
b = boxchart(mean_TC_vec','GroupByColor',categorical(mean_TC_grp',box_order));
%box colors
for ii=1:3
    b(ii).BoxFaceColor =paper_cmap(ii,:);
    b(ii).BoxFaceAlpha = 1;
    b(ii).LineWidth = 1.5;
end

%make black box outline and median line black with overlay boxplot
b2 = boxchart(mean_TC_vec','GroupByColor',categorical(mean_TC_grp',box_order));
%box colors
for ii=1:3
    b2(ii).BoxFaceColor = 'k';
    b2(ii).BoxFaceAlpha = 0.0;
    b2(ii).LineWidth = 1.5;
end

ylabel('Correlation score');

legend(b,'Location','southeast')
xticklabels('A vs. B lap \newline correlation');

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',14, ...
'FontWeight','normal', 'LineWidth', 1.5,'layer','top')



%% Plot bar plots - RUN and NO RUN
fig0 = figure;
fig0.Units = 'centimeters';
fig0.Position(1) = 8;
fig0.Position(2) = 1;
fig0.Position(3) = 24;
fig0.Position(4) = 12;

%tiledLayout (replaces subplot fxn)
gridSize = [2,4];
t0 = tiledlayout(fig0,gridSize(1),gridSize(2),'TileSpacing','compact','Padding','compact','Units','centimeters');
%nest tiles within master tile

% Mean plot
nexttile(t0,1,[2,2])
hold on
axis square
xlim([0 10])
ylim([0 10])
xticks(0:2:10)
yticks(0:2:10)
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
xlabel('AUC/min - A trials')
ylabel('AUC/min - B trials')
title('A, B selective RUN')

%plot connecting lines for same animal points
% for ee=1:size(path_dir,2)
%     plot([ mean_AUC.run.Bsel(ee,1), mean_AUC.run.AB(ee,1), mean_AUC.run.Asel(ee,1), ],...
%         [mean_AUC.run.Bsel(ee,2), mean_AUC.run.AB(ee,2),mean_AUC.run.Asel(ee,2), ] ,'Color',[1 1 1]*0.5,'LineWidth',0.5)
% end

%A selective
for ee=1:11%size(path_dir,2)
    s1 = scatter(mean_AUC.run.Asel(ee,1),mean_AUC.run.Asel(ee,2),'filled','MarkerFaceColor',paper_cmap(1,:));
end

%B selective
for ee=1:11%size(path_dir,2)
    s2 = scatter(mean_AUC.run.Bsel(ee,1),mean_AUC.run.Bsel(ee,2),'filled','MarkerFaceColor',paper_cmap(2,:));
end

for ee=1:11%size(path_dir,2)
    s3 = scatter(mean_AUC.run.AB(ee,1),mean_AUC.run.AB(ee,2),'filled','MarkerFaceColor',paper_cmap(3,:));
end

%plot center line
plot([0 10], [0 10],'Color',[1 1 1]*0,'LineStyle', '--','LineWidth',2)

legend([s1 s2 s3],{'A','B','A&B'},'location','southeast')

clear xData

nexttile(t0,3,[1,2])
hold on;
%axis square
title('Run');
%bar the mean for each group
b = bar(1:3,grouped_means_run,'FaceColor', 'flat');
pause(0.1)
ylabel('AUC/min') 
xlim([0.5 3.5])
ylim([0 10])
%plot the sem for each mean for each group
for ib = 1:numel(b)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    if ib ==1
        xData(1,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==2
        xData(2,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==3
        xData(3,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==4
        xData(4,:) = b(ib).XData + b(ib).XOffset;
    end
    errorbar(xData(ib,:),grouped_means_run(:,ib)',grouped_sem_run(:,ib),'k.','LineWidth',1)
end

%set A group bars to blue
b(1).CData(1:3,:) =  repmat(paper_cmap(1,:),3,1);
%set B group bars to red
b(2).CData(1:3,:) =  repmat(paper_cmap(2,:),3,1);
%set B group bars to red
%b(3).CData(1:3,:) =  repmat(color_mat(3,:),3,1);
%set B group bars to red
%b(4).CData(1:3,:) =  repmat(color_mat(4,:),3,1);

xticks([1 2 3]);
xticklabels({'A sel.','B sel.','A&B'});

legend('A laps','B laps')

set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)

nexttile(t0,7,[1,2])
hold on;
%axis square
title('No run');
%bar the mean for each group
b2 = bar(1:3,grouped_means_norun,'FaceColor', 'flat');
pause(0.1)
ylabel('AUC/min') 
xlim([0.5 3.5])
ylim([0 2])
%plot the sem for each mean for each group
for ib = 1:numel(b2)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    if ib ==1
        xData(1,:) = b2(ib).XData + b2(ib).XOffset;
    elseif ib ==2
        xData(2,:) = b2(ib).XData + b2(ib).XOffset;
    elseif ib ==3
        xData(3,:) = b2(ib).XData + b2(ib).XOffset;
    elseif ib ==4
        xData(4,:) = b2(ib).XData + b2(ib).XOffset;
    end
    errorbar(xData(ib,:),grouped_means_norun(:,ib)',grouped_sem_norun(:,ib),'k.','LineWidth',1)
end

%set A group bars to blue
b2(1).CData(1:3,:) =  repmat(paper_cmap(1,:),3,1);
%set B group bars to red
b2(2).CData(1:3,:) =  repmat(paper_cmap(2,:),3,1);
%set B group bars to red
%b2(3).CData(1:3,:) =  repmat(color_mat(3,:),3,1);
%set B group bars to red
%b(4).CData(1:3,:) =  repmat(color_mat(4,:),3,1);

xticks([1 2 3]);
xticklabels({'A sel.','B sel.','A&B'});
ylabel('AUC/min');
legend('A laps','B laps')

% set(gca,'FontSize',16)
% set(gca,'LineWidth',1.5)

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',14, ...
'FontWeight','normal', 'LineWidth', 1.5,'layer','top')

                

end

