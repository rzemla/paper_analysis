function [outputArg1,outputArg2] = fig3_master_plotter(remap_rate_maps,activity_remap,frac_remapping,remap_prop_figs)
%% Unload the data

%rate maps (common, rate, global, partial, unclassified
first_col_STCs = remap_rate_maps.first_col_STCs;
second_col_STCs = remap_rate_maps.second_col_STCs;
%activity remapping dF/F over time (5s)
mean_activityA = activity_remap.mean_diff_sort_timeA;
mean_activityB = activity_remap.mean_diff_sort_timeB;

%merge into one cell and add activity remapping map here
%common(1) activity/rate(2), activity dF/F(3), global(4), partial (5), unclass (6), 
rasters = [first_col_STCs(1:2); {[mean_activityA, mean_activityB]};...
    first_col_STCs(3); second_col_STCs(1:2)];

%% Paper colors
paper_cmap = return_paper_colormap;

%% Figure 3c Raster rate plots
%still have to make separate titledlayout with colorbar(can't integrate to
%show in only one tile)

fig = figure; % event based STC;
fig.Units = 'centimeters';
fig.Position(1) = 8;
fig.Position(2) = 1;
fig.Position(3) = 20;
fig.Position(4) = 24;

%master layout
gridSize = [3,2];
t1 = tiledlayout(fig,gridSize(1),gridSize(2),'TileSpacing','normal','Padding','compact','Units','centimeters');

%common subplot
s1 = tiledlayout(t1,1,2,'TileSpacing','compact','Padding','compact','Units','centimeters');
s1.Layout.Tile = 1;
s1.Layout.TileSpan = [1,1];
%subplot title
set_subplot_title(s1,'Common', 'Arial', 16,'bold')

%A
nexttile(s1,1)
raster_plot(rasters{1}(:,1:100),'blue')
set_title('A trials',12,'Arial', 'italic',paper_cmap(1,:),0,'west')
ylabel('Neuron #');


%B
nexttile(s1,2)
raster_plot(rasters{1}(:,101:200),'red')
set_title('B trials',12,'Arial', 'italic',paper_cmap(2,:),0,'west')
yticks([]);

%cba.Layout.Location = 'eastoutside';
%set_bar(s1,'west',paper_cmap(2,:))


%global subplot
s2 = tiledlayout(t1,1,2,'TileSpacing','compact','Padding','compact','Units','centimeters');
s2.Layout.Tile = 2;
s2.Layout.TileSpan = [1,1];
%subplot title
set_subplot_title(s2,'Global', 'Arial', 16,'bold')


%A
nexttile(s2)
raster_plot(rasters{4}(:,1:100),'blue')
set_title('A trials',12,'Arial', 'italic',paper_cmap(1,:),'0','west')

%B
nexttile(s2)
raster_plot(rasters{4}(:,101:200),'red')
set_title('B trials',12,'Arial', 'italic',paper_cmap(2,:),0,'west')
yticks([]);


%activity subplot
s3 = tiledlayout(t1,1,2,'TileSpacing','compact','Padding','compact','Units','centimeters');
s3.Layout.Tile = 3;
s3.Layout.TileSpan = [1,1];
set_subplot_title(s3,'Activity', 'Arial', 16,'bold')
%shared x-axis label
xlabel(s3,'Normalized Position','FontName', 'Arial', 'FontSize', 14)

%A
nexttile(s3)
raster_plot(rasters{2}(:,1:100),'blue')
set_title('A trials',12,'Arial', 'italic',paper_cmap(1,:),0,'west')
ylabel('Neuron #');
%B
nexttile(s3)
raster_plot(rasters{2}(:,101:200),'red')
set_title('B trials',12,'Arial', 'italic',paper_cmap(2,:),0,'west')
yticks([]);

%partial subplot
s4 = tiledlayout(t1,1,2,'TileSpacing','compact','Padding','compact','Units','centimeters');
s4.Layout.Tile = 4;
s4.Layout.TileSpan = [1,1];
set_subplot_title(s4,'Partial', 'Arial', 16,'bold')
%A
nexttile(s4)
raster_plot(rasters{5}(:,1:100),'blue')
set_title('A trials',12,'Arial', 'italic',paper_cmap(1,:),0,'west')
%B
nexttile(s4)
raster_plot(rasters{5}(:,101:200),'red')
set_title('B trials',12,'Arial', 'italic',paper_cmap(2,:),0,'west')
yticks([]);

%dF/F activity subplot
s5 = tiledlayout(t1,1,2,'TileSpacing','compact','Padding','compact','Units','centimeters');
s5.Layout.Tile = 5;
s5.Layout.TileSpan = [1,1];
set_subplot_title(s5,'Activity', 'Arial', 16,'bold')
%subtitle properties for this sublayout
s5.Subtitle.String ='\Delta Ca^{\it2+} event magnitude';
s5.Subtitle.FontAngle = 'italic';
xlabel(s5,'Time since event onset [s]','FontName', 'Arial', 'FontSize', 14)
ylabel(s5,'Neuron #');
%A
nexttile(s5)
raster_plot(rasters{3}(:,1:200),'jet')
set_title('A trials',12,'Arial', 'italic',paper_cmap(1,:),0,'west')
%B
nexttile(s5)
raster_plot(rasters{3}(:,201:400),'jet')
set_title('B trials',12,'Arial', 'italic',paper_cmap(2,:),0,'west')
yticks([]);

%unclassified subplot
s6 = tiledlayout(t1,1,2,'TileSpacing','compact','Padding','compact','Units','centimeters');
s6.Layout.Tile = 6;
s6.Layout.TileSpan = [1,1];
set_subplot_title(s6,'Unclassified', 'Arial', 16,'bold')
%shared x-axis label
xlabel(s6,'Normalized Position','FontName', 'Arial', 'FontSize', 14)

%A
nexttile(s6)
raster_plot(rasters{6}(:,1:100),'blue')
set_title('A trials',12,'Arial', 'italic',paper_cmap(1,:),0,'west')
%B
nexttile(s6)
raster_plot(rasters{6}(:,101:200),'red')
set_title('B trials',12,'Arial', 'italic',paper_cmap(2,:),0,'west')
yticks([]);


%set axis font/label and font size
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12, ...
'FontWeight','normal', 'LineWidth', 1.5,'layer','top')


% Set subplot tile title
    function set_title(t_name,t_size,t_font, t_angle,t_color,set_bar,bar_where)
        t = title(t_name);
        t.FontSize = t_size;
        t.FontName = t_font;
        t.FontAngle = t_angle;
        t.Color = t_color;
        if set_bar == 1
            cb = colorbar();
            cb.Layout.Tile = bar_where;
        end
    end
%set subplotbar
%     function set_bar(tx,bar_where,color_map)
%             colormap(tx,color_map)
%             cb = colorbar(tx);
%             cb.Layout.Tile = bar_where;
%     end

    function set_subplot_title(sub,t_string, t_name, t_size,t_weight)
        sub.Title.String = t_string;
        sub.Title.FontName = t_name;
        sub.Title.FontSize = t_size;
        sub.Title.FontWeight = t_weight;
    end

%% Remapping properies plots in Fig. 3

%modifed to used merge data classes

%unload the data
%frac remap bar plot (3d)
class_names = frac_remapping.merge.class_names;
class_mean = frac_remapping.merge.class_mean ;
class_sem = frac_remapping.merge.class_sem;

fig = figure; % event based STC;
fig.Units = 'centimeters';
fig.Position(1) = 8;
fig.Position(2) = 1;
fig.Position(3) = 20;
fig.Position(4) = 24;

%master layout
gridSize = [3,2];
t1 = tiledlayout(fig,gridSize(1),gridSize(2),'TileSpacing','normal','Padding','normal','Units','centimeters');

%3d fractional distribution of each class of remapping neurons
nexttile(t1, 1)
hold on

bar(class_mean,'FaceColor',[139, 0, 139]/255)
%add significance bars
xticks([1:4])
xticklabels({'Common','Rate','Global','Other'})
%sigstar({[1,2], [1,3],[1,4]})
xtickangle(45)
errorbar([1:4],class_mean,class_sem,'k.')
ylim([0 0.6])
ylabel('Fraction of remapping neurons')

%unload data
common_bins_combined = remap_prop_figs.common_dist.common_bins_combined;
B_zone_end = remap_prop_figs.common_dist.B_zone_end; 
A_zone_end = remap_prop_figs.common_dist.A_zone_end;  

%3e common place field distribution (add manual zone labels)
nexttile(t1,2)
hold on
ylim([0 0.2])
yticks([0 0.05 0.1 0.15 0.2])
h1 = histogram(common_bins_combined(1,:),0:10:100,'Normalization','probability');
xticks([0 100])
xticklabels({'0','1'})
xlabel('Normalized position')
ylabel('Normalized density')
h1.FaceColor = [139 0 139]./255;
h1.FaceAlpha = 1;

%place an annotation for each zone (later)

%annotation('textbox',[1 1 1 1],'String','Zone I','FitBoxToText','on','Units','centimeters');

plot([B_zone_end B_zone_end], [0 0.2],'LineStyle','--', 'LineWidth',2,'Color',[220,20,60]./255)
plot([A_zone_end A_zone_end], [0 0.2],'LineStyle','--', 'LineWidth',2,'Color',[65,105,225]./255)
%3x2 subplot for partial 3h

%fig 3f shift of place fields in reward zones between A and B trials
%unload data
zoneI_diff = remap_prop_figs.global_remap_shift.zoneI_diff;
zoneII_diff = remap_prop_figs.global_remap_shift.zoneII_diff;
zoneIII_diff = remap_prop_figs.global_remap_shift.zoneIII_diff;
zone_med = remap_prop_figs.global_remap_shift.zone_med;
zone_sem = remap_prop_figs.global_remap_shift.zone_sem;

nexttile(t1,3)
hold on
xlim([0 4])
ylim([-0.5 0.6])
xticks([1 2 3])
yticks([-0.4, -0.2, 0, 0.2, 0.4])
xticklabels({'Zone I','Zone II','Zone III'})
xlabel('\newline Relative to A place field');
ylabel({'Fraction of B fields before A fields', '- expected fraction'});
dot_size = 14;
scatter(ones(1,size(zoneI_diff,2)),zoneI_diff,dot_size,'filled','MarkerFaceColor',[0.5 0.5 0.5])
scatter(2*ones(1,size(zoneII_diff,2)),zoneII_diff,dot_size,'filled','MarkerFaceColor',[0.5 0.5 0.5])
scatter(3*ones(1,size(zoneIII_diff,2)),zoneIII_diff,dot_size,'filled','MarkerFaceColor',[0.5 0.5 0.5])
bar([1 2 3],zone_med,'FaceColor',[139, 0, 139]/255,'EdgeColor',[0 0 0]/255,'FaceAlpha',0.3,'BarWidth', 0.6)
errorbar([1 2 3],zone_med,zone_sem,'LineStyle','none','Color', [0 0 0])
%0 dash line
%plot([0 4],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',1)

%fig 3g - global remapping between zones
%unload data
frac_zones = remap_prop_figs.global_zones.frac_zones;
mean_fraction_zone_split = remap_prop_figs.global_zones.mean_fraction_zone_split;
sem_fraction_zone_split = remap_prop_figs.global_zones.sem_fraction_zone_split;

nexttile(t1,4)
hold on
bar([1 2 3 4.5 5.5 7 8 9.5 10.5], mean(frac_zones,1),'FaceColor',[139, 0, 139]/255)
xticks([1 2 3 4.5 5.5 7 8 9.5 10.5])
errorbar([1 2 3 4.5 5.5 7 8 9.5 10.5],mean_fraction_zone_split,sem_fraction_zone_split,'k.')
xticklabels({'AI,BI','AII, BII','AIII, BIII','AI, BII','AII, BI',...
    'AII, BIII','AIII, BII','AI, BIII','AIII, BI'})
xtickangle(45)
%zone I,II,III switches
%sigstar({[1,2], [1,3], [2,3]})
%I vs. II switch, II vs. III switch, I vs. III
%sigstar({[4.5 5.5],[7 8],[9.5 10.5]})
yticks([0 0.1 0.2 0.3])
ylabel('Fraction of neurons');

%3h partial place field remapping (main) + inset(histogram dists)
%unload data
partial_A_common_far_input = remap_prop_figs.partial_remap.partial_A_common_far_input;
partial_B_common_far_input = remap_prop_figs.partial_remap.partial_B_common_far_input;
B_zone_end = remap_prop_figs.partial_remap.B_zone_end;
A_zone_end = remap_prop_figs.partial_remap.A_zone_end;

%cdf here
nexttile(t1,5)
hold on
%title('Distribution of partially remapping fields')
[fA,xA] = ecdf(partial_A_common_far_input(:,2));
[fB,xB] = ecdf(partial_B_common_far_input(:,2));
ap = stairs(xA,fA,'LineWidth',2,'Color',[65,105,225]/255);
bp = stairs(xB,fB,'LineWidth',2,'Color',[220,20,60]/255);
xlabel(gca,'Normalized positon')
ylabel('Cumulative probability')

xticks([0 100])
xticklabels({'0','1'})
%set(gca,'LineWidth',1.5)
%set(gca,'FontSize',16)
legend([ap bp],{'A','B'},'Location','northwest','AutoUpdate','off')

plot([B_zone_end B_zone_end], [0 1],'--','Color',[220,20,60]/255, 'LineWidth',2)
plot([A_zone_end A_zone_end], [0 1],'--','Color',[65,105,225]/255, 'LineWidth',2)

%make subplot for partial field position distribution (3h ctd)
%unload data
partial_A_common_far_input = remap_prop_figs.partial_remap.insets.partial_A_common_far_input;
partial_B_common_far_input = remap_prop_figs.partial_remap.insets.partial_B_common_far_input;
B_zone_end = remap_prop_figs.partial_remap.insets.B_zone_end;
A_zone_end = remap_prop_figs.partial_remap.insets.A_zone_end;

s2 = tiledlayout(t1,2,2,'TileSpacing','normal','Padding','normal','Units','centimeters');
s2.Layout.Tile = 6;
s2.Layout.TileSpan = [1,1];

nexttile(s2,1,[1,1])
hold on;
ylabel(s2,'Norm. density')
ylim([0 0.2])
yticks([0 0.1 0.2])
xticks([0 50 100])
xticklabels({'0','0.5','1'})
%ylabel('Norm. \newline density')
%set(gca,'FontSize',16)
histogram(partial_A_common_far_input(:,2),[0:10:100],'Normalization','probability','FaceColor', [65,105,225]/255,'FaceAlpha',1)

plot([B_zone_end B_zone_end], [0 0.2],'--','Color',[220,20,60]/255, 'LineWidth',2)
plot([A_zone_end A_zone_end], [0 0.2],'--','Color',[65,105,225]/255, 'LineWidth',2)

nexttile(s2,3,[1,1])
hold on
ylim([0 0.2])
yticks([0 0.1 0.2])
xticks([0 50 100])
xticklabels({'0','0.5','1'})

xlabel({'Normalized','position'})
%set(gca,'FontSize',16)
histogram(partial_B_common_far_input(:,2),[0:10:100],'Normalization','probability','FaceColor',[220,20,60]/255,'FaceAlpha',1)

plot([B_zone_end B_zone_end], [0 0.2],'--','Color',[220,20,60]/255, 'LineWidth',2)
plot([A_zone_end A_zone_end], [0 0.2],'--','Color',[65,105,225]/255, 'LineWidth',2)

%set axis font/label and font size
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12, ...
'FontWeight','normal', 'LineWidth', 1.5,'layer','top')


%% Export figs to folder (future)

%export rasters fpr first colum
%mkdir(fullfile(path_dir{1},'example_STCs_Fig3D'))
%disp(['Saving raster: 1'])
%export_fig(f ,fullfile(path_dir{1},'example_STCs_Fig3D',[num2str(1),'_300.png']),'-r300')


end

