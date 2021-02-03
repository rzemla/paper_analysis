function [outputArg1,outputArg2] = fig3_master_plotter(remap_rate_maps,activity_remap)
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


%export rasters fpr first colum
%mkdir(fullfile(path_dir{1},'example_STCs_Fig3D'))
%disp(['Saving raster: 1'])
%export_fig(f ,fullfile(path_dir{1},'example_STCs_Fig3D',[num2str(1),'_300.png']),'-r300')


end

