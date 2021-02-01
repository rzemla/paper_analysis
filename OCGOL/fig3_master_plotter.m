function [outputArg1,outputArg2] = fig3_master_plotter(inputArg1,inputArg2)



%% Unload the data






%% Figure 3c rate maps

%red and blue colormaps
cmap_blue=cbrewer('seq', 'Blues', 32);
cmap_red=cbrewer('seq', 'Reds', 32);

%set bottom value to red to bottom value of blue
cmap_red(1,:) = [1 1 1];
cmap_blue(1,:) = [1 1 1];

%% first column plot
ROI_type_order = [1 2; 3 4; 5 6];
subplot_order = [1 3 5; 2 4 6]';
f = figure('Position', [2475 75 372 898]); % event based STC;
for cc =1:3
    
    %subplot(3,2,subplot_order(cc))
    subaxis(3,2,subplot_order(cc,1),'SpacingHorizontal', 0.02,...
        'SpacingVertical',0.01,...
        'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
    
    imagesc(first_col_STCs{cc,1}(:,1:100))
    hold on
    ax1 = gca;
    %set(gca,'FontSize',14)
    %title('5A5B')
    %c1= colorbar;
    %cbar.Label.String = 'Normalized activity';
    %cbar.Ticks = [0 1];
    colormap(ax1,cmap_blue)
    caxis(ax1,[0 1])
    if cc ~= 3
        xticks(ax1,[]);
        %ax1.YTick = [];
        ax1.YAxis.TickLength = [0 0];
    else
        %ax1.YTick = [1,100];
        ax1.XTick = [1,100];
        ax1.XTickLabel = {'0','1'};
        ax1.YAxis.TickLength = [0 0];
    end
    %A/B vertical separator line
    %plot([100 100],[1,size(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1},1)], 'k','LineWidth', 1.5);
    %plot reward zones as dashed lines
    
    %B zone
    plot([32 32],[1,size(first_col_STCs{cc}(:,1:100),1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([74 74],[1,size(first_col_STCs{cc}(:,1:100),1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    
    %B zone
    %plot([132 132],[1,size(first_col_STCs{cc}(:,1:100),1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    %plot([174 174],[1,size(first_col_STCs{cc}(:,1:100),1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    
end

for cc =1:3
    
    subaxis(3,2,subplot_order(cc,2),'SpacingHorizontal', 0.015,...
        'SpacingVertical',0.01,...
        'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
    imagesc(first_col_STCs{cc,1}(:,101:200))
    %title('5A5B')
    hold on
    ax2 = gca;
    %cbar= colorbar;
    %cbar.Label.String = 'Normalized activity';
    %cbar.Ticks = [0 1];
    colormap(ax2,cmap_red)
    caxis(ax2,[0 1])
    if cc ~= 3
        xticks(ax2,[]);
        ax2.YTick = [];
        ax2.YAxis.TickLength = [0 0];
    else
        ax2.YTick = [];%[1,100];
        ax2.XTick = [1 100];%[1,100];
        ax2.XTickLabel = {'0','1'};
        ax2.YAxis.TickLength = [0 0];
    end
    %yticks(ax2,[]);
    %A/B vertical separator line
    %plot([100 100],[1,size(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1},1)], 'k','LineWidth', 1.5);
    %plot reward zones as dashed lines
    
    %B zone
    plot([32 32],[1,size(first_col_STCs{cc}(:,1:100),1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([74 74],[1,size(first_col_STCs{cc}(:,1:100),1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    
    %plot([130 130],[1,size(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    %plot([170 170],[1,size(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %hold off
end

%export rasters fpr first colum
%mkdir(fullfile(path_dir{1},'example_STCs_Fig3D'))
%disp(['Saving raster: 1'])
%export_fig(f ,fullfile(path_dir{1},'example_STCs_Fig3D',[num2str(1),'_300.png']),'-r300')

%% second column plot
ROI_type_order = [1 2; 3 4; 5 6];
subplot_order = [1 3 5; 2 4 6]';
f = figure('Position', [2475 75 372 898]); % event based STC;
for cc =1:2
    
    %subplot(3,2,subplot_order(cc))
subaxis(3,2,subplot_order(cc,1),'SpacingHorizontal', 0.02,...
    'SpacingVertical',0.01,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);

    imagesc(second_col_STCs{cc,1}(:,1:100))
    hold on
    ax1 = gca;
    %title('5A5B')
    %c1= colorbar;
    %cbar.Label.String = 'Normalized activity';
    %cbar.Ticks = [0 1];
    colormap(ax1,cmap_blue)
    caxis(ax1,[0 1])
    if cc ~= 3
        xticks(ax1,[]);
        %ax1.YTick = [];
        ax1.YAxis.TickLength = [0 0];
    else
        %ax1.YTick = [1,100];
        ax1.XTick = [1,100];
        ax1.XTickLabel = {'0','1'};
 ax1.YAxis.TickLength = [0 0];
    end
    %A/B vertical separator line
    %plot([100 100],[1,size(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1},1)], 'k','LineWidth', 1.5);
    %plot reward zones as dashed lines

    %B zone
    plot([32 32],[1,size(second_col_STCs{cc}(:,1:100),1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([74 74],[1,size(second_col_STCs{cc}(:,1:100),1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
        
    
    %plot([130 130],[1,size(STC_norm_trials_AB{ROI_type_order(cc,2)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    %plot([170 170],[1,size(STC_norm_trials_AB{ROI_type_order(cc,2)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
end

for cc =1:2
    
    subaxis(3,2,subplot_order(cc,2),'SpacingHorizontal', 0.015,...
    'SpacingVertical',0.01,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
    imagesc(second_col_STCs{cc,1}(:,101:200))
    %title('5A5B')
    hold on
    ax2 = gca;
    %cbar= colorbar;
    %cbar.Label.String = 'Normalized activity';
    %cbar.Ticks = [0 1];
    colormap(ax2,cmap_red)
    caxis(ax2,[0 1])
    if cc ~= 3
        xticks(ax2,[]);
        ax2.YTick = [];
        ax2.YAxis.TickLength = [0 0];
    else
        ax2.YTick = [];%[1,100];
        ax2.XTick = [1 100];%[1,100];
    ax2.XTickLabel = {'0','1'};
ax2.YAxis.TickLength = [0 0];
    end
    %yticks(ax2,[]);
    %A/B vertical separator line
    %plot([100 100],[1,size(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1},1)], 'k','LineWidth', 1.5);
    %plot reward zones as dashed lines

    %B zone
    plot([32 32],[1,size(second_col_STCs{cc}(:,1:100),1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([74 74],[1,size(second_col_STCs{cc}(:,1:100),1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
     
    
    %plot([130 130],[1,size(STC_norm_trials_AB{ROI_type_order(cc,2)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    %plot([170 170],[1,size(STC_norm_trials_AB{ROI_type_order(cc,2)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    hold off
end
end

