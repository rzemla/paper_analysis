function [outputArg1,outputArg2] = raster_plot(raster_map,map_color)

%% Color maps

%red and blue colormaps
cmap_blue=cbrewer('seq', 'Blues', 32);
cmap_red=cbrewer('seq', 'Reds', 32);

%set bottom value to red to bottom value of blue
cmap_red(1,:) = [1 1 1];
cmap_blue(1,:) = [1 1 1];

%paper colors
paper_cmap = return_paper_colormap;

%% Plot the rasters

imagesc(raster_map)
hold on
ax1 = gca;
%select colormap
switch map_color 
    case 'blue'
        colormap(ax1,cmap_blue)
    case 'red'
        colormap(ax1,cmap_red)
    case 'jet'
        colormap(ax1,'jet')
end

%plotting param for rate vs activity map
switch map_color
    case 'jet'
        caxis(ax1,[0 3])
        ax1.XTick = [1,30,90,150];
        ax1.XTickLabel = {'0','1','3','5'};
        ax1.YAxis.TickLength = [0 0];
    otherwise
        caxis(ax1,[0 1])
        ax1.XTick = [1,100];
        ax1.XTickLabel = {'0','1'};
        ax1.YAxis.TickLength = [0 0];
        %xticks(ax1,[]);
        ax1.YAxis.TickLength = [0 0];
        
        %B zone
        plot([32 32],[1,size(raster_map,1)], 'Color', paper_cmap(2,:), 'LineStyle','--','LineWidth', 1.5);
        %A zone
        plot([74 74],[1,size(raster_map,1)], 'Color', paper_cmap(1,:), 'LineStyle','--','LineWidth', 1.5);
end


end

