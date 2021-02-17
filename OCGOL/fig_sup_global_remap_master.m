function [outputArg1,outputArg2] = fig_sup_global_remap_master(global_pf_dist,global_dist_scatter,reward_zones_all_animal)

%% Unload data

binCenter_data = global_dist_scatter.binCenter_data;

combined_final_global_dist = global_pf_dist.combined_final_global_dist;
median_dist = global_pf_dist.median_dist;
 
A_zone_end = reward_zones_all_animal.A_zone_end;
B_zone_end = reward_zones_all_animal.B_zone_end;

%% Plot master figure

%paper colors
paper_cmap = return_paper_colormap;

fig = figure;
fig.Units = 'centimeters';
fig.Position(1) = 7;
fig.Position(2) = 0;
fig.Position(3) = 24;
fig.Position(4) = 12;

%master layout
gridSize = [1,2];
t1 = tiledlayout(fig,gridSize(1),gridSize(2),'TileSpacing','normal','Padding','compact','Units','centimeters');

nexttile(t1,1)
hold on
axis square
%title('Global')
xlabel('Normalized location of place field A')
ylabel('Normalized location of place field B')
%generate axis for place fields
xticks(0:10:100)
yticks(0:10:100)
xlabels_norm = 0:0.1:1;
x_labels = {};

for xt=1:11
    x_labels{xt} = num2str(xlabels_norm(xt));
end

xticklabels(x_labels)
yticklabels(x_labels)

for ee=1:numel(binCenter_data)
    scatter(binCenter_data{ee}.bin_center.final.global(1,:),binCenter_data{ee}.bin_center.final.global(2,:),12,'filled',...
            'MarkerEdgeColor',[139,0,139]./255,'MarkerFaceColor',[139,0,139]./255)
end

plot([0 100],[0 100],'k','LineWidth',2)

plot([B_zone_end B_zone_end], [0 100],'--', 'LineWidth',2,'Color',paper_cmap(2,:))
plot([A_zone_end A_zone_end], [0 100],'--', 'LineWidth',2,'Color',paper_cmap(1,:))


nexttile(t1,2)
hold on
axis square
title('Global remapping distance A vs. B')
xlabel('A vs. B center of field distance [cm]')
ylabel('Normalized density')
histogram(cell2mat(combined_final_global_dist).*1.96,0:2.5:100,'Normalization','probability','FaceAlpha',1,'FaceColor',[139,0,139]./255)
%plot median line
plot([median_dist median_dist],[0 0.09],'k','LineWidth',1)

yticks([0 0.02 0.04 0.06 0.08 0.1])
set(gca,'XTick',0:10:100)
xticks(0:10:100)
ylim([0 0.09])
set(gca,'TickDir','Out')
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5);

%set axis font/label and font size
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12, ...
    'FontWeight','normal', 'LineWidth', 1.5,'layer','top')


end

