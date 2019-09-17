function [outputArg1,outputArg2] = plot_STC_OCGOL_singleSes_task_remapping(animal_data, tunedLogical,task_remapping_ROIs,path_dir,options)

%% Import variables

%% Define tuned cell combinations across trials

%make conditional here for si or ts tuned neurons
switch options.tuning_criterion
    
    case 'selective_filtered'
        %neurons idxs associated with selective filtering for
        %task-selectivity
        select_filt_ROIs.A = task_selective_ROIs.A.idx;
        select_filt_ROIs.B = task_selective_ROIs.B.idx;
    case 'remapping_filtered'
        %1 - common
        remap_idx{1} = task_remapping_ROIs.common;
        %2 - rate
        remap_idx{2} = task_remapping_ROIs.rate;
        %3 - global near
        remap_idx{3} = task_remapping_ROIs.global_near;
        %4 - global far
        remap_idx{4} = task_remapping_ROIs.global_far;
        %5 - partial
        remap_idx{5} = task_remapping_ROIs.partial;
        %6 - mixed/unclassified
        remap_idx{6} = task_remapping_ROIs.mixed;
end



%% Extract mean STC map in each spatial bin (not normalized and not occupancy divided) (100 bins)
%for each session
%correct only
for ss =1:size(animal_data,2)
    %normalized to self and not across trials for each neuron
    A_STC_sn{ss} = animal_data{ss}.Place_cell{1}.Spatial_tuning_curve;
    B_STC_sn{ss} = animal_data{ss}.Place_cell{2}.Spatial_tuning_curve;

    %dF/F rasters (normalized 0-1) to self
    A_df_sn{ss} = animal_data{ss}.Place_cell{1}.Spatial_tuning_dF;
    B_df_sn{ss} = animal_data{ss}.Place_cell{2}.Spatial_tuning_dF;

    %100 bins, not smoothed, not norm
    A_df_nn{ss} = animal_data{ss}.Place_cell{1}.Spatial_Info.mean_dF_map{8};
    B_df_nn{ss} = animal_data{ss}.Place_cell{2}.Spatial_Info.mean_dF_map{8};
    %100 bins, Gaussian smoothed across all bins (not 2D smoothed), not
    %norm
    A_df_nn_smooth{ss} = animal_data{ss}.Place_cell{1}.Spatial_Info.mean_dF_map_smooth{8};  
    B_df_nn_smooth{ss} = animal_data{ss}.Place_cell{2}.Spatial_Info.mean_dF_map_smooth{8};
    %Gs smoothed, but not normalized (nn) to itself
    A_STC_nn{ss} = animal_data{ss}.Place_cell{1}.Spatial_Info.rate_map_smooth{8};
    B_STC_nn{ss} = animal_data{ss}.Place_cell{2}.Spatial_Info.rate_map_smooth{8};
    
end

%% Normalize eahc STC ROI across both trials in non-norm STCs

for ss =1:size(animal_data,2)
        %get max value for each ROIs between trials
        max_STC_across_trials{ss} = max([A_STC_nn{ss};B_STC_nn{ss}]);
        %normalize each matrix to these values (tn = trial normalized)
        A_STC_tn{ss} = A_STC_nn{ss}./max_STC_across_trials{ss};
        B_STC_tn{ss} = B_STC_nn{ss}./max_STC_across_trials{ss};
    %for max value to normalize to by 1
    %in future, do normalization range (0,1)
end

%% A vs. B on early vs late training (A or B tuned)

%ALL NEURONS in each session that meet criteria - tuned to either A or B
%early day
%for each type fo trials
%for each remapping cell type
for cc=1:size(remap_idx,2)
    for tt=1:2
        for ss =1:size(animal_data,2)
            if tt ==1
                STC_norm_self_AB{cc}{ss}{tt} = [A_STC_sn{ss}(:,remap_idx{cc})', B_STC_sn{ss}(:,remap_idx{cc})'];
                STC_norm_trials_AB{cc}{ss}{tt} = [A_STC_tn{ss}(:,remap_idx{cc})', B_STC_tn{ss}(:,remap_idx{cc})'];
                dF_nonnorm_AB{cc}{ss}{tt} = [A_df_nn{ss}(:,remap_idx{cc})', B_df_nn{ss}(:,remap_idx{cc})'];
                dF_nonnorm_sm_AB{cc}{ss}{tt} = [A_df_nn_smooth{ss}(:,remap_idx{cc})', B_df_nn_smooth{ss}(:,remap_idx{cc})'];
            elseif tt ==2
                STC_norm_self_AB{cc}{ss}{tt} = [A_STC_sn{ss}(:,remap_idx{cc})', B_STC_sn{ss}(:,remap_idx{cc})'];
                STC_norm_trials_AB{cc}{ss}{tt} = [A_STC_tn{ss}(:,remap_idx{cc})', B_STC_tn{ss}(:,remap_idx{cc})'];
                dF_nonnorm_AB{cc}{ss}{tt} = [A_df_nn{ss}(:,remap_idx{cc})', B_df_nn{ss}(:,remap_idx{cc})'];
                dF_nonnorm_sm_AB{cc}{ss}{tt} = [A_df_nn_smooth{ss}(:,remap_idx{cc})', B_df_nn_smooth{ss}(:,remap_idx{cc})'];
            end
            %STC_norm_trials_AB{ss} = [A_STC_tn{ss}(:,neither_tuned{ss})', B_STC_tn{ss}(:,neither_tuned{ss})'];
        end
    end
end


for cc=1:size(remap_idx,2)
    %sort each session by A map
    for tt=1:2
        for ss =1:size(animal_data,2)
            %change sort order depending of A or B selective neurons being
            %looked at
            %maxBin - spatial bin where activity is greatest for each ROI
            if tt==1
                [~,maxBin_all_AB{cc}{ss}{tt}] = max(STC_norm_trials_AB{cc}{ss}{tt}(:,1:100)', [], 1);
                %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
                [~,sortOrder_all_AB{cc}{ss}{tt}] = sort(maxBin_all_AB{cc}{ss}{tt},'ascend');
                %sort dF/F activity rasters
                [max_dff{cc}{ss}{tt},maxBin_all_AB_dFF{cc}{ss}{tt}] = max(dF_nonnorm_sm_AB{cc}{ss}{tt}(:,1:100)', [], 1);
                [~,sortOrder_all_AB_dFF{cc}{ss}{tt}] = sort(max_dff{cc}{ss}{tt},'descend');
            elseif tt==2
                [~,maxBin_all_AB{cc}{ss}{tt}] = max(STC_norm_trials_AB{cc}{ss}{tt}(:,101:200)', [], 1);
                %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
                [~,sortOrder_all_AB{cc}{ss}{tt}] = sort(maxBin_all_AB{cc}{ss}{tt},'ascend');
                %sort dF/F activity rasters
                [max_dff{cc}{ss}{tt},maxBin_all_AB_dFF{cc}{ss}{tt}] = max(dF_nonnorm_sm_AB{cc}{ss}{tt}(:,101:200)', [], 1);
                [~,sortOrder_all_AB_dFF{cc}{ss}{tt}] = sort(maxBin_all_AB_dFF{cc}{ss}{tt},'descend');
            end
        end
    end
end

%for rate remapping difference - sort by diff of max dF/F between A and B
%selected
[~,sortOrder_rate_remap_only_dFF_diff] = sort(abs(max_dff{3}{1}{1} - max_dff{3}{1}{2}),'descend');


%% Plot normalized STCs against each cell using diff colored maps

%red and blue colormaps
cmap_blue=cbrewer('seq', 'Blues', 32);
cmap_red=cbrewer('seq', 'Reds', 32);

%set bottom value to red to bottom value of blue
cmap_red(1,:) = [1 1 1];
cmap_blue(1,:) = [1 1 1];

%first column: common - global near - partial
%second column: rate - global far - mixed

%% first column plot
ROI_type_order = [1 2; 3 4; 5 6];
subplot_order = [1 3 5; 2 4 6]';
f = figure('Position', [2475 75 372 898]); % event based STC;
for cc =1:3
    
    %subplot(3,2,subplot_order(cc))
    subaxis(3,2,subplot_order(cc,1),'SpacingHorizontal', 0.02,...
        'SpacingVertical',0.01,...
        'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
    
    imagesc(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1}(sortOrder_all_AB{ROI_type_order(cc,1)}{1}{1},1:100))
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
    plot([30 30],[1,size(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([70 70],[1,size(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    
    plot([130 130],[1,size(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([170 170],[1,size(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
end

for cc =1:3
    
    subaxis(3,2,subplot_order(cc,2),'SpacingHorizontal', 0.015,...
        'SpacingVertical',0.01,...
        'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
    imagesc(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1}(sortOrder_all_AB{ROI_type_order(cc,1)}{1}{1},101:200))
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
    plot([30 30],[1,size(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([70 70],[1,size(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    
    plot([130 130],[1,size(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([170 170],[1,size(STC_norm_trials_AB{ROI_type_order(cc,1)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    hold off
end

%export rasters fpr first colum
mkdir(fullfile(path_dir{1},'example_STCs_Fig3D'))
disp(['Saving raster: 1'])
export_fig(f ,fullfile(path_dir{1},'example_STCs_Fig3D',[num2str(1),'_300.png']),'-r300')

%% second column plot
ROI_type_order = [1 2; 3 4; 5 6];
subplot_order = [1 3 5; 2 4 6]';
f = figure('Position', [2475 75 372 898]); % event based STC;
for cc =1:3
    
    %subplot(3,2,subplot_order(cc))
subaxis(3,2,subplot_order(cc,1),'SpacingHorizontal', 0.02,...
    'SpacingVertical',0.01,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);

    imagesc(STC_norm_trials_AB{ROI_type_order(cc,2)}{1}{1}(sortOrder_all_AB{ROI_type_order(cc,2)}{1}{1},1:100))
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
    plot([30 30],[1,size(STC_norm_trials_AB{ROI_type_order(cc,2)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([70 70],[1,size(STC_norm_trials_AB{ROI_type_order(cc,2)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    
    plot([130 130],[1,size(STC_norm_trials_AB{ROI_type_order(cc,2)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([170 170],[1,size(STC_norm_trials_AB{ROI_type_order(cc,2)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
end

for cc =1:3
    
    subaxis(3,2,subplot_order(cc,2),'SpacingHorizontal', 0.015,...
    'SpacingVertical',0.01,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
    imagesc(STC_norm_trials_AB{ROI_type_order(cc,2)}{1}{1}(sortOrder_all_AB{ROI_type_order(cc,2)}{1}{1},101:200))
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
    plot([30 30],[1,size(STC_norm_trials_AB{ROI_type_order(cc,2)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([70 70],[1,size(STC_norm_trials_AB{ROI_type_order(cc,2)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    
    plot([130 130],[1,size(STC_norm_trials_AB{ROI_type_order(cc,2)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([170 170],[1,size(STC_norm_trials_AB{ROI_type_order(cc,2)}{1}{1},1)], 'Color', 0*[1 1 1], 'LineStyle','--','LineWidth', 1.5);
    hold off
end

%export rasters fpr first colum
mkdir(fullfile(path_dir{1},'example_STCs_Fig3D'))
disp(['Saving raster: 2'])
export_fig(f ,fullfile(path_dir{1},'example_STCs_Fig3D',[num2str(2),'_300.png']),'-r300')

%% plot side by side; day by day
subplot_order = [1 2; 3 4; 5 6; 7 8; 9 10; 11 12];
f = figure('Position', [2015 84 851 898]); % event based STC;
for cc =1:6
    
    subplot(6,2,subplot_order(cc,1))
    imagesc(STC_norm_trials_AB{cc}{1}{1}(sortOrder_all_AB{cc}{1}{1},:))
    %title('5A5B')
    hold on
    cbar= colorbar;
    cbar.Label.String = 'Normalized activity';
    cbar.Ticks = [0 1];
    colormap('jet')
    caxis([0 1])
    %A/B vertical separator line
    plot([100 100],[1,size(STC_norm_trials_AB{cc}{1}{1},1)], 'k','LineWidth', 1.5);
    %plot reward zones as dashed lines
    %B zone
    plot([30 30],[1,size(STC_norm_trials_AB{cc}{1}{1},1)], 'Color', [1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([70 70],[1,size(STC_norm_trials_AB{cc}{1}{1},1)], 'Color', [1 1 1], 'LineStyle','--','LineWidth', 1.5);
    
    plot([130 130],[1,size(STC_norm_trials_AB{cc}{1}{1},1)], 'Color', [1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([170 170],[1,size(STC_norm_trials_AB{cc}{1}{1},1)], 'Color', [1 1 1], 'LineStyle','--','LineWidth', 1.5);
    
    subplot(6,2,subplot_order(cc,2))
    if cc==3 %ignores and sort by STC order
        %imagesc(dF_nonnorm_sm_AB{cc}{1}{1}(sortOrder_all_AB_dFF{cc}{1}{1},:))
        %imagesc(dF_nonnorm_sm_AB{cc}{1}{1}(sortOrder_rate_remap_only_dFF_diff,:))
        imagesc(dF_nonnorm_sm_AB{cc}{1}{1}(sortOrder_all_AB{cc}{1}{1},:))
    else
        imagesc(dF_nonnorm_sm_AB{cc}{1}{1}(sortOrder_all_AB{cc}{1}{1},:))
    end
    %title('Normalized according to trials')
    hold on
    caxis([0 1.5])
    colormap('jet')
    cbar= colorbar;
    cbar.Label.String = 'dF/F';
    cbar.Ticks = [0 0.5 1];
    ax1 = gca;
    ylabel('Neuron #');
    xlabel('Normalized position');
    ax1.XTick = [1 100 200];
    ax1.XTickLabel = {'0','1','1'};
    %A/B vertical separator line
    plot([100 100],[1,size(STC_norm_trials_AB{cc}{1}{1},1)], 'k','LineWidth', 1.5);
    %plot reward zones as dashed lines
    %B zone
    plot([30 30],[1,size(STC_norm_trials_AB{cc}{1}{1},1)], 'Color', [1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([70 70],[1,size(STC_norm_trials_AB{cc}{1}{1},1)], 'Color', [1 1 1], 'LineStyle','--','LineWidth', 1.5);
    
    plot([130 130],[1,size(STC_norm_trials_AB{cc}{1}{1},1)], 'Color', [1 1 1], 'LineStyle','--','LineWidth', 1.5);
    %A zone
    plot([170 170],[1,size(STC_norm_trials_AB{cc}{1}{1},1)], 'Color', [1 1 1], 'LineStyle','--','LineWidth', 1.5);
    
end


%% figure % event based STC;
%subplot(2,1,1)
% imagesc(STC_norm_self_AB.global{1}{1}(sortOrder_all_AB{1}{1},:))
% %title('5A5B')
% hold on
% colormap('jet')
% caxis([0 1])
% %A/B vertical separator line
% plot([100 100],[1,size(STC_norm_self_AB.global{1}{1},1)], 'k','LineWidth', 1.5);
% %plot reward zones as dashed lines
% %B zone
% plot([30 30],[1,size(STC_norm_self_AB.global{1}{1},1)], 'Color', [1 1 1], 'LineStyle','--','LineWidth', 1.5);
% %A zone
% plot([70 70],[1,size(STC_norm_self_AB.global{1}{1},1)], 'Color', [1 1 1], 'LineStyle','--','LineWidth', 1.5);
% 
% plot([130 130],[1,size(STC_norm_self_AB.global{1}{1},1)], 'Color', [1 1 1], 'LineStyle','--','LineWidth', 1.5);
% %A zone
% plot([170 170],[1,size(STC_norm_self_AB.global{1}{1},1)], 'Color', [1 1 1], 'LineStyle','--','LineWidth', 1.5);


%% Normalized between trials

% %STC normalized across trials
% f= figure('Position', [2090 415 1240 420]);
% subplot(1,2,1)
% imagesc(STC_norm_trials_AB.global{1}{1}(sortOrder_all_AB{1}{1},:))
% %title('Normalized according to trials')
% hold on
% caxis([0 1])
% colormap('jet')
% cbar= colorbar;
% cbar.Label.String = 'Normalized activity';
% cbar.Ticks = [0 0.5 1];
% ax1 = gca;
% ylabel('Neuron #');
% xlabel('Normalized position');
% ax1.XTick = [1 100 200];
% ax1.XTickLabel = {'0','1','1'};
% %A/B vertical separator line
% plot([100 100],[1,size(STC_norm_trials_AB.global{1}{1},1)], 'k','LineWidth', 1.5);
% hold off
% 
% subplot(1,2,2)
% imagesc(STC_norm_trials_AB.global{1}{2}(sortOrder_all_AB{1}{1},:))
% hold on
% caxis([0 1])
% colormap('jet')
% cbar= colorbar;
% cbar.Label.String = 'Normalized activity';
% cbar.Ticks = [0 0.5 1];
% ax1 = gca;
% ylabel('Neuron #');
% xlabel('Normalized position');
% ax1.XTick = [1 100 200];
% ax1.XTickLabel = {'0','1','1'};
% %A/B vertical separator line
% plot([100 100],[1,size(STC_norm_trials_AB.global{1}{1},1)], 'k','LineWidth', 1.5);
% hold off

%% Non normalized dF/F

% %STC normalized across trials
% f= figure('Position', [2090 415 1240 420]);
% %subplot(1,2,1)
% imagesc(dF_nonnorm_sm_AB.global{1}{1}(sortOrder_all_AB{1}{1},:))
% %title('Normalized according to trials')
% hold on
% caxis([0 2])
% colormap('jet')
% cbar= colorbar;
% cbar.Label.String = 'dF/F';
% cbar.Ticks = [0 0.5 1];
% ax1 = gca;
% ylabel('Neuron #');
% xlabel('Normalized position');
% ax1.XTick = [1 100 200];
% ax1.XTickLabel = {'0','1','1'};
% %A/B vertical separator line
% plot([100 100],[1,size(STC_norm_trials_AB.global{1}{1},1)], 'k','LineWidth', 1.5);
% hold off

% subplot(1,2,2)
% imagesc(STC_norm_trials_AB.global{1}{2}(sortOrder_all_AB{1}{2},:))
% hold on
% caxis([0 1])
% colormap('jet')
% cbar= colorbar;
% cbar.Label.String = 'Normalized activity';
% cbar.Ticks = [0 0.5 1];
% ax1 = gca;
% ylabel('Neuron #');
% xlabel('Normalized position');
% ax1.XTick = [1 100 200];
% ax1.XTickLabel = {'0','1','1'};
% %A/B vertical separator line
% plot([100 100],[1,size(STC_norm_trials_AB.global{1}{1},1)], 'k','LineWidth', 1.5);
% hold off

% subplot(2,1,2)
% imagesc(dF_maps_all_AB_early_late{2}(sortOrder_all_AB{2},:))
% hold on
% title('Random AB')
% colormap('jet')
% caxis([0 1])
% %A/B vertical separator line
% plot([100 100],[1,size(dF_maps_all_AB_early_late{2},1)], 'k','LineWidth', 1.5);


%hold off
%% PV correlation plot below
%PVcorr = corr(A_STC_nn{session_nb}(:,tuning_selection{session_nb}),B_STC_nn{session_nb}(:,tuning_selection{session_nb}))
% PVcorr = corr(A_STC_nn{1}',B_STC_nn{1}');
% 
% 
% figure('Position',[1350, 90, 500 860]);
% subplot(2,1,1)
% imagesc(PVcorr)
% hold on
% title('Population vector correlation');
% xlabel('Spatial bin')
% ylabel('Spatial bin')
% colormap('jet')
% caxis([0 1])
% xticks([20 40 60 80 100]);
% yticks([20 40 60 80 100]);
% axis('square')
% cbar2 = colorbar;
% cbar2.Label.String = 'Correlation coefficient';
% cbar2.Ticks = [0 0.5 1];
% 
% subplot(2,1,2)
% hold on
% title('Population vector correlation');
% plot(diag(PVcorr),'k','LineWidth',1.5)
% xlabel('Spatial bin')
% ylabel('Correlation coef.');
% plot([30 30],[0 1],'r--','LineWidth',1.5);
% text([31 31],[0.9 0.9],'Reward zone B','Color','r')
% plot([70 70],[0 1],'b--','LineWidth',1.5);
% text([71 71],[0.3 0.3],'Reward zone A','Color','b')
% plot([10 10],[0 1],'g--','LineWidth',1.5);
% text([11 11],[0.9 0.9],'Odor zone\newline end','Color','g')

%% Make matching ROI list with tuning criteria for both sessions

%tuned to A or B on either sessions
% AorB_idx{1} = find(AorB_tuned{1} ==1);
% AorB_idx{2} = find(AorB_tuned{2} ==1);
% 
% %intersect with
% [tuned_match_idx{1},match_idx{1},~] = intersect(matching_list(:,1),AorB_idx{1},'stable');
% [tuned_match_idx{2},match_idx{2},~] = intersect(matching_list(:,2),AorB_idx{2},'stable');
% 
% %create not logical for nan exclusion from copied matrix assignement below
% include_log{1} = false(1,size(matching_list,1));
% include_log{1}(match_idx{1}) = 1; 
% %session 2 
% include_log{2} = false(1,size(matching_list,1));
% include_log{2}(match_idx{2}) = 1;
% 
% %make copy
% tuned_matching_ROI_list = matching_list;
% %nan first session that are not tuned and last session that are not
% %tuned
% tuned_matching_ROI_list(~include_log{1},1) = nan;
% tuned_matching_ROI_list(~include_log{2},2) = nan;
% 
% %which neurons to remove based on tuning criterion
% keep_ROI = sum(isnan(tuned_matching_ROI_list),2) == 0;
% 
% %retain only tuned and matched ROIs
% tuned_matching_ROI_list(~keep_ROI,:) = [];

%% Generate maps based on tuned matching ROI list
%row - session
%column - trial type
% 
% session_matched_tuned_dF_maps{1,1} = A_STC_sn{1}(:,tuned_matching_ROI_list(:,1))';
% session_matched_tuned_dF_maps{1,2} = B_STC_sn{1}(:,tuned_matching_ROI_list(:,1))';
% session_matched_tuned_dF_maps{2,1} = A_STC_sn{2}(:,tuned_matching_ROI_list(:,2))';
% session_matched_tuned_dF_maps{2,2} = B_STC_sn{2}(:,tuned_matching_ROI_list(:,2))';

% %combined 2x2
% combined_maps_2x2 = cell2mat(session_matched_tuned_dF_maps);
% 
% %single matrix (ses 1 A, B, ses 2 A,B) - in 1 row
% combined_maps_row = cell2mat(reshape(session_matched_tuned_dF_maps,1,4));
% 
% %sort by A trials on session 1
% [~,matched_maxBin_1] = max(session_matched_tuned_dF_maps{1,1}(:,1:100)', [], 1);
% %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
% [~,matched_sortOrder_1] = sort(matched_maxBin_1,'ascend');
% 
% %sort by A trials on session 2
% [~,matched_maxBin_2] = max(session_matched_tuned_dF_maps{2,1}(:,1:100)', [], 1);
% %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
% [~,matched_sortOrder_2] = sort(matched_maxBin_2,'ascend');

%session 1 above and session 2 below
% figure;
% subplot(2,2,1)
% imagesc(session_matched_tuned_dF_maps{1,1}(matched_sortOrder_1,:))
% colormap('jet')
% caxis([0 1]);
% subplot(2,2,2)
% imagesc(session_matched_tuned_dF_maps{1,2}(matched_sortOrder_1,:))
% colormap('jet')
% caxis([0 1]);
% subplot(2,2,3)
% imagesc(session_matched_tuned_dF_maps{2,1}(matched_sortOrder_1,:))
% colormap('jet')
% caxis([0 1]);
% subplot(2,2,4)
% imagesc(session_matched_tuned_dF_maps{2,2}(matched_sortOrder_1,:))
% colormap('jet')
% caxis([0 1]);

% %horizontal and vertical separation lines
% plot([100 100],[1,size(session_matched_tuned_dF_maps{1,1},1)*2], 'k','LineWidth', 1.5);
% plot([100 100],[1,size(dF_maps_all_AB_early_late{1},1)], 'k','LineWidth', 1.5);
% 

%% Extract STCs with tuned ROIs - in nontuned neurons will scale the weakest signal to 1 regardless of tuning

%definition of spatial tuning curve
%Gaussian smoothed onset rate map / spatial bin occupancy time (sec)
%Normalization from (0-1) for each ROI (ROI-by-ROI)

%these contain NaNs
% STC_A = animal_data{1}.Place_cell{1}.Spatial_tuning_curve;
% STC_B = animal_data{1}.Place_cell{2}.Spatial_tuning_curve;
% 
% A_STC_both = STC_A(:,AandB_tuned);
% B_STC_both = STC_B(:,AandB_tuned);
% 
% A_STC_onlyA = STC_A(:,onlyA_tuned);
% B_STC_onlyA = STC_B(:,onlyA_tuned);



end

