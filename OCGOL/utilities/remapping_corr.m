function [cutoffs_95] = remapping_corr(path_dir)

%% Load in remapping idx data

for aa=1:size(path_dir,2)
    load_data_path{aa} = fullfile(path_dir{aa},'cumul_analysis','remap_corr_idx.mat');
    remap_corr_idx{aa} = load(string(load_data_path{aa}));
end

%% Load each ROI STC data
for aa=1:size(path_dir,2)
    load_STC_path{aa} = fullfile(path_dir{aa},'cumul_analysis','STC.mat');
    STC{aa} = load(string(load_STC_path{aa}));
end

%% Load in place field separation distances for common neurons
for aa=1:size(path_dir,2)
    load_distance_path{aa} = fullfile(path_dir{aa},'cumul_analysis','place_field_common_distances.mat');
    distances{aa} = load(string(load_distance_path{aa}));
end

%% Load in bin centers (for partial remapper sorting)
for aa=1:size(path_dir,2)
    load_bin_center_path{aa} = fullfile(path_dir{aa},'cumul_analysis','place_field_centers_remap.mat');
    bin_centers{aa} = load(string(load_bin_center_path{aa}));
end

%% Load mean dF/F
%dF/F maps for rate remapping neurons
for aa=1:size(path_dir,2)
    load_dFF_remap_path{aa} = fullfile(path_dir{aa},'cumul_analysis','rate_mean_dFF.mat');
    rate_mean_dFF{aa} = load(string(load_dFF_remap_path{aa}));
    %load out of struct
    rate_mean_dFF_no_struct{aa} = rate_mean_dFF{aa}.rate_mean_dFF;
end

rate_dFF_combined = cell2mat(rate_mean_dFF_no_struct');

%load(fullfile(path_dir{1},'cumul_analysis','rate_mean_dFF.mat'),'rate_mean_dFF');


%save(fullfile(path_dir{1},'cumul_analysis','place_field_centers_remap.mat'),'bin_center');
%save(fullfile(path_dir{1},'cumul_analysis','place_field_common_distances.mat'),'common_bin_conv_diff','common_pf_distance_metric');

%% Combine distances and compute the 95th percentile cutoff

for aa=1:size(path_dir,2)
    bin_dist{aa} = distances{aa}.common_bin_conv_diff;
end

%get 95% bin cutoff for fields that are considered common (for defining
%partial neurons)
bin_cutoff_95 = quantile(cell2mat(bin_dist),0.95);
cm_cutoff_95 = bin_cutoff_95*1.96;
ang_cutoff_95 = rad2deg((bin_cutoff_95./100)*2*pi);

%export cutoffs
cutoffs_95.bin = bin_cutoff_95;
cutoffs_95.cm = cm_cutoff_95;
cutoffs_95.ang = ang_cutoff_95;

%plot histogram (supplement)
figure
hold on
histogram(cell2mat(bin_dist),30)

%% Combined into single matrix
for ss=1:size(remap_corr_idx,2)
    %non_global_idx{ss} = remap_corr_idx{ss}.remapping_corr_idx.non_global_idx;
    %global
    global_idx{ss} = remap_corr_idx{ss}.remapping_corr_idx.final.global;
    %rate remapping
    rate_idx{ss} = remap_corr_idx{ss}.remapping_corr_idx.final.rate_remap_all;
    %common
    common_idx{ss} = remap_corr_idx{ss}.remapping_corr_idx.final.common;
    %only signifcance for rate and nothing else (speed or interaction)
    rate_only_idx{ss} = remap_corr_idx{ss}.remapping_corr_idx.final.rate_remap_grp_only;
    
    %partial remapping
    partial_idx{ss} = remap_corr_idx{ss}.remapping_corr_idx.final.partial;
    
    %unclassified
    unclass_idx{ss} = remap_corr_idx{ss}.remapping_corr_idx.final.unclass;
    
    %use this as only interested (group signifance regardless of speed or
    %interaction)
    
    %common
    %common_idx{ss} = remap_corr_idx{ss}.remapping_corr_idx.final.common;
    %common_idx{ss} = setdiff(non_global_idx{ss},rate_only_idx{ss});
end

%plot the STCs associated with each class
%non_global_comb = cell2mat(non_global_idx');
global_comb = cell2mat(global_idx');
rate_only_comb = cell2mat(rate_idx);

%common neurons 
common_comb = cell2mat(common_idx');

%number of neurons that are only task category modulated
cell2mat(rate_only_idx);

%partial
partial_comb = cell2mat(partial_idx);

%unclassfied com
unclass_comb = cell2mat(unclass_idx);



%% Extract each class STC
for aa=1:size(path_dir,2)
    %A vs B side by side
    STC_tn_global{aa} = [STC{aa}.STC_export.A_STC_tn{1}(:,global_idx{aa})', STC{aa}.STC_export.B_STC_tn{1}(:,global_idx{aa})'];
    STC_tn_common{aa} = [STC{aa}.STC_export.A_STC_tn{1}(:,common_idx{aa})', STC{aa}.STC_export.B_STC_tn{1}(:,common_idx{aa})'];
    STC_tn_rate{aa} = [STC{aa}.STC_export.A_STC_tn{1}(:,rate_idx{aa})', STC{aa}.STC_export.B_STC_tn{1}(:,rate_idx{aa})'];
    STC_tn_partial{aa} = [STC{aa}.STC_export.A_STC_tn{1}(:,partial_idx{aa})', STC{aa}.STC_export.B_STC_tn{1}(:,partial_idx{aa})'];
    STC_tn_unclass{aa} = [STC{aa}.STC_export.A_STC_tn{1}(:,unclass_idx{aa})', STC{aa}.STC_export.B_STC_tn{1}(:,unclass_idx{aa})'];
end

%% Sort partial idx's by position of common bin in A
for aa=1:size(path_dir,2)
    [~,sortOrder_partial{aa}] = sort(bin_centers{aa}.bin_center.partial_com(1,:),'ascend');
    %combined all A centroid positions into cells and sort below
    bin_centers_partial_com{aa} = bin_centers{aa}.bin_center.partial_com(1,:)
end

%sort all indices
[~,A_partial_com_sort_all] = sort(cell2mat(bin_centers_partial_com),'ascend');

%% Sort STCs for each session
%global sort
for aa=1:size(path_dir,2)
    %sort global
    %input: ROI x bin (concatenated from both A and B trials)
    sortOrder_global{aa} = sortSTC(STC_tn_global{aa},1:100);
    
    sortOrder_common{aa} = sortSTC(STC_tn_common{aa},1:100);
    
    sortOrder_rate{aa} = sortSTC(STC_tn_rate{aa},1:100);
    
    sortOrder_unclass{aa} = sortSTC(STC_tn_unclass{aa},1:100);
    
    %use the common bins for A to sort each (above)
    %sortOrder_partial{aa} = sortSTC(STC_tn_partial{aa},1:100);
    
end

%sort each global STC
for aa=1:size(path_dir,2)
    STC_tn_global_sorted{aa} = STC_tn_global{aa}(sortOrder_global{aa},:);
    
    STC_tn_common_sorted{aa} = STC_tn_common{aa}(sortOrder_common{aa},:);
    
    STC_tn_rate_sorted{aa} = STC_tn_rate{aa}(sortOrder_rate{aa},:);
    
    STC_tn_partial_sorted{aa} = STC_tn_partial{aa}(sortOrder_partial{aa},:);
    
    STC_tn_unclass_sorted{aa} = STC_tn_unclass{aa}(sortOrder_unclass{aa},:);
    
end

%% Sort all neurons combined
%global
STC_tn_global_all = cell2mat(STC_tn_global_sorted');
STC_tn_global_all_sorted = STC_tn_global_all(sortSTC(STC_tn_global_all,1:100),:);

%common
STC_tn_common_all = cell2mat(STC_tn_common_sorted');
STC_tn_common_all_sorted = STC_tn_common_all(sortSTC(STC_tn_common_all,1:100),:);

%rate
STC_tn_rate_all = cell2mat(STC_tn_rate_sorted');
STC_tn_rate_all_sorted = STC_tn_rate_all(sortSTC(STC_tn_rate_all,1:100),:);


%unclass
STC_tn_unclass_all = cell2mat(STC_tn_unclass_sorted');
STC_tn_unclass_all_sorted = STC_tn_unclass_all(sortSTC(STC_tn_unclass_all,1:100),:);

%partial (sorted by each animal)
STC_tn_partial_all = cell2mat(STC_tn_partial_sorted');

%
STC_tn_partial_nonsorted = cell2mat(STC_tn_partial');

%use common idx of A tuned neurons to sort
STC_tn_partial_all_sorted = STC_tn_partial_nonsorted(A_partial_com_sort_all,:);

%% dF/F STC plot for each neuron


%% Gaussian kernel smoothing filter
%using convolution with custom window
%size of window
%rename option variable to use for input for define_Gaussian_kernel fxn
options.sigma_filter = 3;
gaussFilter = define_Gaussian_kernel(options);

%% Smooth extended_rate_map

%return sort order for rate remapping neurons and plot associated dF/F
rate_dFF_mean_sorted_all = rate_dFF_combined(sortSTC(rate_dFF_combined,1:100),:);

%split into 2 raster and smooth individually
rate_dFF_mean_sorted.A = rate_dFF_mean_sorted_all(:,1:100);
rate_dFF_mean_sorted.B = rate_dFF_mean_sorted_all(:,101:200);

%smooth extended rate map for each ROI
for rr = 1:64
    %smooth A
    gauss_smoothed_dFF.A(rr,:) = conv(rate_dFF_mean_sorted.A(rr,:),gaussFilter, 'same');
    %smooth B
    gauss_smoothed_dFF.B(rr,:) = conv(rate_dFF_mean_sorted.B(rr,:),gaussFilter, 'same');
end

if 0
figure
for ii=1:64
    hold on
    plot(gauss_smoothed_dFF.A(ii,1:100),'b')
    plot(gauss_smoothed_dFF.B(ii,1:100),'r')
    pause
    clf
end
end

%sort by peak difference of dF/F
 [~,sort_max_dff_diff] = sort(max(gauss_smoothed_dFF.A,[],2) - max(gauss_smoothed_dFF.B,[],2),'ascend')

 if 0
figure
hold on
idx_step = 0;

for ii=sort_max_dff_diff(51:60)'
    plot(gauss_smoothed_dFF.A(ii,:)+idx_step,'b','LineWidth',1)
    plot(gauss_smoothed_dFF.B(ii,:)+idx_step,'r','LineWidth',1)
    idx_step = idx_step +2;
end


figure
subplot(1,2,1)
imagesc(gauss_smoothed_dFF.A)
hold on
caxis([0 2])
colormap('jet')
subplot(1,2,2)
imagesc(gauss_smoothed_dFF.B)
hold on
caxis([0 2])
colormap('jet')
 end
 
%% Compare partial remapping neurons curve by curve

if 0
figure
for rr=1:427
    hold on
    plot(STC_tn_partial_all_sorted(rr,1:100),'b')
    plot(STC_tn_partial_all_sorted(rr,101:200),'r')
    %pause
    clf
end
end

%% Plot individually sorted global, common, rate STCs

% figure
% imagesc(cell2mat(STC_tn_common_sorted'))
% hold on
% title('Common')
% caxis([0 1])
% colormap('jet')
% 
% figure
% imagesc(cell2mat(STC_tn_global_sorted'))
% hold on
% title('Global remappers')
% caxis([0 1])
% colormap('jet')
% 
% figure
% imagesc(cell2mat(STC_tn_rate_sorted'))
% hold on
% title('Rate')
% caxis([0 1])
% colormap('jet')
% 
% figure
% imagesc(cell2mat(STC_tn_partial_sorted'))
% hold on
% title('Partial')
% caxis([0 1])
% colormap('jet')



%% Plot combined all sorted global, common, rate STCs

figure
imagesc(STC_tn_common_all_sorted)
hold on
title('Common')
caxis([0 1])
colormap('jet')

figure
imagesc(STC_tn_global_all_sorted)
hold on
title('Global remappers')
caxis([0 1])
colormap('jet')

figure
imagesc(STC_tn_rate_all_sorted)
hold on
title('Rate')
caxis([0 1])
colormap('jet')

figure
imagesc(STC_tn_partial_all_sorted)
hold on
title('Partial')
caxis([0 1])
colormap('jet')

figure
imagesc(STC_tn_unclass_all_sorted)
hold on
title('Unclassified')
caxis([0 1])
colormap('jet')

%% Plot normalized STCs against each cell using diff colored maps

STC_tn_common_all_sorted
STC_tn_global_all_sorted
STC_tn_rate_all_sorted
STC_tn_partial_all_sorted
STC_tn_unclass_all_sorted

%red and blue colormaps
cmap_blue=cbrewer('seq', 'Blues', 32);
cmap_red=cbrewer('seq', 'Reds', 32);

%set bottom value to red to bottom value of blue
cmap_red(1,:) = [1 1 1];
cmap_blue(1,:) = [1 1 1];

%first column: common - rate - global 
first_col_STCs = {STC_tn_common_all_sorted; STC_tn_rate_all_sorted; STC_tn_global_all_sorted};
second_col_STCs = {STC_tn_partial_all_sorted; STC_tn_unclass_all_sorted};
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

%export rasters fpr first colum
%mkdir(fullfile(path_dir{1},'example_STCs_Fig3D'))
%disp(['Saving raster: 2'])
%export_fig(f ,fullfile(path_dir{1},'example_STCs_Fig3D',[num2str(2),'_300.png']),'-r300')


end

