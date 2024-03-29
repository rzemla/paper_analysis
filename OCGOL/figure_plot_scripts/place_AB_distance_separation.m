function [combined_distances,combined_dist_metric,global_pf_dist] = place_AB_distance_separation(path_dir)

%% Load in place field distance data

for aa=1:size(path_dir,2)
    load_data_path{aa} = fullfile(path_dir{aa},'cumul_analysis','place_field_AB_distances.mat');
    distances{aa} = load(string(load_data_path{aa}));
end

%load(fullfile(path_dir{1},'cumul_analysis','place_field_AB_distances.mat'),'ts_bin_conv_diff','ts_bin_conv_diff');

%% Combined into single matrix

for ss=1:size(distances,2)
    combined_distances{ss} = distances{ss}.ts_bin_conv_diff;
    combined_dist_metric{ss} = distances{ss}.pf_distance_metric_ts;
    combined_final_global_dist{ss} = distances{ss}.final_global_remap_bin_diff;   
end

%% Export the data loaded in from all animals

%each bin ~ 1.96 cm (for 100 bins) 
cm_conv_factor = 1.96;

%generate the cm distance separation for each animal
for aa=1:size(path_dir,2)
    global_dist_cm{aa} = combined_final_global_dist{aa}.*cm_conv_factor;
end

%plot the histogram for each animal and choose based on visual inspection
%animal to use data from; display number of global remappers for each
%animal
figure('Position', [2058 101 586 865])
for aa=1:size(path_dir,2)
    subplot(6,2,aa)
    hold on
    ylim([0 15])
    xlabel('Distance [cm]')
    ylabel('Counts')
    title(num2str(aa))
    histogram(global_dist_cm{aa},0:5:100)
    txt = ['n=' num2str(size(global_dist_cm{aa},2))];
    text(80,9,txt)
end

%% Plot histogram 
%(this was generated to see the distribution between the A and B centers
%and decide how to split the data)
if 0
    figure
    subplot(2,1,1)
    hold on
    title('A&B Bin distance')
    histogram(cell2mat(combined_distances),30,'Normalization','probability')
    
    subplot(2,1,2)
    hold on
    title('A&B Bin distance metric')
    histogram(cell2mat(combined_dist_metric),30,'Normalization','probability')
end

%% Data for prism to get descriptive statistics

separation_cm_global = cell2mat(combined_final_global_dist).*1.96;


%% Calculate median of distribution and plot
median_dist = median(separation_cm_global);

separation_cm_global = separation_cm_global';

%% This is the histogram used in the supplement
f = figure;
set(f,'renderer','Painters')
set(f, 'DefaultAxesFontName', 'Arial');
set(f, 'DefaultTextFontName', 'Arial');
hold on
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
%export rasters fpr first colum

%export_fig(f ,fullfile('G:\Google_drive\task_selective_place_paper\input_figures_to_illustrator\Figure_3_figures','global_AB_dist.eps'))

%% Export place field separation data for global remappers
global_pf_dist.combined_final_global_dist = combined_final_global_dist;
global_pf_dist.median_dist = median_dist;

end

