%% Define experiment core directories
%screened and good datasets for analysis
%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_120618_1'}; %
%field rate error (after PSEM silencing
path_dir{1} = {'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1'};
path_dir{2} = {'G:\Figure_2_3_selective_remap\I42R_AB_d1_032118_1'};

path_dir{3} = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1'};
path_dir{4} = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_2'};

path_dir{5} = {'G:\Figure_2_3_selective_remap\I53LT_AB_sal_113018_1'}; 
path_dir{6} = {'G:\Figure_2_3_selective_remap\I56_RTLS_AB_prePost_sal_042419_1'};

path_dir{7} = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_113018_1'};
path_dir{8} = {'G:\Figure_2_3_selective_remap\I57_RTLS_AB_prePost_792_042519_1'};

path_dir{9} = {'G:\Figure_2_3_selective_remap\I45_RT_AB_d1_062018_1'};
path_dir{10} = {'G:\Figure_2_3_selective_remap\I46_AB_d1_062018_1'};
path_dir{11} = {'G:\Figure_2_3_selective_remap\I57_LT_ABrand_no_punish_042119_1'};

%path_dir = {'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I42R_AB_d1_032118_1'};

%path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_2'};
 
%path_dir = {'G:\Figure_2_3_selective_remap\I53LT_AB_sal_113018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I56_RTLS_AB_prePost_sal_042419_1'};

%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_113018_1'};  
%path_dir = {'G:\Figure_2_3_selective_remap\I57_RTLS_AB_prePost_792_042519_1'};

%path_dir = {'G:\Figure_2_3_selective_remap\I45_RT_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I46_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I57_LT_ABrand_no_punish_042119_1'};

%N =10; FOV = 11 (2 FOV from single animal)

%% Load the reward endzones calculated from script below (remapping_centroids_updated) into variable

%load in struct with reward zone start and end for each animal
load(fullfile('G:\Figure_2_3_selective_remap\cumulative_data_output','reward_zones_all_animals.mat'));


%% Source data export struct - Fig 2 data

source_data_task_sel_remap = struct();

%% Fraction of place cell (A sel, B sel, A&B remapping,neither) - bar chart and pie chart
%Figure 2C
[fraction_place_plot_data,source_data_task_sel_remap] = fraction_place_cells(path_dir,source_data_task_sel_remap);

%% A/B/A&B/N cell SI/TS score distributions plots (Figure 2C/D)
%export these values and feed into master plotter

si_ts_score_dist_data = si_ts_score_distributions(path_dir);

%% Cumulative task-selective STCs

task_sel_STC_data = cumulative_task_sel_STC(path_dir);

%% Centroid distribution for A or B selective neurons - 
%how many bins to display
%1 - 20 (every 5) relative to 100 bins
%2 - 25 (every 4) relative to 100 bins
%3 - 50 (every 2) relative to 100 bins
options.bin_choose = 2;

[centroid_dist_data,source_data_task_sel_remap] = ...
    centroid_dist_selective(path_dir,reward_zones_all_animal,options,source_data_task_sel_remap);

%% AUC/min scatterplots of A vs B neurons for each animal 

[AUC_data,source_data_task_sel_remap] = auc_scatterplots(path_dir,source_data_task_sel_remap);

%% PV/TC correlation

[tc_corr_sel_data,source_data_task_sel_remap] = tc_pv_correlation_task_sel(path_dir,source_data_task_sel_remap);

%% Master plotter for Figure 2

fig2_master_plotter(fraction_place_plot_data,si_ts_score_dist_data,...
        task_sel_STC_data,centroid_dist_data,AUC_data, tc_corr_sel_data);

%% Export source data Figure 2

%navigate to matlab summary stats directory
cd('G:\Google_drive\task_selective_place_paper\matlab_data')

save('source_data_fig2.mat', 'source_data_task_sel_remap','-v7.3');


%% Source data export struct - Fig 3

source_data_AB_remap = struct();

%% Source data export struct - supplement data/figs for Fig 2 and Fig 3

source_data_sup_2_3 = struct();

%% FIGURE 3
%% Place field A&B distance separation (global remapper distance distribution - supplement)

[combined_distances,combined_dist_metric,global_pf_dist] = place_AB_distance_separation(path_dir);

%export the lap speed data and event speed data for cumulative analysis
save(fullfile('G:\Figure_2_3_selective_remap\cumulative_data_output','pf_distances_A&B_neurons.mat'),'combined_distances', 'combined_dist_metric');

%distances between global remapping neurons for sup figure data
source_data_sup_2_3.global_pf_dist= global_pf_dist;

%% Remapping ROIs based on correlation map criteria and return quantile for common field separation (Fig 3c)
%STCs of the different categories of neurons are here (Figure 3C)
%use this as definition for the separation of fields (common and partial
%field) for partially remapping neurons
[common_cutoffs_95,remap_rate_maps] = remapping_corr(path_dir);

save(fullfile('G:\Figure_2_3_selective_remap\cumulative_data_output','common_cutoff.mat'),'common_cutoffs_95');

%% Generate mean/SEM traces for rate remapping neurons - Figure 3 rate remapping (supplement activity remap and Fig 3c data)
%also calculate A-B/A+B index and compare against common neurons
%extract supplement data
[activity_remap,remap_sup_plot_data] = rate_remap_traces(path_dir);

%% Activity remap supplement figure

fig_sup_activity_remap_master(remap_sup_plot_data)

%% Fractional remapping using updated criteria - fraction of each class

[frac_remapping,source_data_AB_remap] =frac_remapping_neurons_corr_criteria(path_dir,source_data_AB_remap);

%% Generate scatterplot with correlation against p-value for supplement (supplement data)
%return global correlation scores for each animal/FOV

[r_global] = corr_score_scatter(path_dir);

%export this value
save(fullfile('G:\Figure_2_3_selective_remap\cumulative_data_output','r_global.mat'),'r_global');

%% Modification to original code (updated using established criteria - majority of Figure 3 data is here)
%extract scatterplot points for for global remappers from here

%return bin position of the respective reward zones
[reward_zones_all_animal,partial_idx_by_animal_zone,remap_prop_figs,global_dist_scatter,source_data_AB_remap] = ...
remapping_centroids_updated(path_dir,source_data_AB_remap);

save(fullfile('G:\Figure_2_3_selective_remap\cumulative_data_output','partial_idx_by_animal.mat'),'partial_idx_by_animal_zone');

%save the reward zone endpoints (100 bins) into a common file for use by functions
%above
save(fullfile('G:\Figure_2_3_selective_remap\cumulative_data_output','reward_zones_all_animals.mat'),'reward_zones_all_animal');


%scatterplot data for global remappers (bin position (100)
source_data_sup_2_3.global_dist_scatter = global_dist_scatter;


%% Export source data Figure 3

%navigate to matlab summary stats directory
cd('G:\Google_drive\task_selective_place_paper\matlab_data')

save('source_data_fig3.mat', 'source_data_AB_remap','-v7.3');



%% Master plotter for Figure 3
%load frac remapping data
fig3_master_plotter(remap_rate_maps,activity_remap,frac_remapping,remap_prop_figs)

%% Event vs speed scatterplot for task selective neurons for supplement (Sup 3 - distributions) (Ex fig 4)

event_speed_plot = event_speed_scatter(path_dir);

%% Place field analysis (width and number for selective neurons) - supplement place field data (Ex fig 4)

pf_prop_data = place_field_analysis(path_dir);


%% Supplement Figure 4 speed of animal in place cells and place field props

fig_sup_speed_place_field_master(event_speed_plot,pf_prop_data)

%% Supplement plotter global remappers

fig_sup_global_remap_master(global_pf_dist,global_dist_scatter,reward_zones_all_animal)

%% Speed analysis for each animal used in Figure 2/3 (mean bin speed on A/B laps) (Sup 2)

lap_speed_by_animal(path_dir)

%save(fullfile(path_dir{1},'cumul_analysis','lap_and_event_speed.mat'),'mean_bin_speed', 'lap_bin_split','mean_event_speed');


%% Centroid difference for A&B tuned neurons and centroid diff as fxn of max bin (OLD)
%scatterplot of centroid difference as a function of center between
%centroid of max place field - not used
%centroid_difference(path_dir)

%% Remapping centroids (OLD)
%Figure 3E/F ( all Figure 3 code is here - organize this)
%options.lowPVcorr = [6 7 8];
%organize this 
%remapping_centroids(path_dir,options)

%% Fractional distribution of remapping neuron subtypes - Figure 3C (OLD)

%frac_remapping_neurons(path_dir)
