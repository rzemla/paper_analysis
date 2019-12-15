function [combined_distances,combined_dist_metric] = place_AB_distance_separation(path_dir)

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
    
end

%% Export the data loaded in from all animals



%% Plot histogram
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

