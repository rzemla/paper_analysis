function [outputArg1,outputArg2] = remapping_corr(path_dir)

%% Load in place field distance data

for aa=1:size(path_dir,2)
    load_data_path{aa} = fullfile(path_dir{aa},'cumul_analysis','remap_corr_idx.mat');
    remap_corr_idx{aa} = load(string(load_data_path{aa}));
end

%load(fullfile(path_dir{1},'cumul_analysis','place_field_AB_distances.mat'),'ts_bin_conv_diff','ts_bin_conv_diff');

%% Combined into single matrix
for ss=1:size(remap_corr_idx,2)
    non_global_idx{ss} = remap_corr_idx{ss}.remapping_corr_idx.non_global_idx;
    global_idx{ss} = remap_corr_idx{ss}.remapping_corr_idx.global_idx;
end

non_global_comb = cell2mat(non_global_idx');
global_comb = cell2mat(global_idx');

%total neurons from all animals
total_ct = length(non_global_comb) + length(global_comb);

%fraction global
length(global_comb)./total_ct;

%save the neurons parsed by correlation remapping criteria
%save(fullfile(path_dir{1},'cumul_analysis','remap_corr_idx.mat'),'remapping_corr_idx'); 

end

