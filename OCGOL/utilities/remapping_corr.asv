function [outputArg1,outputArg2] = remapping_corr(path_dir)

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

%% Load in STCs
%load(fullfile(path_dir{1},'cumul_analysis','place_field_AB_distances.mat'),'ts_bin_conv_diff','ts_bin_conv_diff');

%% Combined into single matrix
for ss=1:size(remap_corr_idx,2)
    non_global_idx{ss} = remap_corr_idx{ss}.remapping_corr_idx.non_global_idx;
    global_idx{ss} = remap_corr_idx{ss}.remapping_corr_idx.global_idx;
    rate_only_idx{ss} = remap_corr_idx{ss}.remapping_corr_idx.rate_remap_grp_only;
    rate_part_idx{ss} = remap_corr_idx{ss}.remapping_corr_idx.rate_remap_all;
    common_idx{ss} = setdiff(non_global_idx{ss},rate_only_idx{ss});
end

%plot the STCs associated with each class
non_global_comb = cell2mat(non_global_idx');
global_comb = cell2mat(global_idx');
rate_only_comb = cell2mat(rate_only_idx);

%common neurons 
common_comb = cell2mat(common_idx');

%total neurons from all animals
total_ct = length(non_global_comb) + length(global_comb);

%fraction global
length(global_comb)./total_ct;

%save the neurons parsed by correlation remapping criteria
%save(fullfile(path_dir{1},'cumul_analysis','remap_corr_idx.mat'),'remapping_corr_idx'); 

%% Extract each class STC
for aa=1:size(path_dir,2)
    %A vs B side by side
    STC_tn_global{aa} = [STC{aa}.STC_export.A_STC_tn{1}(:,global_idx{aa})', STC{aa}.STC_export.B_STC_tn{1}(:,global_idx{aa})'];
    STC_tn_common{aa} = [STC{aa}.STC_export.A_STC_tn{1}(:,common_idx{aa})', STC{aa}.STC_export.B_STC_tn{1}(:,common_idx{aa})'];
    STC_tn_rate_only{aa} = [STC{aa}.STC_export.A_STC_tn{1}(:,rate_only_idx{aa})', STC{aa}.STC_export.B_STC_tn{1}(:,rate_only_idx{aa})'];
end

%% Sort STCs for each session
%global sort
for aa=1:size(path_dir,2)
    %sort global
    %input: ROI x bin (concatenated from both A and B trials)
    sortOrder_global{aa} = sortSTC(STC_tn_global{aa},1:100);
    
    sortOrder_common{aa} = sortSTC(STC_tn_common{aa},1:100);
    
    sortOrder_rate_only{aa} = sortSTC(STC_tn_rate_only{aa},1:100);
    
end

%sort each global STC
for aa=1:size(path_dir,2)
    STC_tn_global_sorted{aa} = STC_tn_global{aa}(sortOrder_global{aa},:);
    
    STC_tn_common_sorted{aa} = STC_tn_common{aa}(sortOrder_common{aa},:);
    
    STC_tn_rate_sorted{aa} = STC_tn_rate_only{aa}(sortOrder_rate_only{aa},:);
    
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


%% Plot individually sorted global, common, rate STCs

figure
imagesc(cell2mat(STC_tn_common_sorted'))
hold on
title('Common')
caxis([0 1])
colormap('jet')

figure
imagesc(cell2mat(STC_tn_global_sorted'))
hold on
title('Global remappers')
caxis([0 1])
colormap('jet')

figure
imagesc(cell2mat(STC_tn_rate_sorted'))
hold on
title('Rate')
caxis([0 1])
colormap('jet')

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

end

