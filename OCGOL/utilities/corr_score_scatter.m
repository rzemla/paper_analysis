function [outputArg1,outputArg2] = corr_score_scatter(path_dir)


%% Load place cells defined according to correlation sig criteria (updated)
%correctly identified neurons will be in the final sub-struct
for ee=1:size(path_dir,2)
    load_remapping_ROI{ee} = fullfile(path_dir{ee},'cumul_analysis','remap_corr_idx.mat');
    remapping_ROI_data{ee} = load(string(load_remapping_ROI{ee}));
end

%% Load r-scores and p-values from all animals

%%%export the r values and p values associated with each rate map correlation
for ee=1:size(path_dir,2)
    load_corr_scores{ee} = fullfile(path_dir{ee},'cumul_analysis','tun_curve_corr.mat');
    r_scores_p_val{ee} = load(string(load_corr_scores{ee}));
end

%% Split the variables for remapping neuron idxs and r_scores/p_val into shorter names

%idxs of each class of remapping neurons
for ee=1:size(path_dir,2)
    common_idx{ee} = remapping_ROI_data{ee}.remapping_corr_idx.final.common;
    global_idx{ee} = remapping_ROI_data{ee}.remapping_corr_idx.final.global;
    partial_idx{ee} = remapping_ROI_data{ee}.remapping_corr_idx.final.partial;
    rate_idx{ee} = remapping_ROI_data{ee}.remapping_corr_idx.final.rate_remap_all;
end

%r values and p values for each class
for ee=1:size(path_dir,2)
    r_val{ee} = r_scores_p_val{ee}.tun_curve_corr.r;
    p_val{ee} = r_scores_p_val{ee}.tun_curve_corr.p_val;
    
end

%% Pair common, global, and partial with respective index from each animal

%r correlation values and p values for each of 4 classes

for ee=1:size(path_dir,2)
    %common
    r_common{ee} = r_val{ee}(common_idx{ee});
    p_common{ee} = p_val{ee}(common_idx{ee});
    
    %global
    r_global{ee} = r_val{ee}(global_idx{ee});
    p_global{ee} = p_val{ee}(global_idx{ee});
    
    %partial
    r_partial{ee} = r_val{ee}(partial_idx{ee});
    p_partial{ee} = p_val{ee}(partial_idx{ee});
    
    %rate
    r_rate{ee} = r_val{ee}(rate_idx{ee});
    p_rate{ee} = p_val{ee}(rate_idx{ee});
end

%% Merge into 1 vector r and p values

%global merge
r_global_merge = cell2mat(r_global);
p_global_merge = cell2mat(p_global);

%common merge
r_common_merge = cell2mat(r_common);
p_common_merge = cell2mat(p_common);

%partial merge
r_partial_merge = cell2mat(r_partial);
p_partial_merge = cell2mat(p_partial);

%rate merge
r_rate_merge = cell2mat(r_rate);
p_rate_merge = cell2mat(p_rate);


%% Plot scatter of correlation value vs -log(10) p value
%raise the values to the power of -1 to get the negative to the log10 when
%you use 'yscale' log

%darkorchid - common
%orchid - global

cmap_remapping_type = [147,112,219;...
                    255,0,255]./255;
% -log10 of p value
p_log10_thres = -log10(0.05);

marker_size =8;
%try log scale without logging the values
figure
%subplot(1,3,3)
axis square
hold on
xlim([-0.7 1.05])
ylim([-5 55])
%ylim([10^(-2) 10^60])
xlabel('Correlation score')
ylabel('-log_1_0(p)')
s1 = scatter(r_common_merge,-log10(p_common_merge),marker_size,'filled','MarkerFaceColor',cmap_remapping_type(1,:));
s2 = scatter(r_global_merge,-log10(p_global_merge),marker_size,'filled','MarkerFaceColor',cmap_remapping_type(2,:));
%plot p-threshold
plot([-1 2],[p_log10_thres, p_log10_thres],'k--','LineWidth',1)
%plot 0 correlation threshold
plot([0 0],[-5 55],'k--')
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)
%define legend
legend([s1 s2],{'Common','Global'},'location','northwest')


%scatter(r_common_merge,p_common_merge,marker_size,'filled','MarkerFaceColor',cmap_remapping_type(1,:))
%scatter(r_common_merge,p_common_merge.^(-1),marker_size,'filled','MarkerFaceColor',cmap_remapping_type(1,:))
%scatter(r_rate_merge,p_rate_merge.^(-1),marker_size,'filled','MarkerFaceColor',cmap_remapping_type(3,:))
%scatter(r_partial_merge,-log10(p_partial_merge),14,'filled')

%set(gca,'yscale','log')




% 
% figure
% hold on
% xlim([-0.7 1.1])
% %ylim([-2 50])
% xlabel('Correlation score')
% ylabel('-log10(p)')
% scatter(r_common_merge,(p_common_merge),14,'filled')
% scatter(r_global_merge,(p_global_merge),14,'filled')
% scatter(r_rate_merge,(p_rate_merge),14,'filled')
% 

end

