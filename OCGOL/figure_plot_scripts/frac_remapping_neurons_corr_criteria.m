function [frac_remapping] = frac_remapping_neurons_corr_criteria(path_dir)



%% Read relevant data - A&B tuned (TS)
for ee=1:size(path_dir,2)
    load_data_path_select_ROI{ee} = fullfile(path_dir{ee},'cumul_analysis','select_ROI_criteria.mat');
    select_ROI_data{ee} = load(string(load_data_path_select_ROI{ee}));
end

%% Load place cells defined according to correlation sig criteria (updated)
for ee=1:size(path_dir,2)
    load_remapping_ROI{ee} = fullfile(path_dir{ee},'cumul_analysis','remap_corr_idx.mat');
    remapping_ROI_data{ee} = load(string(load_remapping_ROI{ee}));
end

remapping_ROI_data{1, 1}.remapping_corr_idx.final.unclass

%% Get fraction relative to total A&B neurons for animal

%order: common, rate, global, partial, unclassified

%total remapping neurons for each animal (matrix)
for ee=1:size(path_dir,2)
    total_remap(ee,1) = size(remapping_ROI_data{ee}.remapping_corr_idx.final.common,1);
    total_remap(ee,2) = size(remapping_ROI_data{ee}.remapping_corr_idx.final.rate_remap_all,2);
    total_remap(ee,3) = size(remapping_ROI_data{ee}.remapping_corr_idx.final.global,1);
    total_remap(ee,4) = size(remapping_ROI_data{ee}.remapping_corr_idx.final.partial,2);
    total_remap(ee,5) = size(remapping_ROI_data{ee}.remapping_corr_idx.final.unclass,2);
    %total_remap(ee,6) = size(select_ROI_data{ee}.task_remapping_ROIs.mixed,2);
end

%total counts of each class of neurons
total_remap_count = sum(total_remap,2);
%fraction of total for each animal
fraction_remap = total_remap./total_remap_count;
%total of each class of neurons
total_class_count = sum(total_remap,1);
%total neurons
total_neurons = sum(sum(total_remap,1),2);

%% Mean, std, and sem of classes above
class_mean = mean(fraction_remap,1);
class_std = std(fraction_remap,[],1);
class_sem = class_std./sqrt(size(total_remap,1));

%% Export data

frac_remapping.class_names = {'Common','Activity','Global','Partial', 'Unclassified'};
frac_remapping.class_mean = class_mean;
frac_remapping.class_sem = class_sem;

%% Plot bars

figure()
hold on
ylabel('Fraction of total remapping neurons')
bar(class_mean,'FaceColor',[139, 0, 139]/255)
%add significance bars
xticks([1:5])
xticklabels({'Common','Rate','Global','Partial','Unclassified'})
%sigstar({[1,2], [1,3],[1,4]})
xtickangle(45)
set(gca,'FontSize',16)
errorbar([1:5],class_mean,class_sem,'k.')
ylim([0 0.4])


%% Rearrange order to match that of color scheme below
%from: near, far rate, common, partial, mixed
%to: common, rate, near, far, partial,mixed
% class_mean_ordered = class_mean;
% class_std_ordered = class_std;
% class_sem_ordered = class_sem;
% %mean re-order
% class_mean_ordered(1) = class_mean(4);
% class_mean_ordered(2) = class_mean(3);
% class_mean_ordered(3) = class_mean(1);
% class_mean_ordered(4) = class_mean(2);
% %std re-order
% class_std_ordered(1) = class_std(4);
% class_std_ordered(2) = class_std(3);
% class_std_ordered(3) = class_std(1);
% class_std_ordered(4) = class_std(2);
% %sem re-order
% class_sem_ordered(1) = class_sem(4);
% class_sem_ordered(2) = class_sem(3);
% class_sem_ordered(3) = class_sem(1);
% class_sem_ordered(4) = class_sem(2);

%% Color scheme for classes
% %forest green - common
% [34,139,34];
% %rate - orange red
% [255,69,0];
% %light cyan (medium turquoise) - near remap
% [72,209,204];
% %dark cyan - global remap
% [0 139 139];
% %partial - metallic gold
% [212,175,55];
% %dark purple (indigo) - unclassified
% [75,0,130];
% 
% %combine into one matrix
% cmap = [34,139,34; 255,69,0; 72,209,204; 0 139 139; 212,175,55; 75,0,130];

%order common - rate - near - far - partial - mixed

%% Perform Kruskall Wallis - refer to prism for this

end

