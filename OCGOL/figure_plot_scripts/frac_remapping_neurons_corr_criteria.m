function [outputArg1,outputArg2] = frac_remapping_neurons_corr_criteria(path_dir)



%% Read relevant data - A&B tuned (TS)
for ee=1:size(path_dir,2)
    load_data_path_select_ROI{ee} = fullfile(path_dir{ee},'cumul_analysis','select_ROI_criteria.mat');
    select_ROI_data{ee} = load(string(load_data_path_select_ROI{ee}));
end

%% Get fraction relative to total A&B neurons for animal

%order: common, rate, global, partial, unclassified

%total remapping neurons for each animal (matrix)
for ee=1:size(path_dir,2)
    total_remap(ee,1) = size(select_ROI_data{ee}.task_remapping_ROIs.global_near,2);
    total_remap(ee,2) = size(select_ROI_data{ee}.task_remapping_ROIs.global_far,2);
    total_remap(ee,3) = size(select_ROI_data{ee}.task_remapping_ROIs.rate,2);
    total_remap(ee,4) = size(select_ROI_data{ee}.task_remapping_ROIs.common,2);
    total_remap(ee,5) = size(select_ROI_data{ee}.task_remapping_ROIs.partial,2);
    total_remap(ee,6) = size(select_ROI_data{ee}.task_remapping_ROIs.mixed,2);
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

%% Rearrange order to match that of color scheme below
%from: near, far rate, common, partial, mixed
%to: common, rate, near, far, partial,mixed
class_mean_ordered = class_mean;
class_std_ordered = class_std;
class_sem_ordered = class_sem;
%mean re-order
class_mean_ordered(1) = class_mean(4);
class_mean_ordered(2) = class_mean(3);
class_mean_ordered(3) = class_mean(1);
class_mean_ordered(4) = class_mean(2);
%std re-order
class_std_ordered(1) = class_std(4);
class_std_ordered(2) = class_std(3);
class_std_ordered(3) = class_std(1);
class_std_ordered(4) = class_std(2);
%sem re-order
class_sem_ordered(1) = class_sem(4);
class_sem_ordered(2) = class_sem(3);
class_sem_ordered(3) = class_sem(1);
class_sem_ordered(4) = class_sem(2);

%% Color scheme for classes
%forest green - common
[34,139,34];
%rate - orange red
[255,69,0];
%light cyan (medium turquoise) - near remap
[72,209,204];
%dark cyan - global remap
[0 139 139];
%partial - metallic gold
[212,175,55];
%dark purple (indigo) - unclassified
[75,0,130];

%combine into one matrix
cmap = [34,139,34; 255,69,0; 72,209,204; 0 139 139; 212,175,55; 75,0,130];

%order common - rate - near - far - partial - mixed

%% Plot bars

figure('Position',[2704 336 641 500])
hold on
ylabel('Fraction of total remapping neurons')
bar(class_mean_ordered,'FaceColor',[139, 0, 139]/255)
xticks([1:6])
xticklabels({'Common','Rate','Global near','Global far','Partial','Mixed'})
xtickangle(45)
set(gca,'FontSize',16)
errorbar([1:6],class_mean_ordered,class_sem_ordered,'k.')
ylim([0 0.4])


end

