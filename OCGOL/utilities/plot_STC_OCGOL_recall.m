function [outputArg1,outputArg2] = plot_STC_OCGOL_recall(animal_data, tunedLogical,registered,options)

%% Import variables


%% Get list of matching ROIs across sessions

%matched on all days
%matching_list = registered.multi.assigned_all;
%matched on any day (including nans)
%matching_list = registered.multi.assigned;
%use matches that were additionally manually filtered for mismatches
matching_list = registered.multi.assigned_filtered;

%% Define tuned cell combinations across trials

%make conditional here for si or ts tuned neurons
switch options.tuning_criterion
    case 'si' %spatial information
        %for each session
        for ss =options.sessionSelect%1:size(animal_data,2)
            %spatial information criterion
            Atuned{ss} = tunedLogical(ss).si.Atuned;
            Btuned{ss} = tunedLogical(ss).si.Btuned;
            
            AandB_tuned{ss} =  tunedLogical(ss).si.AandB_tuned;
            AorB_tuned{ss} = tunedLogical(ss).si.AorB_tuned;
            onlyA_tuned{ss} = tunedLogical(ss).si.onlyA_tuned;
            onlyB_tuned{ss} = tunedLogical(ss).si.onlyB_tuned;
        end
   case 'ts' %spatial information 
        for ss =options.sessionSelect%1:size(animal_data,2)
            %spatial information criterion
            Atuned{ss} = tunedLogical(ss).ts.Atuned;
            Btuned{ss} = tunedLogical(ss).ts.Btuned;
            
            AandB_tuned{ss} =  tunedLogical(ss).ts.AandB_tuned;
            AorB_tuned{ss} = tunedLogical(ss).ts.AorB_tuned;
            onlyA_tuned{ss} = tunedLogical(ss).ts.onlyA_tuned;
            onlyB_tuned{ss} = tunedLogical(ss).ts.onlyB_tuned;
        end
        
end

%% Extract mean STC map in each spatial bin (not normalized and not occupancy divided) (100 bins)
%for each session
for ss = options.sessionSelect%1:size(animal_data,2)
    A_df{ss} = animal_data{ss}.Place_cell{1}.Spatial_tuning_curve;
    B_df{ss} = animal_data{ss}.Place_cell{2}.Spatial_tuning_curve;
    
    A_df_both{ss} = A_df{ss}(:,AandB_tuned{ss});
    B_df_both{ss} = B_df{ss}(:,AandB_tuned{ss});
    
    A_df_onlyA{ss} = A_df{ss}(:,onlyA_tuned{ss});
    B_df_onlyA{ss} = B_df{ss}(:,onlyA_tuned{ss});
end

%% Compare session number (temporary)
%session to compare against
ses_comp = 2;

%% A vs. B on early vs late training (A or B tuned)

%ALL NEURONS in each session that meet criteria - tuned to either A or B
%early day
for ss =options.sessionSelect%1:size(animal_data,2)
    dF_maps_all_AB_early_late{ss} = [A_df{ss}(:,AorB_tuned{ss})', B_df{ss}(:,AorB_tuned{ss})'];
end

%sort each session by A map
for ss = options.sessionSelect %1:size(animal_data,2)
    %maxBin - spatial bin where activity is greatest for each ROI
    [~,maxBin_all_AB{ss}] = max(dF_maps_all_AB_early_late{ss}(:,1:100)', [], 1);
    %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
    [~,sortOrder_all_AB{ss}] = sort(maxBin_all_AB{ss},'ascend');
end

%plot side by side; day by day
figure;
subplot(2,1,1)
imagesc(dF_maps_all_AB_early_late{1}(sortOrder_all_AB{1},:))
title('5A5B')
hold on
colormap('jet')
caxis([0 1])
%A/B vertical separator line
plot([100 100],[1,size(dF_maps_all_AB_early_late{1},1)], 'k','LineWidth', 1.5);


hold off

subplot(2,1,2)
imagesc(dF_maps_all_AB_early_late{ses_comp}(sortOrder_all_AB{ses_comp},:))
hold on
title('Random AB')
colormap('jet')
caxis([0 1])
%A/B vertical separator line
plot([100 100],[1,size(dF_maps_all_AB_early_late{ses_comp},1)], 'k','LineWidth', 1.5);


hold off

%% Plot STCs from ROIs matching 1 and any chosen session thereafter
ses_comp = 7;

%get the matching ROIs from 2 ses
matching_ses_ROI_idxs = matching_list(:,[1,ses_comp]);
%get rid of nan values and get non-nan ROIs on both ses
matching_ses_ROI_idxs_nonan = matching_ses_ROI_idxs(find(sum(isnan(matching_ses_ROI_idxs),2)==0),:);

%generate matching A STC
A_comp_STCs = [A_df{1}(:,matching_ses_ROI_idxs_nonan(:,1))', A_df{ses_comp}(:,matching_ses_ROI_idxs_nonan(:,2))'];
%generate matching B STC
B_comp_STCs = [B_df{1}(:,matching_ses_ROI_idxs_nonan(:,1))', B_df{ses_comp}(:,matching_ses_ROI_idxs_nonan(:,2))'];

%sort by A on ses 1 and apply sort downstream
%sort by A trials on session 1
[~,matched_maxBin_A] = max(A_comp_STCs(:,1:100)', [], 1);
%[~,matched_maxBin_A] = max(A_comp_STCs(:,101:200)', [], 1);
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,matched_sortOrder_A] = sort(matched_maxBin_A,'ascend');

[~,matched_maxBin_B] = max(B_comp_STCs(:,1:100)', [], 1);
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,matched_sortOrder_B] = sort(matched_maxBin_B,'ascend');

%plot
figure('Position',[2630 50 490 890])
subplot(2,1,1)
imagesc(A_comp_STCs(matched_sortOrder_A,:))
hold on
colormap('jet')
caxis([0 1]);

subplot(2,1,2)
imagesc(B_comp_STCs(matched_sortOrder_B,:))
hold on
colormap('jet')
caxis([0 1]);

%% Make matching ROI list with tuning criteria for both sessions


%tuned to A or B on either sessions
AorB_idx{1} = find(AorB_tuned{1} ==1);
AorB_idx{ses_comp} = find(AorB_tuned{ses_comp} ==1);

%intersect with
[tuned_match_idx{1},match_idx{1},~] = intersect(matching_list(:,1),AorB_idx{1},'stable');
[tuned_match_idx{ses_comp},match_idx{ses_comp},~] = intersect(matching_list(:,ses_comp),AorB_idx{ses_comp},'stable');

%create not logical for nan exclusion from copied matrix assignement below
include_log{1} = false(1,size(matching_list,1));
include_log{1}(match_idx{1}) = 1; 
%session 2 
include_log{2} = false(1,size(matching_list,1));
include_log{2}(match_idx{ses_comp}) = 1;

%make copy
tuned_matching_ROI_list = matching_list;
%nan first session that are not tuned and last session that are not
%tuned
tuned_matching_ROI_list(~include_log{1},1) = nan;
tuned_matching_ROI_list(~include_log{2},ses_comp) = nan;

%which neurons to remove based on tuning criterion
keep_ROI = sum(isnan(tuned_matching_ROI_list),2) == 0;

%retain only tuned and matched ROIs
tuned_matching_ROI_list(~keep_ROI,:) = [];

%% Generate maps based on tuned matching ROI list
%row - session
%column - trial type

session_matched_tuned_dF_maps{1,1} = A_df{1}(:,tuned_matching_ROI_list(:,1))';
session_matched_tuned_dF_maps{1,2} = B_df{1}(:,tuned_matching_ROI_list(:,1))';
session_matched_tuned_dF_maps{2,1} = A_df{ses_comp}(:,tuned_matching_ROI_list(:,ses_comp))';
session_matched_tuned_dF_maps{2,2} = B_df{ses_comp}(:,tuned_matching_ROI_list(:,ses_comp))';

%combined 2x2
combined_maps_2x2 = cell2mat(session_matched_tuned_dF_maps);

%single matrix (ses 1 A, B, ses 2 A,B) - in 1 row
combined_maps_row = cell2mat(reshape(session_matched_tuned_dF_maps,1,4));

%sort by A trials on session 1
[~,matched_maxBin_1] = max(session_matched_tuned_dF_maps{1,1}(:,1:100)', [], 1);
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,matched_sortOrder_1] = sort(matched_maxBin_1,'ascend');

%sort by A trials on session 2
[~,matched_maxBin_2] = max(session_matched_tuned_dF_maps{2,1}(:,1:100)', [], 1);
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,matched_sortOrder_2] = sort(matched_maxBin_2,'ascend');

%session 1 above and session 2 below
figure;
subplot(2,2,1)
imagesc(session_matched_tuned_dF_maps{1,1}(matched_sortOrder_1,:))
colormap('jet')
caxis([0 1]);
subplot(2,2,2)
imagesc(session_matched_tuned_dF_maps{1,2}(matched_sortOrder_1,:))
colormap('jet')
caxis([0 1]);
subplot(2,2,3)
imagesc(session_matched_tuned_dF_maps{2,1}(matched_sortOrder_1,:))
colormap('jet')
caxis([0 1]);
subplot(2,2,4)
imagesc(session_matched_tuned_dF_maps{2,2}(matched_sortOrder_1,:))
colormap('jet')
caxis([0 1]);

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

