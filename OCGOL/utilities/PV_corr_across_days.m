function [outputArg1,outputArg2] = PV_corr_across_days(animal_data, tunedLogical,registered,options)


%% Get list of matching ROIs across sessions

%matched on all days
%matching_list = registered.multi.assigned_all;
%matched on any day (including nans)
%matching_list = registered.multi.assigned;
%use matches that were additionally manually filtered for mismatches
matching_list = registered.multi.assigned_filtered;


%% Extract mean STC map in each spatial bin (not normalized and not occupancy divided) (100 bins)
%for each session
for ss = options.sessionSelect%1:size(animal_data,2)
    %normalized to A or B trials independently
    A_STC{ss} = animal_data{ss}.Place_cell{1}.Spatial_tuning_curve;
    B_STC{ss} = animal_data{ss}.Place_cell{2}.Spatial_tuning_curve;
    
    %A_STC{ss} =
    %B_STC{ss} = 
    
end

%% STC of neurons matching across all sessions

x=1


%% Plot STCs from ROIs matching 1 and any chosen session thereafter
ses_comp = 7;

%get the matching ROIs from 2 ses
matching_ses_ROI_idxs = matching_list(:,[1,ses_comp]);
%get rid of nan values and get non-nan ROIs on both ses
matching_ses_ROI_idxs_nonan = matching_ses_ROI_idxs(find(sum(isnan(matching_ses_ROI_idxs),2)==0),:);

%generate matching A STC
A_comp_STCs = [A_STC{1}(:,matching_ses_ROI_idxs_nonan(:,1))', A_STC{ses_comp}(:,matching_ses_ROI_idxs_nonan(:,2))'];
%generate matching B STC
B_comp_STCs = [B_STC{1}(:,matching_ses_ROI_idxs_nonan(:,1))', B_STC{ses_comp}(:,matching_ses_ROI_idxs_nonan(:,2))'];

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


end

