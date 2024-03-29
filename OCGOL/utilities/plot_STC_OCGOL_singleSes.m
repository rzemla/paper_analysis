function [outputArg1,outputArg2] = plot_STC_OCGOL_singleSes(animal_data, tunedLogical,options)

%% Import variables

%% Define tuned cell combinations across trials

%make conditional here for si or ts tuned neurons
switch options.tuning_criterion
    case 'si' %spatial information
        %for each session
        for ss =1:size(animal_data,2)
            %spatial information criterion
            Atuned{ss} = tunedLogical(ss).si.Atuned;
            Btuned{ss} = tunedLogical(ss).si.Btuned;
            
            AandB_tuned{ss} =  tunedLogical(ss).si.AandB_tuned;
            AorB_tuned{ss} = tunedLogical(ss).si.AorB_tuned;
            onlyA_tuned{ss} = tunedLogical(ss).si.onlyA_tuned;
            onlyB_tuned{ss} = tunedLogical(ss).si.onlyB_tuned;
            
            neither_tuned{ss} = tunedLogical(ss).si.neither;
        end
   case 'ts' %spatial information 
        for ss =1:size(animal_data,2)
            %spatial information criterion
            Atuned{ss} = tunedLogical(ss).ts.Atuned;
            Btuned{ss} = tunedLogical(ss).ts.Btuned;
            
            AandB_tuned{ss} =  tunedLogical(ss).ts.AandB_tuned;
            AorB_tuned{ss} = tunedLogical(ss).ts.AorB_tuned;
            onlyA_tuned{ss} = tunedLogical(ss).ts.onlyA_tuned;
            onlyB_tuned{ss} = tunedLogical(ss).ts.onlyB_tuned;
        end
    case 'selective_filtered'
        %neurons idxs associated with selective filtering for
        %task-selectivity
        select_filt_ROIs.A = task_selective_ROIs.A.idx;
        select_filt_ROIs.B = task_selective_ROIs.B.idx;
        
end

%% Extract mean STC map in each spatial bin (not normalized and not occupancy divided) (100 bins)
%for each session
%correct only
for ss =1:size(animal_data,2)
    A_df{ss} = animal_data{ss}.Place_cell{1}.Spatial_tuning_curve;
    B_df{ss} = animal_data{ss}.Place_cell{2}.Spatial_tuning_curve;
    
    %Gs smoothed, but not normalized (nn) to itself
    A_STC_nn{ss} = animal_data{ss}.Place_cell{1}.Spatial_Info.rate_map_smooth{8};
    B_STC_nn{ss} = animal_data{ss}.Place_cell{2}.Spatial_Info.rate_map_smooth{8};
    
    %A_df_both{ss} = A_df{ss}(:,AandB_tuned{ss});
    %B_df_both{ss} = B_df{ss}(:,AandB_tuned{ss});
    
    %A_df_onlyA{ss} = A_df{ss}(:,onlyA_tuned{ss});
    %B_df_onlyA{ss} = B_df{ss}(:,onlyA_tuned{ss});
end

%% Normalize eahc STC ROI across both trials in non-norm STCs

for ss =1:size(animal_data,2)
        %get max value for each ROIs between trials
        max_STC_across_trials{ss} = max([A_STC_nn{ss};B_STC_nn{ss}]);
        %normalize each matrix to these values (tn = trial normalized)
        A_STC_tn{ss} = A_STC_nn{ss}./max_STC_across_trials{ss};
        B_STC_tn{ss} = B_STC_nn{ss}./max_STC_across_trials{ss};
    %for max value to normalize to by 1
    %in future, do normalization range (0,1)
end

%% A vs. B on early vs late training (A or B tuned)

%ALL NEURONS in each session that meet criteria - tuned to either A or B
%early day
for ss =1:size(animal_data,2)
    STC_norm_self_AB{ss} = [A_df{ss}(:,AorB_tuned{ss})', B_df{ss}(:,AorB_tuned{ss})'];
    STC_norm_trials_AB{ss} = [A_STC_tn{ss}(:,AorB_tuned{ss})', B_STC_tn{ss}(:,AorB_tuned{ss})'];
    %STC_norm_trials_AB{ss} = [A_STC_tn{ss}(:,neither_tuned{ss})', B_STC_tn{ss}(:,neither_tuned{ss})'];
end

%sort each session by A map
for ss =1:size(animal_data,2)
    %maxBin - spatial bin where activity is greatest for each ROI
    [~,maxBin_all_AB{ss}] = max(STC_norm_trials_AB{ss}(:,1:100)', [], 1);
    %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
    [~,sortOrder_all_AB{ss}] = sort(maxBin_all_AB{ss},'ascend');
end

%plot side by side; day by day
figure;
%subplot(2,1,1)
imagesc(STC_norm_self_AB{1}(sortOrder_all_AB{1},:))
%title('5A5B')
hold on
colormap('jet')
caxis([0 1])
%A/B vertical separator line
plot([100 100],[1,size(STC_norm_self_AB{1},1)], 'k','LineWidth', 1.5);

%STC normalized across trials
f= figure;
%subplot(2,1,1)
imagesc(STC_norm_trials_AB{1}(sortOrder_all_AB{1},:))
%title('Normalized according to trials')
hold on
caxis([0 1])
colormap('jet')
cbar= colorbar;
cbar.Label.String = 'Normalized activity';
cbar.Ticks = [0 0.5 1];
ax1 = gca;
ylabel('Neuron #');
xlabel('Normalized position');
ax1.XTick = [1 100 200];
ax1.XTickLabel = {'0','1','1'};
%A/B vertical separator line
plot([100 100],[1,size(STC_norm_trials_AB{1},1)], 'k','LineWidth', 1.5);

hold off

% subplot(2,1,2)
% imagesc(dF_maps_all_AB_early_late{2}(sortOrder_all_AB{2},:))
% hold on
% title('Random AB')
% colormap('jet')
% caxis([0 1])
% %A/B vertical separator line
% plot([100 100],[1,size(dF_maps_all_AB_early_late{2},1)], 'k','LineWidth', 1.5);


%hold off

%PVcorr = corr(A_STC_nn{session_nb}(:,tuning_selection{session_nb}),B_STC_nn{session_nb}(:,tuning_selection{session_nb}))
PVcorr = corr(A_STC_nn{1}',B_STC_nn{1}');


figure('Position',[1350, 90, 500 860]);
subplot(2,1,1)
imagesc(PVcorr)
hold on
title('Population vector correlation');
xlabel('Spatial bin')
ylabel('Spatial bin')
colormap('jet')
caxis([0 1])
xticks([20 40 60 80 100]);
yticks([20 40 60 80 100]);
axis('square')
cbar2 = colorbar;
cbar2.Label.String = 'Correlation coefficient';
cbar2.Ticks = [0 0.5 1];

subplot(2,1,2)
hold on
title('Population vector correlation');
plot(diag(PVcorr),'k','LineWidth',1.5)
xlabel('Spatial bin')
ylabel('Correlation coef.');
plot([30 30],[0 1],'r--','LineWidth',1.5);
text([31 31],[0.9 0.9],'Reward zone B','Color','r')
plot([70 70],[0 1],'b--','LineWidth',1.5);
text([71 71],[0.3 0.3],'Reward zone A','Color','b')
plot([10 10],[0 1],'g--','LineWidth',1.5);
text([11 11],[0.9 0.9],'Odor zone\newline end','Color','g')

%% Make matching ROI list with tuning criteria for both sessions

%tuned to A or B on either sessions
% AorB_idx{1} = find(AorB_tuned{1} ==1);
% AorB_idx{2} = find(AorB_tuned{2} ==1);
% 
% %intersect with
% [tuned_match_idx{1},match_idx{1},~] = intersect(matching_list(:,1),AorB_idx{1},'stable');
% [tuned_match_idx{2},match_idx{2},~] = intersect(matching_list(:,2),AorB_idx{2},'stable');
% 
% %create not logical for nan exclusion from copied matrix assignement below
% include_log{1} = false(1,size(matching_list,1));
% include_log{1}(match_idx{1}) = 1; 
% %session 2 
% include_log{2} = false(1,size(matching_list,1));
% include_log{2}(match_idx{2}) = 1;
% 
% %make copy
% tuned_matching_ROI_list = matching_list;
% %nan first session that are not tuned and last session that are not
% %tuned
% tuned_matching_ROI_list(~include_log{1},1) = nan;
% tuned_matching_ROI_list(~include_log{2},2) = nan;
% 
% %which neurons to remove based on tuning criterion
% keep_ROI = sum(isnan(tuned_matching_ROI_list),2) == 0;
% 
% %retain only tuned and matched ROIs
% tuned_matching_ROI_list(~keep_ROI,:) = [];

%% Generate maps based on tuned matching ROI list
%row - session
%column - trial type
% 
% session_matched_tuned_dF_maps{1,1} = A_df{1}(:,tuned_matching_ROI_list(:,1))';
% session_matched_tuned_dF_maps{1,2} = B_df{1}(:,tuned_matching_ROI_list(:,1))';
% session_matched_tuned_dF_maps{2,1} = A_df{2}(:,tuned_matching_ROI_list(:,2))';
% session_matched_tuned_dF_maps{2,2} = B_df{2}(:,tuned_matching_ROI_list(:,2))';

% %combined 2x2
% combined_maps_2x2 = cell2mat(session_matched_tuned_dF_maps);
% 
% %single matrix (ses 1 A, B, ses 2 A,B) - in 1 row
% combined_maps_row = cell2mat(reshape(session_matched_tuned_dF_maps,1,4));
% 
% %sort by A trials on session 1
% [~,matched_maxBin_1] = max(session_matched_tuned_dF_maps{1,1}(:,1:100)', [], 1);
% %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
% [~,matched_sortOrder_1] = sort(matched_maxBin_1,'ascend');
% 
% %sort by A trials on session 2
% [~,matched_maxBin_2] = max(session_matched_tuned_dF_maps{2,1}(:,1:100)', [], 1);
% %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
% [~,matched_sortOrder_2] = sort(matched_maxBin_2,'ascend');

%session 1 above and session 2 below
% figure;
% subplot(2,2,1)
% imagesc(session_matched_tuned_dF_maps{1,1}(matched_sortOrder_1,:))
% colormap('jet')
% caxis([0 1]);
% subplot(2,2,2)
% imagesc(session_matched_tuned_dF_maps{1,2}(matched_sortOrder_1,:))
% colormap('jet')
% caxis([0 1]);
% subplot(2,2,3)
% imagesc(session_matched_tuned_dF_maps{2,1}(matched_sortOrder_1,:))
% colormap('jet')
% caxis([0 1]);
% subplot(2,2,4)
% imagesc(session_matched_tuned_dF_maps{2,2}(matched_sortOrder_1,:))
% colormap('jet')
% caxis([0 1]);

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

