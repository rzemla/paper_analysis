function [outputArg1,outputArg2] = plot_STC_transient_rate(animal_data, tunedLogical,registered,field_event_rates, pf_vector,options)

%% Import variables


%% Get list of matching ROIs across sessions

matching_list = registered.multi.assigned_all;

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
        
end

%% Extract mean STC map in each spatial bin
%for each session
for ss =1:size(animal_data,2)
    %(normalized to neuron peak heigh, Gs smoothed (100 bins)
    A_STC{ss} = animal_data{ss}.Place_cell{4}.Spatial_tuning_curve;
    B_STC{ss} = animal_data{ss}.Place_cell{5}.Spatial_tuning_curve;
    
    %Gs smoothed, but not normalized (nn) to itself
    A_STC_nn{ss} = animal_data{ss}.Place_cell{4}.Spatial_Info.rate_map_smooth{8};
    B_STC_nn{ss} = animal_data{ss}.Place_cell{5}.Spatial_Info.rate_map_smooth{8};
    
    A_STC_both{ss} = A_STC{ss}(:,AandB_tuned{ss});
    B_STC_both{ss} = B_STC{ss}(:,AandB_tuned{ss});
    
    A_STC_onlyA{ss} = A_STC{ss}(:,onlyA_tuned{ss});
    B_STC_onlyA{ss} = B_STC{ss}(:,onlyA_tuned{ss});
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

%% Create a sort lists based on transient rate
%for each each peak find the idx of the one with highest in-field transient rate
for ss=1:size(animal_data,2)
    %for each trial (A or B) regardless if correct
    for tt=4:5
        %for each ROI
        for rr =1:size(field_event_rates{ss}{tt},2)
            %if more than 1 field
            if size(field_event_rates{ss}{tt}{rr},2) > 1
                [~,max_transient_peak{ss}{tt}(rr)] = max(field_event_rates{ss}{tt}{rr});

            elseif size(field_event_rates{ss}{tt}{rr},2) == 1 %if 1
                max_transient_peak{ss}{tt}(rr) = 1;
            else
                max_transient_peak{ss}{tt}(rr) = nan;
            end
        end
    end
end

%% Get centroid for each neuron (bin of centroid)

%two approaches (use 1 for now) - can set up as conditional option
%1) use global vector for single field and place field vector
%for neurons with more than 1 vector
%2) use max field vectors regardless if 1 or multiple fields

%get tuning vectors from tuning specificity calculation
for ss=1:size(animal_data,2)
    %for each trial (A or B) regardless if correct
    for tt=4:5
        %tuning vectors for each ROI
        tun_vectors{ss}{tt} = animal_data{ss}.Place_cell{tt}.Tuning_Specificity.tuning_vector;
        %original sum vector
        tun_vector{ss}{tt} = animal_data{ss}.Place_cell{tt}.Tuning_Specificity.tuning_vector_specificity;
    end
end

%convert centroid to bin location
for ss=1:size(animal_data,2)
    %for each trial (A or B) regardless if correct
    for tt=4:5
        %for each ROI
        for rr=1:size(pf_vector{ss}{tt},2)
            %if not empty
            if ~isempty(pf_vector{ss}{tt}{rr})
                %if single field, use original vector
                if size(pf_vector{ss}{tt}{rr},2) == 1
                    sort_vector{ss}{tt}(rr) = tun_vector{ss}{tt}(rr);
                else %if more than 1
                    sort_vector{ss}{tt}(rr) = pf_vector{ss}{tt}{rr}(max_transient_peak{ss}{tt}(rr));
                end
            else
                sort_vector{ss}{tt}(rr) = nan;
            end
        end
    end
end

%get the bin location of sort tuning vector
%get degree angle of each and classify in 1-100 bin
for ss=1:size(animal_data,2)
    %for each trial (A or B) regardless if correct
    for tt=4:5
        sort_vec_angle{ss}{tt} = angle(sort_vector{ss}{tt});
    end
end

for ss=1:size(animal_data,2)
    %for each trial (A or B) regardless if correct
    for tt=4:5
        %convert rad 2 deg and assign spatial bins
        deg_angles{ss}{tt} = rad2deg(sort_vec_angle{ss}{tt});
        %convert to degree range from 0-360
        neg_angles_logical{ss}{tt} = deg_angles{ss}{tt} < 0;
        %add 360 to negative angles
        deg_angles{ss}{tt}(neg_angles_logical{ss}{tt}) = 360 + deg_angles{ss}{tt}(neg_angles_logical{ss}{tt});
        %convert angles to bin position
        bins{ss}{tt} = discretize(deg_angles{ss}{tt},0:360/100:360);
    end
end

%% Create sort list based on tuning and centroid position (non-matching)

session_nb = 2;

%global indices to use based on trial tuning criteria
tuning_selection = AandB_tuned;

%sort order = 4 - A all; 5 = B - all
trial_sort = 4;

%Get list of A tuned and non-tuned neurons by SI
A_tuned_select = A_STC_tn{session_nb}(:,tuning_selection{session_nb});
B_tuned_select = B_STC_tn{session_nb}(:,tuning_selection{session_nb});
  

%yield index arrangement - sort by A
[centroid_vals,centroid_sortOrder] = sort(bins{session_nb}{trial_sort}(tuning_selection{session_nb}),'ascend');

%remove indices with nans
centroid_sortOrder(isnan(centroid_vals)) = [];

%plot sorted by A centroids
figure;
subplot(1,2,1)
imagesc(A_tuned_select(:,centroid_sortOrder)')
hold on
title(['Session ' num2str(session_nb)])
caxis([0 1])
colormap('jet')

subplot(1,2,2)
imagesc(B_tuned_select(:,centroid_sortOrder)')
hold on
caxis([0 1])
colormap('jet')

%% Generate match list combined with tuning criteria
%global indices to use based on trial tuning criteria
tuning_selection = AorB_tuned;

%tuned to A or B on either sessions
select_match_idx{1} = find(tuning_selection{1} ==1);
select_match_idx{2} = find(tuning_selection{2} ==1);

%intersect with
[tuned_match_idx{1},match_idx{1},~] = intersect(matching_list(:,1),select_match_idx{1},'stable');
[tuned_match_idx{2},match_idx{2},~] = intersect(matching_list(:,2),select_match_idx{2},'stable');

%create not logical for nan exclusion from copied matrix assignement below
include_log{1} = false(1,size(matching_list,1));
include_log{1}(match_idx{1}) = 1; 
%session 2 
include_log{2} = false(1,size(matching_list,1));
include_log{2}(match_idx{2}) = 1;

%make copy
tuned_matching_ROI_list = matching_list;
%nan first session that are not tuned and last session that are not
%tuned
tuned_matching_ROI_list(~include_log{1},1) = nan;
tuned_matching_ROI_list(~include_log{2},2) = nan;

%which neurons to remove based on tuning criterion
keep_ROI = sum(isnan(tuned_matching_ROI_list),2) == 0;

%retain only tuned and matched ROIs
tuned_matching_ROI_list(~keep_ROI,:) = [];

%% Plot matching list

%extract non-sorted maps for each session and trial type
%session 1
session_matched_tuned_STC_tn_maps{1,1} = A_STC_tn{1}(:,tuned_matching_ROI_list(:,1))';
session_matched_tuned_STC_tn_maps{1,2} = B_STC_tn{1}(:,tuned_matching_ROI_list(:,1))';
%session 2
session_matched_tuned_STC_tn_maps{2,1} = A_STC_tn{2}(:,tuned_matching_ROI_list(:,2))';
session_matched_tuned_STC_tn_maps{2,2} = B_STC_tn{2}(:,tuned_matching_ROI_list(:,2))';

%sort by centroid
tuned_matching_ROI_list(:,1)

%get centroids of tuned/matching neurons
%A
centroids_matching{1}{4} = bins{1{4}  
%B
centroids_matching{1}{5} 

%yield index arrangement - sort by A
[centroid_vals,centroid_sortOrder] = sort(bins{session_nb}{trial_sort}(tuning_selection{session_nb}),'ascend');

%remove indices with nans
centroid_sortOrder(isnan(centroid_vals)) = [];


%combined 2x2
combined_maps_2x2 = cell2mat(session_matched_tuned_STC_tn_maps);

%single matrix (ses 1 A, B, ses 2 A,B) - in 1 row
combined_maps_row = cell2mat(reshape(session_matched_tuned_STC_tn_maps,1,4));

%sort by A trials on session 1
[~,matched_maxBin_1] = max(session_matched_tuned_STC_tn_maps{1,1}(:,1:100)', [], 1);
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,matched_sortOrder_1] = sort(matched_maxBin_1,'ascend');

%sort by A trials on session 2
[~,matched_maxBin_2] = max(session_matched_tuned_STC_tn_maps{2,1}(:,1:100)', [], 1);
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,matched_sortOrder_2] = sort(matched_maxBin_2,'ascend');

%session 1 above and session 2 below
figure;
subplot(2,2,1)
imagesc(session_matched_tuned_STC_tn_maps{1,1}(matched_sortOrder_1,:))
colormap('jet')
caxis([0 1]);
subplot(2,2,2)
imagesc(session_matched_tuned_STC_tn_maps{1,2}(matched_sortOrder_1,:))
colormap('jet')
caxis([0 1]);
subplot(2,2,3)
imagesc(session_matched_tuned_STC_tn_maps{2,1}(matched_sortOrder_1,:))
colormap('jet')
caxis([0 1]);
subplot(2,2,4)
imagesc(session_matched_tuned_STC_tn_maps{2,2}(matched_sortOrder_1,:))
colormap('jet')
caxis([0 1]);




%% A vs. B on early vs late training (A or B tuned)

%ALL NEURONS in each session that meet criteria - tuned to either A or B
%early day
for ss =1:size(animal_data,2)
    dF_maps_all_AB_early_late{ss} = [A_STC{ss}(:,AorB_tuned{ss})', B_STC{ss}(:,AorB_tuned{ss})'];
end

%sort each session by A map
for ss =1:size(animal_data,2)
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
imagesc(dF_maps_all_AB_early_late{2}(sortOrder_all_AB{2},:))
hold on
title('Random AB')
colormap('jet')
caxis([0 1])
%A/B vertical separator line
plot([100 100],[1,size(dF_maps_all_AB_early_late{2},1)], 'k','LineWidth', 1.5);


hold off

%% Make matching ROI list with tuning criteria for both sessions

%tuned to A or B on either sessions
AorB_idx{1} = find(AorB_tuned{1} ==1);
AorB_idx{2} = find(AorB_tuned{2} ==1);

%intersect with
[tuned_match_idx{1},match_idx{1},~] = intersect(matching_list(:,1),AorB_idx{1},'stable');
[tuned_match_idx{2},match_idx{2},~] = intersect(matching_list(:,2),AorB_idx{2},'stable');

%create not logical for nan exclusion from copied matrix assignement below
include_log{1} = false(1,size(matching_list,1));
include_log{1}(match_idx{1}) = 1; 
%session 2 
include_log{2} = false(1,size(matching_list,1));
include_log{2}(match_idx{2}) = 1;

%make copy
tuned_matching_ROI_list = matching_list;
%nan first session that are not tuned and last session that are not
%tuned
tuned_matching_ROI_list(~include_log{1},1) = nan;
tuned_matching_ROI_list(~include_log{2},2) = nan;

%which neurons to remove based on tuning criterion
keep_ROI = sum(isnan(tuned_matching_ROI_list),2) == 0;

%retain only tuned and matched ROIs
tuned_matching_ROI_list(~keep_ROI,:) = [];

%% Generate maps based on tuned matching ROI list
%row - session
%column - trial type

session_matched_tuned_dF_maps{1,1} = A_STC{1}(:,tuned_matching_ROI_list(:,1))';
session_matched_tuned_dF_maps{1,2} = B_STC{1}(:,tuned_matching_ROI_list(:,1))';
session_matched_tuned_dF_maps{2,1} = A_STC{2}(:,tuned_matching_ROI_list(:,2))';
session_matched_tuned_dF_maps{2,2} = B_STC{2}(:,tuned_matching_ROI_list(:,2))';

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

