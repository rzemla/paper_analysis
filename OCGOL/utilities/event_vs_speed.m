function [outputArg1,outputArg2] = event_vs_speed(session_vars, task_selective_ROIs,ROI_idx_tuning_class,...
                                    select_fields,max_transient_peak,mean_bin_speed,lap_bin_split,options)

%% Define inputs variables

selectTrial = options.selectTrial;
sessionSelect = options.sessionSelect;

%% Define/load variables for each session

%% 

%run ones for each set of trials
split_run_ones = session_vars{1}.Behavior_split_lap.run_ones;

%all run events by lap
run_onset_matrix_each_lap = session_vars{1}.Events_split_lap.Run.run_onset_binary;

%order of correct and incorrect A/B trials
trialOrder = session_vars{1}.Behavior.performance.trialOrder;

%get lap idxs associated with the correct A and B trials
corr_laps_idx.A = find(trialOrder == 2);
corr_laps_idx.B = find(trialOrder ==3);

for ll = 1:size(corr_laps_idx.A,1)
    bin_aligned_onsets.A{ll} = run_onset_matrix_each_lap{corr_laps_idx.A(ll)}(logical(split_run_ones{corr_laps_idx.A(ll)}),:);
end

%number of ROIs in session
nbROI = size(bin_aligned_onsets.A{1, 1},2);

%overlay bin map on top of onsets for each ROI
nb_laps.A = size(corr_laps_idx.A,1);


%repmat each bin assign to number of ROIs
for ll=1:size(lap_bin_split.A,2)
    lap_bin_split_ROI_ex.A{ll} = repmat(lap_bin_split.A{ll},1,nbROI);
end

%bin assign for each ROI
for ll=1:size(lap_bin_split.A,2)
    event_run_bins_mat.A{ll} = lap_bin_split_ROI_ex.A{ll}.*bin_aligned_onsets.A{ll};
end

%for each ROI, get bin activations on each lap
for rr=1:nbROI
    %for each laps
    for ll=1:size(event_run_bins_mat.A,2)
        event_run_bin_extract{rr}{ll} = event_run_bins_mat.A{ll}((event_run_bins_mat.A{ll}(:,rr) ~=0),rr); 
    end
    
end

%place field edges for correct A laps
place_field_bin_edges.A = session_vars{1}.Place_cell{1}.placeField.edge;

%get edges associated with the max transient peak
for rr=1:nbROI
    if~isnan(max_transient_peak{1}{1}(rr))
    max_place_field_bin_edges.A{rr} = place_field_bin_edges.A{rr}(max_transient_peak{1}{1}(rr),:);
    else
        max_place_field_bin_edges.A{rr} = nan;
    end
end

for rr=1:nbROI
    %if place field exists
    if ~isnan(max_place_field_bin_edges.A{rr})
        
    if (max_place_field_bin_edges.A{rr}(1) < max_place_field_bin_edges.A{rr}(2))
        for 
    else %if split field along track end
        
    end
    end
end

%extract the events associated with the correct place field for each ROI


%get corresponding bin for each event - extract from speed script` script


%mean bin speed - mean speed in each bin on each lap





end

