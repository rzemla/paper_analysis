function [event_run_pf_speed, event_run_bin_extract_pf_sel] = extract_event_speeds(corr_laps_idx, split_run_ones, run_onset_matrix_each_lap,lap_bin_split,...
                                    session_vars,max_transient_peak,mean_bin_speed,trialType)

for ll = 1:size(corr_laps_idx,1)
    bin_aligned_onsets{ll} = run_onset_matrix_each_lap{corr_laps_idx(ll)}(logical(split_run_ones{corr_laps_idx(ll)}),:);
end

%number of ROIs in session
nbROI = size(bin_aligned_onsets{1, 1},2);

%overlay bin map on top of onsets for each ROI
nb_laps = size(corr_laps_idx,1);


%repmat each bin assign to number of ROIs
for ll=1:size(lap_bin_split,2)
    lap_bin_split_ROI_ex{ll} = repmat(lap_bin_split{ll},1,nbROI);
end

%bin assign for each ROI
for ll=1:size(lap_bin_split,2)
    event_run_bins_mat{ll} = lap_bin_split_ROI_ex{ll}.*bin_aligned_onsets{ll};
end

%for each ROI, get bin activations on each lap
for rr=1:nbROI
    %for each laps
    for ll=1:size(event_run_bins_mat,2)
        event_run_bin_extract{rr}{ll} = event_run_bins_mat{ll}((event_run_bins_mat{ll}(:,rr) ~=0),rr); 
    end   
end

%place field edges for correct A laps
place_field_bin_edges = session_vars{1}.Place_cell{trialType}.placeField.edge;

%get edges associated with the max transient peak
for rr=1:nbROI
    if~isnan(max_transient_peak{1}{trialType}(rr))
    max_place_field_bin_edges{rr} = place_field_bin_edges{rr}(max_transient_peak{1}{trialType}(rr),:);
    else
        max_place_field_bin_edges{rr} = nan;
    end
end

%extract the events associated with max transient place field
for rr=1:nbROI
    %if place field exists
    if ~isnan(max_place_field_bin_edges{rr})
        
        if (max_place_field_bin_edges{rr}(1) < max_place_field_bin_edges{rr}(2))
            %define edges of place field
            bin_start = max_place_field_bin_edges{rr}(1);
            bin_end = max_place_field_bin_edges{rr}(2);
            
            for ll=1:nb_laps
                event_idx_temp = find(event_run_bin_extract{rr}{ll} >= bin_start &  event_run_bin_extract{rr}{ll} <= bin_end);
                %extract the event within the field on each lap
                event_run_bin_extract_pf_sel{rr}{ll} = event_run_bin_extract{rr}{ll}(event_idx_temp);
            end
        else %if split field along track end
            %split into 2 finds
            %far edge (near end of track) (from x - 100)
            bin_start = max_place_field_bin_edges{rr}(1);
            %near edge (toward start of track) (from 1 -x)
            bin_end = max_place_field_bin_edges{rr}(2);
            %iterate through each lap
            for ll=1:nb_laps
                end_idx = find(event_run_bin_extract{rr}{ll} >= bin_start &  event_run_bin_extract{rr}{ll} <= 100);
                start_idx = find(event_run_bin_extract{rr}{ll} >= 1 &  event_run_bin_extract{rr}{ll} <= bin_end);
                
                %combine the indices above
                event_idx_temp = sort([start_idx,end_idx]);
                %extract the event within the field on each lap
                event_run_bin_extract_pf_sel{rr}{ll} = event_run_bin_extract{rr}{ll}(event_idx_temp);
            end
        end
    else
        event_run_bin_extract_pf_sel{rr} = [];
    end
end

%extract bin speed associated with each event
for rr=1:nbROI
    if ~isempty(event_run_bin_extract_pf_sel{rr})
        for ll=1:nb_laps
            event_run_pf_speed{rr}{ll} = mean_bin_speed(ll,event_run_bin_extract_pf_sel{rr}{ll});
        end
    else
         event_run_pf_speed{rr}{ll} = [];
    end
end



%extract task-selective speeds
% Asel_speed = event_run_pf_speed(task_selective_ROIs.A.idx);
% 
% %convert each ROI speed to mat
% for rr=1:size(Asel_speed,2)
%     Asel_speed_mat{rr} = cell2mat(Asel_speed{rr});
%     
% end

end

