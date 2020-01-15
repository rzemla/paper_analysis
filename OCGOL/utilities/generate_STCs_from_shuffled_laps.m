function [STC_nonNorm,rate_map_sm] = generate_STCs_from_shuffled_laps(lap_idx,run_events_by_lap,run_bin_by_lap,n_shuffle)

%number of shuffles
%n_shuffle = 100;

%preallocate cells for storing idxs
split_lap_set = cell(2,n_shuffle);

for ii=1:n_shuffle
    %permute the indices
    idx_permute = randperm(size(lap_idx,1));
    %find the middle index
    mid_idx = ceil(size(lap_idx,1)/2);
    
    %split the numeric indices into two cells
    split_lap_set{1,ii} = lap_idx(idx_permute(1:mid_idx));
    split_lap_set{2,ii} = lap_idx(idx_permute((mid_idx+1):end));
end

%preallocate cells for storing STCs from each split of laps
STC_nonNorm = cell(2,n_shuffle);
rate_map_sm = cell(2,n_shuffle);

for ii=1:n_shuffle
%generate STCs for first half
    [STC_nonNorm{1,ii},~,rate_map_sm{1,ii}] = ...
        return_lap_class_STC(split_lap_set{1,ii},run_events_by_lap,run_bin_by_lap);
%generate STCs for second half
    [STC_nonNorm{2,ii},~,rate_map_sm{2,ii}] = ...
        return_lap_class_STC(split_lap_set{2,ii},run_events_by_lap,run_bin_by_lap);
end

end

