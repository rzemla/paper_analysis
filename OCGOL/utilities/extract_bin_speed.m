function [mean_bin_speed,corr_lap_bin_split] = extract_bin_speed(corr_lap_idx,session_vars,speed_each_lap,run_epoch_each_lap,trialType)

% correct A laps (100 bins)
corr_lap_bins = session_vars{1}.Place_cell{trialType}.Bin{8};

%speed for each frame on corr A laps
speed_corr = speed_each_lap(corr_lap_idx);

%break run bins into cells (by each corr A lap)
%get run epoch on/off first
run_epoch_corr = run_epoch_each_lap(corr_lap_idx);

%extract 1 only indices (RUN EPOCH ON)
for ll=1:size(run_epoch_corr,2)
    run_only_corr{ll} = ll*ones(size(find(run_epoch_corr{ll}==1),1),1);
end
%convert to vector
run_only_corr_vec = cell2mat(run_only_corr');

%get indices for extracting bins as laps
for ll=1:size(run_only_corr,2)
    lap_bin_idxs(ll,1) = find(run_only_corr_vec == ll,1,'first');
    lap_bin_idxs(ll,2) = find(run_only_corr_vec == ll,1,'last');
end    

%split bins into laps - correct A
for ll=1:size(lap_bin_idxs,1)
    corr_lap_bin_split{ll} = corr_lap_bins(lap_bin_idxs(ll,1):lap_bin_idxs(ll,2));
end

%generate matching speed values for bins for each lap
for ll=1:size(lap_bin_idxs,1)
    speed_match_bin{ll} = speed_corr{ll}(logical(run_epoch_corr{ll}));
end

%get mean speed in each bin on each lap
%for each lap
for ll=1:size(lap_bin_idxs,1)
    %for each bin
    for bb=1:100
        speed_bin_lap{ll}{bb} = speed_match_bin{ll}(find(corr_lap_bin_split{ll} == bb));
    end 
end

%take mean speed in each bin on each lap
for ll=1:size(lap_bin_idxs,1)
    mean_bin_speed(ll,:) = cellfun(@mean,speed_bin_lap{ll});
end


end

