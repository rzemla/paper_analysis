function [outputArg1,outputArg2] = corr_incorr_correlation(session_vars)

%% Input variables

%order of correct and incorrect trials
trialOrder = session_vars{1}.Behavior.performance.trialOrder;

%run and no-run epochs included
position_by_lap = session_vars{1}.Behavior_split_lap.position;

%run ones by lap
run_ones_by_lap = session_vars{1}.Behavior_split_lap.run_ones;

%run ones
run_ones_all_laps = session_vars{1}.Behavior.run_ones;

%significant run events (run and no run epoch included)
lap_run_events = session_vars{1}.Events_split_lap.Run.run_onset_binary;

%normalized position 
norm_position = session_vars{1}.Behavior.resampled.normalizedposition;

%lap number by frame 
lapNb = session_vars{1}.Behavior.resampled.lapNb;


%% Get normalized position and lap assignment during run epochs

run_norm_position = norm_position(logical(run_ones_all_laps));

run_lapNb = lapNb(logical(run_ones_all_laps));

%% Bin normalized run position into 100 bins (across all laps
%assign each position to a bin from 1-100

[N,edges,run_bin] = histcounts(run_norm_position,100);


%% Get run and all normalized position by lap

%all position
%for each lap
for ll=1:max(lapNb)
    norm_position_by_lap{ll} = norm_position(find(lapNb == ll));
end

%run position
for ll=1:max(run_lapNb)
    run_norm_position_by_lap{ll} = run_norm_position(find(run_lapNb == ll));
end

%split run bin position by lap
for ll=1:max(run_lapNb)
    run_bin_by_lap{ll} = run_bin(find(run_lapNb == ll));
end

%extract only the run epochs for the run events
for ll=1:max(run_lapNb)
    run_events_by_lap{ll} = lap_run_events{ll}(logical(run_ones_by_lap{ll}),:);
end

%% Get indices into correct and incorrect laps

A.corr_lap_idx = find(trialOrder == 2);
A.incorr_lap_idx = find(trialOrder == 20);

B.corr_lap_idx = find(trialOrder == 3);
B.incorr_lap_idx = find(trialOrder == 30);

%% Split run event and run bin cells into correct and incorrect laps

A.corr_run_events = run_events_by_lap(A.corr_lap_idx);
A.corr_run_bins = run_bin_by_lap(A.corr_lap_idx);

A.incorr_run_events = run_events_by_lap(A.incorr_lap_idx);
A.incorr_run_bins = run_bin_by_lap(A.incorr_lap_idx);

%% Merge correct run events and correct run bins

%correct
A.corr_events_merged = cell2mat(A.corr_run_events');
A.corr_run_bins_merged = cell2mat(A.corr_run_bins');

%incorrect
A.incorr_events_merged = cell2mat(A.incorr_run_events');
A.incorr_run_bins_merged = cell2mat(A.incorr_run_bins');

%% Sum events for each spatial bin

%for each bin
for bb=1:100
    %correct A
    A.corr.event_count(:,bb) = sum(A.corr_events_merged(find(A.corr_run_bins_merged == bb),:),1);

    %incorrect A
    A.incorr.event_count(:,bb) = sum(A.incorr_events_merged(find(A.incorr_run_bins_merged == bb),:),1);
    
end

%% Get frame occupancy for each bin and convert to seconds

for bb=1:100
    %A correct
    A.corr.occup_fr(bb) = size(find(A.corr_run_bins_merged == bb),1);
    A.corr.occup_s(bb) = A.corr.occup_fr(bb)*session_vars{1}.Imaging.dt;
    
    %A incorrect
    A.incorr.occup_fr(bb) = size(find(A.incorr_run_bins_merged == bb),1);
    A.incorr.occup_s(bb) = A.incorr.occup_fr(bb)*session_vars{1}.Imaging.dt;    
end

%% RESUME HERE

return_smoothed_rate_map

%% Get position during run epochs only for position binning

norm_position_by_lap{ll}(logical(run_ones_by_lap{ll}))


%% Construct non-normalized STCs for correct and incorrect trials for A and B

size(norm_position,1)

cell2mat(lap_run_events');

end

