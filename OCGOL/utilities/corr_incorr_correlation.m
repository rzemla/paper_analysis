function [A,B,trial_counts_tbl] = corr_incorr_correlation(session_vars)



%% Input variables

%order of correct and incorrect trials
trialOrder = session_vars.Behavior.performance.trialOrder;

%run and no-run epochs included
position_by_lap = session_vars.Behavior_split_lap.position;

%run ones by lap
run_ones_by_lap = session_vars.Behavior_split_lap.run_ones;

%run ones
run_ones_all_laps = session_vars.Behavior.run_ones;

%significant run events (run and no run epoch included)
lap_run_events = session_vars.Events_split_lap.Run.run_onset_binary;

%normalized position 
norm_position = session_vars.Behavior.resampled.normalizedposition;

%lap number by frame 
lapNb = session_vars.Behavior.resampled.lapNb;

%period of imaging at 30 Hz (sec)
dt = 0.033422246999976;


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

%% Get # of trials/laps for A corr, A incorr, B corr, B incorr

trial_counts(1) = size(A.corr_lap_idx,1);
trial_counts(2) = size(A.incorr_lap_idx,1);
trial_counts(3) = size(B.corr_lap_idx,1);
trial_counts(4) = size(B.incorr_lap_idx,1);

trial_counts_tbl = array2table(trial_counts,'VariableNames',{'A_correct','A_incorrect','B_correct','B_incorrect'});

%% Generate the 3 types of STCs for each of the 4 classes of neurons
%QC checked

%run for each class if at least 1 lap
%A correct STCs
if trial_counts(1) ~= 0
    [A.corr.STC_nonNorm,A.corr.STC_norm,A.corr.rate_map_sm] = ...
        return_lap_class_STC(A.corr_lap_idx,run_events_by_lap,run_bin_by_lap);
else %set maps to empty
    A.corr.STC_nonNorm = [];
    A.corr.STC_norm = [];
    A.corr.rate_map_sm = [];        
end

%A incorrect STCs
if trial_counts(2) ~= 0
    [A.incorr.STC_nonNorm,A.incorr.STC_norm,A.incorr.rate_map_sm] = ...
        return_lap_class_STC(A.incorr_lap_idx,run_events_by_lap,run_bin_by_lap);
else %set maps to empty
    A.incorr.STC_nonNorm = [];
    A.incorr.STC_norm = [];
    A.incorr.rate_map_sm = [];      
end

%B correct STCs
if trial_counts(3) ~= 0
    [B.corr.STC_nonNorm,B.corr.STC_norm,B.corr.rate_map_sm] = ...
        return_lap_class_STC(B.corr_lap_idx,run_events_by_lap,run_bin_by_lap);
else %set maps to empty
    B.corr.STC_nonNorm = [];
    B.corr.STC_norm = [];
    B.corr.rate_map_sm = [];
end

%B incorrect STCs
if trial_counts(4) ~= 0
    [B.incorr.STC_nonNorm,B.incorr.STC_norm,B.incorr.rate_map_sm] = ...
        return_lap_class_STC(B.incorr_lap_idx,run_events_by_lap,run_bin_by_lap);
else %set maps to empty
    B.incorr.STC_nonNorm = [];
    B.incorr.STC_norm = [];
    B.incorr.rate_map_sm = [];    
end


%% Split any input set of laps into 2 random sets and generate STC
%code for shuffling maps - for use later
if 0
    %split STC from any category (insert lap indices) and generate
    
    % # of shuffles
    n_shuffle = 50;
    [STC_nonNorm,rate_map_sm] = generate_STCs_from_shuffled_laps(A.corr_lap_idx,run_events_by_lap,run_bin_by_lap,n_shuffle);
end

%% DEVELOPMENT CODE BELOW %%
% 
% %% Split run event and run bin cells into correct and incorrect laps
% 
% A.corr_run_events = run_events_by_lap(A.corr_lap_idx);
% A.corr_run_bins = run_bin_by_lap(A.corr_lap_idx);
% 
% A.incorr_run_events = run_events_by_lap(A.incorr_lap_idx);
% A.incorr_run_bins = run_bin_by_lap(A.incorr_lap_idx);
% 
% %% Merge correct run events and correct run bins
% 
% %correct
% A.corr_events_merged = cell2mat(A.corr_run_events');
% A.corr_run_bins_merged = cell2mat(A.corr_run_bins');
% 
% %incorrect
% A.incorr_events_merged = cell2mat(A.incorr_run_events');
% A.incorr_run_bins_merged = cell2mat(A.incorr_run_bins');
% 
% %% Sum events for each spatial bin
% 
% %for each bin
% for bb=1:100
%     %correct A
%     A.corr.event_count(:,bb) = sum(A.corr_events_merged(find(A.corr_run_bins_merged == bb),:),1);
% 
%     %incorrect A
%     A.incorr.event_count(:,bb) = sum(A.incorr_events_merged(find(A.incorr_run_bins_merged == bb),:),1);
%     
% end
% 
% %% Get frame occupancy for each bin and convert to seconds
% 
% for bb=1:100
%     %A correct
%     A.corr.occup_fr(bb) = size(find(A.corr_run_bins_merged == bb),1);
%     A.corr.occup_s(bb) = A.corr.occup_fr(bb)*dt;
%     
%     %A incorrect
%     A.incorr.occup_fr(bb) = size(find(A.incorr_run_bins_merged == bb),1);
%     A.incorr.occup_s(bb) = A.incorr.occup_fr(bb)*dt;    
% end
% 
% %% Generate STCs for each subgroup of trials
% 
% %generate literature reported STCs
% %A.corr.rate_map_sm = return_smoothed_rate_map(A.corr.event_count,A.corr.occup_s);
% 
% %generate 3 maps: STC (non-norm), STC (norm), rate_map_sm used in paper
% [A.corr.STC_nonNorm,A.corr.STC_norm,A.corr.rate_map_sm] = return_all_STCs(A.corr.event_count,A.corr.occup_s);
% 
% 
% %% QC check of STCs
% if 0
%     %scroll through each STC 1 by 1
%     figure
%     hold on
%     for rr=1:496
%         subplot(3,1,1)
%         plot(A.corr.rate_map_sm(:,rr))
%         subplot(3,1,2)
%         plot(A.corr.STC_nonNorm(:,rr))
%         subplot(3,1,3)
%         plot(A.corr.STC_norm(:,rr))
%         pause
%         clf
%     end
% end
% 
% %% Get position during run epochs only for position binning
% 
% norm_position_by_lap{ll}(logical(run_ones_by_lap{ll}))
% 
% 
% %% Construct non-normalized STCs for correct and incorrect trials for A and B
% 
% size(norm_position,1)
% 
% cell2mat(lap_run_events');

end

