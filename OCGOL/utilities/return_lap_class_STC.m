function [STC_nonNorm,STC_norm,rate_map_sm] = return_lap_class_STC(lap_idx,run_events_by_lap,run_bin_by_lap)
%%
%INPUTS
%1 - indices of laps to generate STC out of
%2 - matrix of binary onsets for each ROI on each bin by lap during run
%epoch
%3 - run bins corresponding to the onsets for each 

%OUTPUTS
%1 - nonNorm STC as used in the paper (Gauss sigma=3 smoothed
%events/non-normalized occupancy)
%2 - same STC above, but normalized
%3 - events smoothed (Gauss sigma=3)/ occupancy smoothed (Gauss sigma =3)

%data must be input spatially binned into 100 bins during run epoch

%% Imaging period for 30 Hz (to generate occupancy time)
dt = 0.033422246999976;

%% Split run event and run bin cells into correct and incorrect laps
run_events = run_events_by_lap(lap_idx);
run_bins = run_bin_by_lap(lap_idx);

%% Merge correct run events and correct run bins
events_merged = cell2mat(run_events');
run_bins_merged = cell2mat(run_bins');

%% Sum events for each spatial bin

%for each bin
for bb=1:100
    event_count(:,bb) = sum(events_merged(find(run_bins_merged == bb),:),1);
end

%% Get frame occupancy for each bin and convert to seconds

for bb=1:100
    %occupancy in frames
    occup_fr(bb) = size(find(run_bins_merged == bb),1);
    %occupancy in seconds
    occup_s(bb) = occup_fr(bb)*dt;
      
end

%% Generate STCs for each subgroup of trials

%generate 3 maps: STC (non-norm), STC (norm), rate_map_sm used in paper
[STC_nonNorm,STC_norm,rate_map_sm] = return_all_STCs(event_count,occup_s);



end

