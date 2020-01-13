function [matching_tun_curves] = export_day_matched_STCs(session_vars,session_vars_append,registered)
%export matching STCs for A and B neurons

%go down pipeline and remove ROIs with less than 5 run events

%% Use matching ROIs list 

match_list = registered.multi.assigned_filtered;

min_run_events = 5;

%% Only take neurons with minimum 5 events during run (global)

%GLOBAL minimum of 5 events
%for each session, set neurons to NaNs which have less than 5 Run events
for ss=1:size(match_list,2)
    neuron_idx_less_than_5{ss} = find(session_vars_append{ss}.Events.Run.properties.nb_events < min_run_events);
end

%% Neurons with minimum 5 events during run on A and B trials
%EACH trial type minimum of 5 events
for ss=1:size(match_list,2)
     A_less_than_5_events = find(session_vars{ss}.Events_split{4}.Run.properties.nb_events < min_run_events);
     B_less_than_5_events = find(session_vars{ss}.Events_split{5}.Run.properties.nb_events < min_run_events);
     
     %intersect the 2 sets of ROIs
     neuron_idx_less_than_5_AandB{ss} = intersect(A_less_than_5_events, B_less_than_5_events); 
end

%% Filter out neurons from each column that have less than 5 run events

[match_list_run_event_filtered_5_event_global] = filter_matching_list(match_list,neuron_idx_less_than_5);

%% Filter out neurons from each column that have less than 5 run events on A and B trials

[match_list_run_event_filtered_5_event_both_A_B] = filter_matching_list(match_list,neuron_idx_less_than_5_AandB);


%% Get run event map here and occupancy for A and B trials

%for each session extract these maps
for ss=1:size(session_vars,2)
    %A trials
    %rate map non-smoothed
    event_map{ss}.A = session_vars{ss}.Place_cell{4}.Spatial_Info.rate_map{8}';
    %occupancy in seconds
    occupancy_s{ss}.A = session_vars{ss}.Place_cell{4}.Spatial_Info.occupancy_map{8};
    
    %B trial
    %rate map non-smoothed
    event_map{ss}.B = session_vars{ss}.Place_cell{5}.Spatial_Info.rate_map{8}';
    %occupancy in seconds
    occupancy_s{ss}.B = session_vars{ss}.Place_cell{5}.Spatial_Info.occupancy_map{8};
end


%for each session, get smoothed maps (smoothed rate map/smoothed occupancy)
for ss=1:size(session_vars,2)
    %return rate map that has both events and occupancy time smoothed for A
    %trials
    [rate_map_ev_sm_oc_sm{ss}.A] = return_smoothed_rate_map(event_map{ss}.A,occupancy_s{ss}.A);
    %same for B trials
    [rate_map_ev_sm_oc_sm{ss}.B] = return_smoothed_rate_map(event_map{ss}.B,occupancy_s{ss}.B);
end

%% Extract existing used STC (normalized and non-normalized)

%select which STCs to use
%get for each sessions
for ss=1:size(session_vars,2)
    STC_norm{ss}.A = session_vars{ss}.Place_cell{4}.Spatial_tuning_curve;
    STC_norm{ss}.B = session_vars{ss}.Place_cell{5}.Spatial_tuning_curve;
end

%non-normalized
for ss=1:size(session_vars,2)
    STC_nonNorm{ss}.A = session_vars{ss}.Place_cell{4}.Spatial_tuning_curve_no_norm;
    STC_nonNorm{ss}.B = session_vars{ss}.Place_cell{5}.Spatial_tuning_curve_no_norm;
end

%% Extract matched normalized STCs - at least 5 global run events, normalized STC
%each column is a session
%each row is a different trial type

%STC input here: bin x ROI index
[session_STC_norm_min_5_events_global] = multi_session_STC_extract(match_list_run_event_filtered_5_event_global,STC_norm);

%% Extract matched non-normalized STCs - at least 5 global events

[session_STC_nonNorm_min_5_events_global] = multi_session_STC_extract(match_list_run_event_filtered_5_event_global,STC_nonNorm);

%% Extract matched smoothed rate maps - at least 5 global events (smoothed events over smoothed occupancy)

[session_sm_rate_map_min_5_events_global] = multi_session_STC_extract(match_list_run_event_filtered_5_event_global,rate_map_ev_sm_oc_sm);

%% Extract matched normalized STCs - at least 5 run event on A and B trials

[session_STC_norm_min_5_events_both_A_B] = multi_session_STC_extract(match_list_run_event_filtered_5_event_both_A_B,STC_norm);

%% Extract matched non-normalized STCS - at least 5 run events on A and B trials

[session_STC_nonNorm_min_5_events_both_A_B] = multi_session_STC_extract(match_list_run_event_filtered_5_event_both_A_B,STC_nonNorm);

%% Extract matched smoothed rate maps - at least 5 run events on A and B trials (smoothed events over smoothed occupancy)

[session_sm_rate_map_min_5_events_both_A_B] = multi_session_STC_extract(match_list_run_event_filtered_5_event_both_A_B,rate_map_ev_sm_oc_sm);

%% Arrange as a struct for export purposes

%minimum 5 run events on any lap
matching_tun_curves.min_5_global_events.STC_nonNorm = session_STC_nonNorm_min_5_events_global;
matching_tun_curves.min_5_global_events.STC_norm = session_STC_norm_min_5_events_global;
matching_tun_curves.min_5_global_events.sm_rate_map = session_sm_rate_map_min_5_events_global;

%minimum 5 run events on A and B laps
matching_tun_curves.min_5_events_both_A_B.STC_nonNorm = session_STC_nonNorm_min_5_events_both_A_B;
matching_tun_curves.min_5_events_both_A_B.STC_norm = session_STC_norm_min_5_events_both_A_B;
matching_tun_curves.min_5_events_both_A_B.sm_rate_map = session_sm_rate_map_min_5_events_both_A_B;




end

