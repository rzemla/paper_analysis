function [remapping_corr_idx,remap_idx_traces,com_idx_traces,AUC_remappers] = speed_AUC_comparison(task_remapping_ROIs, remapping_corr_idx, lap_bin_split, session_vars, max_transient_peak, STC_export,event_5_min_and_occup_filtered_ROI,ABtuned_all_si_ts,options)

%QC checked

%rate remapping neurons
rate_ROI_idx = task_remapping_ROIs.rate;

%non global ROIs from rate map correlations
non_global_ROIs = remapping_corr_idx.non_global_idx;

%partial_ROI (near/common field defined from threshold and no limit on far
%field)
partial_ROIs = task_remapping_ROIs.partial;

%session indices
sessionSelect = options.sessionSelect;
%which trials to use
selectTrial = options.selectTrial;

%% Define/load variables for each session

%for each session
for ii = sessionSelect
    % behavior and imaging related variables
    Behavior_split_lap{ii} = session_vars{ii}.Behavior_split_lap;
    Events_split_lap{ii} = session_vars{ii}.Events_split_lap;
    Behavior_split{ii} = session_vars{ii}.Behavior_split;
    Event_split{ii} = session_vars{ii}.Events_split;
    Imaging_split{ii} = session_vars{ii}.Imaging_split;
    Place_cell{ii} = session_vars{ii}.Place_cell;
    Behavior_full{ii} = session_vars{ii}.Behavior;
    
    %all within run domain
    position{ii}  = Behavior_split_lap{ii}.Run.position;
    time{ii} = Behavior_split_lap{ii}.Run.time;
    events_full{ii} = Events_split_lap{ii}.Run.run_onset_binary;
    run_intervals{ii} = Behavior_split_lap{ii}.run_ones;
    
    %global trial type order across restricted laps
    trialOrder{ii} = Behavior_full{ii}.performance.trialOrder;
end

%for each session
for ss=sessionSelect
    %for each lap
    for ii=1:size(run_intervals{ss},2)
        events{ss}{ii} = events_full{ss}{ii}(logical(run_intervals{ss}{ii}),:);
    end
end

%% Convert lap bins into continuous vector

A_run_bins = cell2mat(lap_bin_split.A');
B_run_bins = cell2mat(lap_bin_split.B');


%% Get events in run epochs

for ss=sessionSelect
    %for each ROI
    for rr=1:size(events{ss}{1},2) %get the ROI size
        %get events in run epoch only in order to extract associated run
        %bin position
        %get entire interval first
        run_onsets_in_run_epoch{ss}.A{rr} = Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr);
        %extract only indices corresponding to run interval
        run_onsets_in_run_epoch{ss}.A{rr} = run_onsets_in_run_epoch{ss}.A{rr}(logical(session_vars{ss}.Behavior_split{selectTrial(1)}.run_ones));
        
        %event traces
        event_traces{ss}.A{rr} = Event_split{ss}{selectTrial(1)}.Run.properties.trace{rr};  
        
        %event bin positions
        event_run_bins{ss}.A{rr} = A_run_bins(find(run_onsets_in_run_epoch{ss}.A{rr} ==1));
        
        %get run onsets in run only interval
        event_run_idx_whole{ss}.A{rr} = find(run_onsets_in_run_epoch{ss}.A{rr} ==1);
        
        %time of significant run events in A
        event_norm_time{ss}.A{rr} = Imaging_split{ss}{selectTrial(1)}.time_restricted(find(Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr) == 1))/60;
        %normalizesd position of significant run events in A
        event_norm_pos_run{ss}.A{rr} = Behavior_split{ss}{selectTrial(1)}.resampled.position_norm(find(Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr) == 1));
        %lap assignment for A
        event_lap_idx{ss}.A{rr} = Behavior_split{ss}{selectTrial(1)}.resampled.lapNb(logical(Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr)));
        %event speed on each lap
        event_speed{ss}.A{rr} = Behavior_split{ss}{selectTrial(1)}.speed(find(Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr) == 1));
        
        %get AUC events here
        event_AUC{ss}.A{rr} = Event_split{ss}{selectTrial(1)}.Run.properties.AUC{rr};
        
        
        %bin position
        %get entire interval first
        run_onsets_in_run_epoch{ss}.B{rr} = Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr);
        %extract only indices corresponding to run interval
        run_onsets_in_run_epoch{ss}.B{rr} = run_onsets_in_run_epoch{ss}.B{rr}(logical(session_vars{ss}.Behavior_split{selectTrial(2)}.run_ones));
        
        %event traces
        event_traces{ss}.B{rr} = Event_split{ss}{selectTrial(2)}.Run.properties.trace{rr};  
        
        %event bin positions
        event_run_bins{ss}.B{rr} = B_run_bins(find(run_onsets_in_run_epoch{ss}.B{rr} ==1));
        
        %get run onsets in run only interval
        event_run_idx_whole{ss}.B{rr} = find(run_onsets_in_run_epoch{ss}.B{rr} ==1);        
        
        %time of significant run events in B
        event_norm_time{ss}.B{rr} = Imaging_split{ss}{selectTrial(2)}.time_restricted(find(Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr) == 1))/60;
        %normalizesd position of significant run events in B
        event_norm_pos_run{ss}.B{rr} = Behavior_split{ss}{selectTrial(2)}.resampled.position_norm(find(Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr) == 1));
        %lap assignment for B
        event_lap_idx{ss}.B{rr} = Behavior_split{ss}{selectTrial(2)}.resampled.lapNb(logical(Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr)));
        %event speed on each lap
        event_speed{ss}.B{rr} = Behavior_split{ss}{selectTrial(2)}.speed(find(Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr) == 1));        
        
        %get AUC events here
        event_AUC{ss}.B{rr} = Event_split{ss}{selectTrial(2)}.Run.properties.AUC{rr};
        
    end
end

%% Extract max transient rate place fields

%extract max transient fields
for ss=sessionSelect
    %max fields
    for rr=1:size(events{ss}{1},2)
        %for correct A trials, extract max transient place field
        if ~isnan(max_transient_peak{ss}{1}(rr))
            max_pf_edges.A{ss}{rr} = session_vars{ss}.Place_cell{1}.placeField.edge{rr}(max_transient_peak{ss}{1}(rr),:);
        else
            max_pf_edges.A{ss}{rr} = nan;
        end
        
        %for correct B trials, extract max transient place field
        if ~isnan(max_transient_peak{ss}{2}(rr))
            max_pf_edges.B{ss}{rr} = session_vars{ss}.Place_cell{2}.placeField.edge{rr}(max_transient_peak{ss}{2}(rr),:);
        else
            max_pf_edges.B{ss}{rr} = nan;
        end

    end
end

%% Extract indices of events that are within the max transient rate place field

%extract bins that fall within max field (for A laps)
for ss=sessionSelect
    for rr=1:size(max_pf_edges.A{ss},2)
        %check if empty
        if ~isnan(max_pf_edges.A{ss}{rr})
            %if field not crossing the edge of track
            if max_pf_edges.A{ss}{rr}(2) > max_pf_edges.A{ss}{rr}(1)
                
                select_field_idx.A{ss}{rr} = find(event_run_bins{ss}.A{rr} >=max_pf_edges.A{ss}{rr}(1) & event_run_bins{ss}.A{rr} <=max_pf_edges.A{ss}{rr}(2));
                
            else
                disp(rr)
                idx_temp_end = find(event_run_bins{ss}.A{rr} >=max_pf_edges.A{ss}{rr}(1) & event_run_bins{ss}.A{rr} <=100);
                idx_temp_start = find(event_run_bins{ss}.A{rr} >=1  & event_run_bins{ss}.A{rr} <=max_pf_edges.A{ss}{rr}(2));
                
                select_field_idx.A{ss}{rr} = sort([idx_temp_end; idx_temp_start]);
            end
        else
            select_field_idx.A{ss}{rr} = nan;
        end
        
    end
end

%extract bins that fall within max field (for B laps)
for ss=sessionSelect
    for rr=1:size(max_pf_edges.B{ss},2)
        %check if empty
        if ~isnan(max_pf_edges.B{ss}{rr})
            %if field not crossing the edge of track
            if max_pf_edges.B{ss}{rr}(2) > max_pf_edges.B{ss}{rr}(1)
                %take start to end range
                select_field_idx.B{ss}{rr} = find(event_run_bins{ss}.B{rr} >=max_pf_edges.B{ss}{rr}(1) & event_run_bins{ss}.B{rr} <=max_pf_edges.B{ss}{rr}(2));
                
            else
                %split into 2 sets of edges and concatenate
                idx_temp_end = find(event_run_bins{ss}.B{rr} >=max_pf_edges.B{ss}{rr}(1) & event_run_bins{ss}.B{rr} <=100);
                idx_temp_start = find(event_run_bins{ss}.B{rr} >=1  & event_run_bins{ss}.B{rr} <=max_pf_edges.B{ss}{rr}(2));
                
                select_field_idx.B{ss}{rr} = sort([idx_temp_end; idx_temp_start]);
            end
        else
            select_field_idx.B{ss}{rr} = nan;
        end
        
    end
end

%% Extract speeds, traces, and AUC based on discovered events within place field above
for ss=sessionSelect
    %for all ROIs
    for rr=1:size(max_pf_edges.A{ss},2)
        %combined as neighboring cells
        if ~isnan(select_field_idx.A{ss}{rr})
            %speed
            field_speed_events{ss}{rr,1} = event_speed{ss}.A{rr}(select_field_idx.A{ss}{rr});
            %AUC
            field_AUC_events{ss}{rr,1} = event_AUC{ss}.A{rr}(select_field_idx.A{ss}{rr});
            %trace
            field_trace_events{ss}{rr,1} = event_traces{ss}.A{rr}(select_field_idx.A{ss}{rr});
            
        else
            field_speed_events{ss}{rr,1} = nan;
            field_AUC_events{ss}{rr,1} = nan;
            field_trace_events{ss}{rr,1} = nan;
        end
        
        if ~isnan(select_field_idx.B{ss}{rr})
            field_speed_events{ss}{rr,2} = event_speed{ss}.B{rr}(select_field_idx.B{ss}{rr});
            %AUC
            field_AUC_events{ss}{rr,2} = event_AUC{ss}.B{rr}(select_field_idx.B{ss}{rr});
            %trace
            field_trace_events{ss}{rr,2} = event_traces{ss}.B{rr}(select_field_idx.B{ss}{rr});
        else
            field_speed_events{ss}{rr,2} = nan;
            field_AUC_events{ss}{rr,2} = nan;
            field_trace_events{ss}{rr,2} = nan;
        end
        
    end
end


%% Plot sample traces for rate remapping neurons

% figure
% subplot(2,1,1)
% hold on
% ylim([0 3])
% xlim([0 250])
% for ii=1:9
%     plot(field_trace_events{1, 1}{278, 1}{ii},'b')
% end
% 
% subplot(2,1,2)
% hold on
% ylim([0 3])
% xlim([0 250])
% for ii=1:9
%     plot(field_trace_events{1, 1}{278, 2}{ii},'r')
% end


%% Plot scatter of speed vs AUC for A events and B events

figure
for rr=251%rate_ROI_idx
    hold on
    xlim([0 25])
    xlabel('Speed [cm/s]')
    ylim([0 7])
    ylabel('AUC [dF/F*s]')
    %speed -x; AUC - y
    %A trials
    scatter(field_speed_events{1}{rr,1},field_AUC_events{1}{rr,1},20,'filled','MarkerFaceColor','b' )
    %B trials
    scatter(field_speed_events{1}{rr,2},field_AUC_events{1}{rr,2},20,'filled','MarkerFaceColor','r' )
    %pause
    %clf
end


%% Run 2-way anova on AUC values and speed values
%(Y1, Y2, X1, X2)
% Y1: AUC of group 1
% Y2: AUC of group 2
% X1: Speed of group 1
% X2: Speed of group 2
% for rr=rate_ROI_idx
%     [pGroup(rr), pAll(:,rr), table] = jz_anova2(field_AUC_events{1}{rr,1},field_AUC_events{1}{rr,2}, field_speed_events{1}{rr,1},field_speed_events{1}{rr,2} )
% end

%preallocate vectors
pGroup = nan(1,size(events{ss}{1},2));
pAll = nan(3,size(events{ss}{1},2));

for rr=sort(non_global_ROIs)'
    [pGroup(rr), pAll(:,rr), table] = jz_anova2(field_AUC_events{1}{rr,1},field_AUC_events{1}{rr,2}, field_speed_events{1}{rr,1},field_speed_events{1}{rr,2} );
end

%get all signifcant p-values
allSig = pAll < 0.05;
%either state for other terms
otherTermsOr = allSig(2,:) | allSig(3,:);

%find ROI with group effect regardless of other interactions
rate_remap_ROI_group_all = find(pGroup <0.05);

%find ROI with group effect only (cell identity explains AUC)
rate_remap_group_only = setdiff(rate_remap_ROI_group_all,intersect(rate_remap_ROI_group_all,find(otherTermsOr ==1)));

%% Get event onset/offset in frames
%onset and offset frame of event
Event_split{1, 1}{1, 1}.Run.run_onset_offset{1, 1};
%trace of event
Event_split{1, 1}{1, 1}.Run.properties.trace;

event_run_bins;

%% Run wilcoxon for all ROIs - check if any sig
for ss=sessionSelect
    %for all ROIs
    for rr=1:size(max_pf_edges.A{ss},2)
        if ~(sum(isnan(field_speed_events{ss}{rr,1})) | sum(isnan(field_speed_events{ss}{rr,2})))
        [speed_pval(rr),speed_true(rr),~] = ranksum(field_speed_events{ss}{rr,1},field_speed_events{ss}{rr,2});
        else
            speed_pval(rr) = 0;
            speed_true(rr) = 0;
        end
    end
end

%% Plot example STCS (trial normalized)
if 0
    %generate combined STC for selected ROIs
    STC_tn_partial = [STC_export.A_STC_tn{1}(:,partial_ROIs)', STC_export.B_STC_tn{1}(:,partial_ROIs)'];
    
    %center bins - get sort order from here
    [~,com_partial_idx_sort] = sort(bin_center.partial_com(1,:),'ascend');
    
    %plot the sorted partial neurons
    figure
    imagesc(STC_tn_partial(com_partial_idx_sort,:))
    hold on
    caxis([0 1])
    colormap('jet')
end

%% Filter global events based on occupancy metric from remapping_categorize code

global_idx_occup_filt = intersect(remapping_corr_idx.global_idx,event_5_min_and_occup_filtered_ROI);


%% Get number of unclassified neurons

%common tuned neurons
common_idx = setdiff(non_global_ROIs,rate_remap_ROI_group_all);

%rate remapping
rate_remap_ROI_group_all;
%global remapping
global_idx_occup_filt;
%all neurons A&B tuned
ABtuned_all_si_ts;

%define unclassified category
%remove common
common_rem = setdiff(ABtuned_all_si_ts,common_idx);
%remove global
com_glob_rem = setdiff(common_rem,global_idx_occup_filt);
%remove rate
com_glob_rate_rem = setdiff(com_glob_rem,rate_remap_ROI_group_all);
%remove partial
unclass_idx = setdiff(com_glob_rate_rem,partial_ROIs);

%% Parse each class based on remapping category - finish here
%extract speed true for remapping states
speed_true(rate_ROI_idx)
field_speed_events{ss}(rate_ROI_idx,:)

%% Export rate remapping ROIs

%rate remapping neurons - by only category group
remapping_corr_idx.final.rate_remap_grp_only = rate_remap_group_only;
remapping_corr_idx.final.rate_remap_all = rate_remap_ROI_group_all;


%define common in relation remove the rate remapping neurons from
%non_global_ROIs
remapping_corr_idx.final.common = common_idx;

%global idx (occupancy filtered)
remapping_corr_idx.final.global = global_idx_occup_filt;

%partial idx
remapping_corr_idx.final.partial = partial_ROIs;

%unclassified idx
remapping_corr_idx.final.unclass = unclass_idx;

%all A&B tuned by SI or TS criteria
remapping_corr_idx.final.AB_tuned_si_ts = ABtuned_all_si_ts;

%% Extract AUC values for A and B trials for common and remapping neurons

remapping_corr_idx.final.rate_remap_all
remapping_corr_idx.final.common

%extract AUC of each event for common neurons
for rr=1:size(remapping_corr_idx.final.common,1)
    %A AUC events
    common_AUC{rr,1} = field_AUC_events{1}{remapping_corr_idx.final.common(rr),1};
    %B AUC events
    common_AUC{rr,2} = field_AUC_events{1}{remapping_corr_idx.final.common(rr),2};
end

%extract AUC of each event for common neurons
for rr=1:size(remapping_corr_idx.final.rate_remap_all,2)
    %A AUC events
    remap_AUC{rr,1} = field_AUC_events{1}{remapping_corr_idx.final.rate_remap_all(rr),1};
    %B AUC events
    remap_AUC{rr,2} = field_AUC_events{1}{remapping_corr_idx.final.rate_remap_all(rr),2};
end

%use for export
AUC_remappers.common = common_AUC;
AUC_remappers.remap = remap_AUC;

% x = abs(cellfun(@sum,common_AUC(:,1)) - cellfun(@sum,common_AUC(:,2)))
% 
% y = abs(cellfun(@sum,remap_AUC(:,1)) - cellfun(@sum,remap_AUC(:,2)))
% 
% figure
% hold on
% histogram(x,10)
% histogram(y,10)

%% Extract and export traces associated with each remapping neuron
nb_remap_idx = size(rate_remap_ROI_group_all,2);

idx_start = 1;
for rr=rate_remap_ROI_group_all
    [traces_combined_remap{idx_start}.A] = extract_traces(field_trace_events{1}{rr,1});
    [traces_combined_remap{idx_start}.B] = extract_traces(field_trace_events{1}{rr,2});
    idx_start = idx_start + 1;
end

%export this variable
remap_idx_traces = traces_combined_remap;

%% Extract and export traces associated with each common neuron (for histogram comparison against remapping neurons)
nb_common_idx = size(common_idx,1);

idx_start = 1;
for rr=common_idx'
    [traces_combined_com{idx_start}.A] = extract_traces(field_trace_events{1}{rr,1});
    [traces_combined_com{idx_start}.B] = extract_traces(field_trace_events{1}{rr,2});
    idx_start = idx_start + 1;
end

%export this variable
com_idx_traces = traces_combined_com;

%%
end

