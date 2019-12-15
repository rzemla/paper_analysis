function [outputArg1,outputArg2] = speed_comparison(task_remapping_ROIs, lap_bin_split, session_vars, max_transient_peak,options)

%rate remapping neurons
rate_ROI_idx = task_remapping_ROIs.rate;

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
        
        %bin position
        %get entire interval first
        run_onsets_in_run_epoch{ss}.B{rr} = Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr);
        %extract only indices corresponding to run interval
        run_onsets_in_run_epoch{ss}.B{rr} = run_onsets_in_run_epoch{ss}.B{rr}(logical(session_vars{ss}.Behavior_split{selectTrial(2)}.run_ones));
        
        
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

%% Extract speeds based on discovered events within place field above
for ss=sessionSelect
    %for all ROIs
    for rr=1:size(max_pf_edges.A{ss},2)
        %combined as neighboring cells
        if ~isnan(select_field_idx.A{ss}{rr})
            field_speed_events{ss}{rr,1} = event_speed{ss}.A{rr}(select_field_idx.A{ss}{rr});
        else
            field_speed_events{ss}{rr,1} = nan;
        end
        
        if ~isnan(select_field_idx.B{ss}{rr})
            field_speed_events{ss}{rr,2} = event_speed{ss}.B{rr}(select_field_idx.B{ss}{rr});
        else
            field_speed_events{ss}{rr,2} = nan;
        end
        
        
    end
end

%% Run wilcoxon for all ROIs - check if any sig
for ss=sessionSelect
    %for all ROIs
    for rr=1:size(max_pf_edges.A{ss},2)
        if ~(sum(isnan(field_speed_events{ss}{rr,1})) | sum(isnan(field_speed_events{ss}{rr,2})))
        [speed_pval(rr),speed_true(rr),~] = ranksum(field_speed_events{ss}{rr,1},field_speed_events{ss}{rr,2});
        else
            speed_pval(rr) = nan;
            speed_true(rr) = nan;
        end
    end
end

%% Parse each class based on remapping category - finish here
%extract speed true for remapping states
speed_true(rate_ROI_idx)
field_speed_events{ss}(rate_ROI_idx,:)


%%
end

