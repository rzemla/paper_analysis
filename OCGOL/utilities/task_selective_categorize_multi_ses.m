function [task_selective_ROIs] = task_selective_categorize_multi_ses(tuned_logical,session_vars, max_transient_peak, options)
%split mutually tuned neurons by remapping category: 
%common (less than certain centroid difference between max
%tuned_log = tunedLogical.ts.AandB_tuned;

%QC checked
%commented out speed

%% Define inputs variables

selectTrial = options.selectTrial;
sessionSelect = options.sessionSelect;

%% Get ROI indices of A selective, B selective, A&B selective neurons

% %choose if SI or TS tuned
% switch options.tuning_criterion
%     case 'si' %spatial information
%         AandB_tuned_idx = find(tuned_logical.si.AandB_tuned == 1);
%         Aonly_tuned_idx = find(tuned_logical.si.onlyA_tuned == 1);
%         Bonly_tuned_idx = find(tuned_logical.si.onlyB_tuned == 1);
%

%for each session
for ss=sessionSelect
    %case 'ts' %tuning specificity
    %both A and B tuned by TS
    AandB_tuned_idx{ss} = find(tuned_logical(ss).ts.AandB_tuned ==1);
    %only A tuned by TS
    Aonly_tuned_idx{ss} = find(tuned_logical(ss).ts.onlyA_tuned == 1);
    %only B tuned by tS
    Bonly_tuned_idx{ss} = find(tuned_logical(ss).ts.onlyB_tuned == 1);
    
    %SI only tuned in either
    Aonly_tuned_idx_si{ss} = find(tuned_logical(ss).si.onlyA_tuned == 1);
    Bonly_tuned_idx_si{ss} = find(tuned_logical(ss).si.onlyB_tuned == 1);
    
    %all A tuned by SI
    A_tuned_si_idx{ss} = find(tuned_logical(ss).si.Atuned == 1);
    %all B tuned by SI
    B_tuned_si_idx{ss} = find(tuned_logical(ss).si.Btuned == 1);
    
    %all A tuned by TS
    A_tuned_ts_idx{ss} = find(tuned_logical(ss).ts.Atuned == 1);
    %all B tuned by TS
    B_tuned_ts_idx{ss} = find(tuned_logical(ss).ts.Btuned == 1);
    
    %only A tuned by TS, but also not SI B tuned
    Aonly_notSIb_idx{ss} = setdiff(Aonly_tuned_idx{ss},B_tuned_si_idx{ss});
    
    %only B tuned by TS, but also not SI A tuned
    Bonly_notSIa_idx{ss} = setdiff(Bonly_tuned_idx{ss},A_tuned_si_idx{ss});
    
    %only A tuned by SI, but also not TS B tuned
    Aonly_notTSb_idx{ss} = setdiff(Aonly_tuned_idx_si{ss},B_tuned_ts_idx{ss});
    
    %only B tuned by SI, but also not TS A tuned
    Bonly_notTSa_idx{ss} = setdiff(Bonly_tuned_idx_si{ss},A_tuned_ts_idx{ss});
    
    %merge A and B chosen by exclusive
    A_only_idx{ss} = union(Aonly_notSIb_idx{ss}, Aonly_notTSb_idx{ss});
    B_only_idx{ss} = union(Bonly_notSIa_idx{ss},Bonly_notTSa_idx{ss});
    %end
end

%% Spatial crtiteria input ROIs - will be used as initial preselected ROIs
%for all downstream processing below
%TS tuned in either A exclusive or B exclusive and also no opposing trial
%typing by SI criterion
switch options.tuning_criterion
    case 'ts'
        input_idx_Aonly = Aonly_notSIb_idx;
        input_idx_Bonly = Bonly_notSIa_idx;
        
    case 'both'
        %SI only and TS only tuned for either A or B trials
        input_idx_Aonly = A_only_idx;
        input_idx_Bonly = B_only_idx;
end

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

%% Get lap indices for each lap in all B or only correct B trials

for ss=sessionSelect
    %get unique lap indices
    lapA_idxs{ss} = unique(Behavior_split{ss}{selectTrial(1)}.resampled.lapNb);
    lapB_idxs{ss} = unique(Behavior_split{ss}{selectTrial(2)}.resampled.lapNb);
    
    %get lap start and end indices for all A or B trials
    %all A
    for ll=1:size(lapA_idxs{ss},1)
        lap_idxs{ss}.A(ll,1) = find(Behavior_split{ss}{selectTrial(1)}.resampled.lapNb == lapA_idxs{ss}(ll),1,'first');
        lap_idxs{ss}.A(ll,2) = find(Behavior_split{ss}{selectTrial(1)}.resampled.lapNb == lapA_idxs{ss}(ll),1,'last');
    end
    
    %all B
    for ll=1:size(lapB_idxs{ss},1)
        lap_idxs{ss}.B(ll,1) = find(Behavior_split{ss}{selectTrial(2)}.resampled.lapNb == lapB_idxs{ss}(ll),1,'first');
        lap_idxs{ss}.B(ll,2) = find(Behavior_split{ss}{selectTrial(2)}.resampled.lapNb == lapB_idxs{ss}(ll),1,'last');
    end
    
end

%% Event onsets in run interval
%use trial definition set by trialSelect
%get speed data for each event here

for ss=sessionSelect
    %for each ROI
    for rr=1:size(events{ss}{1},2) %get the ROI size
        %time of significant run events in A
        event_norm_time{ss}.A{rr} = Imaging_split{ss}{selectTrial(1)}.time_restricted(find(Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr) == 1))/60;
        %normalizesd position of significant run events in A
        event_norm_pos_run{ss}.A{rr} = Behavior_split{ss}{selectTrial(1)}.resampled.position_norm(find(Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr) == 1));
        %lap assignment for A
        event_lap_idx{ss}.A{rr} = Behavior_split{ss}{selectTrial(1)}.resampled.lapNb(logical(Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr)));
        %event speed on each lap
        %event_speed{ss}.A{rr} = Behavior_split{ss}{selectTrial(1)}.speed(find(Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr) == 1));
        
        %time of significant run events in B
        event_norm_time{ss}.B{rr} = Imaging_split{ss}{selectTrial(2)}.time_restricted(find(Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr) == 1))/60;
        %normalizesd position of significant run events in B
        event_norm_pos_run{ss}.B{rr} = Behavior_split{ss}{selectTrial(2)}.resampled.position_norm(find(Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr) == 1));
        %lap assignment for B
        event_lap_idx{ss}.B{rr} = Behavior_split{ss}{selectTrial(2)}.resampled.lapNb(logical(Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr)));
        %event speed on each lap
        %event_speed{ss}.B{rr} = Behavior_split{ss}{selectTrial(2)}.speed(find(Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr) == 1));        
    end
end

%% Plot the run epochs, corresponding position bin edges of place field

%notes on data input
%overlay area of max place field bin
%spatial bin assignment for each run-epoch frame (100 bins)
%A trials
%session_vars{ss}.Place_cell{1}.Bin{8};
%B trials
%session_vars{ss}.Place_cell{2}.Bin{8};

%get the run epoch binaries for each set of trials - A session - 5339 run
%frames
%plot the run epochs as patches
%session_vars{ss}.Behavior_split{1}.run_ones;

%get the run epoch binaries for each set of trials - B session - 4547 run
%frames
%session_vars{ss}.Behavior_split{2}.run_ones;

%run events - binary onset across entire run/no run interval
%Event_split{ss}{1}.Run.run_onset_ones;

%get bins
for ss=sessionSelect
    %get edges for corresponding bins!! - find place in spatial info where
    %correct A trials
    run_position_norm{ss}{selectTrial(1)} = Behavior_split{ss}{selectTrial(1)}.resampled.run_position_norm;
    %correct B trials
    run_position_norm{ss}{selectTrial(2)} = Behavior_split{ss}{selectTrial(2)}.resampled.run_position_norm;
    
    %Bin running position in 100 bins and get edges for each set of laps:
    %for each number of bins, bin the normalized position during run epochs
    
    %for correct A trials
    [count_bin{ss}{selectTrial(1)},edges{ss}{selectTrial(1)},bin{ss}{selectTrial(1)}] = histcounts(run_position_norm{ss}{selectTrial(1)}, 100);
    %for correct B trials
    [count_bin{ss}{selectTrial(2)},edges{ss}{selectTrial(2)},bin{ss}{selectTrial(2)}] = histcounts(run_position_norm{ss}{selectTrial(2)}, 100);
end

%% Remove ROIs idx's without a id'd place field
for ss=sessionSelect
    %edges of all identified
    %correct A
    placeFieldEdges{ss}{selectTrial(1)} = Place_cell{ss}{selectTrial(1)}.placeField.edge(input_idx_Aonly{ss});
    %correct B
    placeFieldEdges{ss}{selectTrial(2)} = Place_cell{ss}{selectTrial(2)}.placeField.edge(input_idx_Bonly{ss});

    %in A trials
    idx_wo_placeFields{ss}{selectTrial(1)} = input_idx_Aonly{ss}(find(cellfun(@isempty,placeFieldEdges{ss}{selectTrial(1)}) ==1));
    %in B trials
    idx_wo_placeFields{ss}{selectTrial(2)} = input_idx_Bonly{ss}(find(cellfun(@isempty,placeFieldEdges{ss}{selectTrial(2)}) ==1));
   
    %remove ROIs for A and B that do not have place fields
    %for A trials
    if ~isempty(idx_wo_placeFields{ss}{selectTrial(1)})
        %copy
        Aonly_field_filtered{ss} = setdiff(input_idx_Aonly{ss},idx_wo_placeFields{ss}{selectTrial(1)});
        
    else %keep the previous ROIs (copy only)
        Aonly_field_filtered{ss} = input_idx_Aonly{ss};
    end
    %for B trials
    if ~isempty(idx_wo_placeFields{ss}{selectTrial(2)})
        %copy
        Bonly_field_filtered{ss} = setdiff(input_idx_Bonly{ss},idx_wo_placeFields{ss}{selectTrial(2)});
        
    else %keep the previous ROIs (copy only)
        Bonly_field_filtered{ss} = input_idx_Bonly{ss};
    end
    
end

%% Determine the equivalent normalized position range of the max place field

for ss=sessionSelect
    %extract the max place field index for A and B trial using filtered A tuned
    %and B tuned ROIs
    max_field_idx{ss}{selectTrial(1)} = max_transient_peak{ss}{selectTrial(1)}(Aonly_field_filtered{ss});
    max_field_idx{ss}{selectTrial(2)} = max_transient_peak{ss}{selectTrial(2)}(Bonly_field_filtered{ss});
    
    %get the edges of the max transient place field for each set of idxs
    %edges of all identified
    %correct A
    placeField_filtered{ss}{selectTrial(1)} = Place_cell{ss}{selectTrial(1)}.placeField.edge(Aonly_field_filtered{ss});
    %correct B
    placeField_filtered{ss}{selectTrial(2)} = Place_cell{ss}{selectTrial(2)}.placeField.edge(Bonly_field_filtered{ss});
    
    %check which field has more than 1 field and select edges of the one with
    %higher transient rate
    
    %for each trial (A and B )
    for tt=selectTrial
        for rr=1:size(placeField_filtered{ss}{tt},2)
            if size(placeField_filtered{ss}{tt}{rr},1) > 1
                placeField_filtered_max{ss}{tt}{rr} = placeField_filtered{ss}{tt}{rr}(max_field_idx{ss}{tt}(rr),:);
            else
                placeField_filtered_max{ss}{tt}{rr} = placeField_filtered{ss}{tt}{rr};
            end
        end
    end
   
    %convert edges from relevant place field to normalized postion edges
    %trial for edges can be either - same thing
    for tt=selectTrial
        for rr=1:size(placeField_filtered_max{ss}{tt},2)
            %start position of PF
            placeField_filtered_max_posnorm{ss}{tt}{rr}(1) = edges{ss}{selectTrial(1)}(placeField_filtered_max{ss}{tt}{rr}(1));
            %end position of PF
            placeField_filtered_max_posnorm{ss}{tt}{rr}(2) = edges{ss}{selectTrial(1)}(placeField_filtered_max{ss}{tt}{rr}(2)+1)-0.01;
        end
    end
end

%% Get normalized position distance converstion factor
for ss=sessionSelect
    %get median lap length based on the registered length of each lap
    median_track_len(ss) = median(Behavior_full{ss}.position_lap(:,2));
    
    %conversion factor (norm_pos/cm length) - 100 bins
    norm_conv_factor(ss) = median_track_len(ss)/100;
    %median_track_len/1;
end

%% Filter out neurons that do not have at least 5 sig events in max place field in at least 5 distinct laps

%calcium event position and absolute restrict time
%all neuron idx space
% event_norm_pos_run.A;
% event_norm_time.A;  
% event_lap_idx.A;
% 
% Aonly_field_filtered;
% Bonly_field_filtered;

%mopdify this for discontinuous place fields

%for each session
for ss=sessionSelect
    %find events occuring within max place field for each ROI
    for tt=selectTrial
        for rr=1:size(placeField_filtered_max_posnorm{ss}{tt},2)
            %only 1 max place field here so not need to check 
            %get idxs of events with max place field            
            if tt == selectTrial(1) %correct A trials or all A
                %check if first edge is less than 2nd field - correct
                %detection
                if placeField_filtered_max_posnorm{ss}{tt}{rr}(1) < placeField_filtered_max_posnorm{ss}{tt}{rr}(2)
                events_in_field{ss}{tt}{rr} = find(event_norm_pos_run{ss}.A{Aonly_field_filtered{ss}(rr)} >= placeField_filtered_max_posnorm{ss}{tt}{rr}(1) & ...
                    event_norm_pos_run{ss}.A{Aonly_field_filtered{ss}(rr)} <= placeField_filtered_max_posnorm{ss}{tt}{rr}(2));
                else %do two event finds along the split field along the start/edge of track 
                    events_in_field_temp_1 = find(event_norm_pos_run{ss}.A{Aonly_field_filtered{ss}(rr)} >= 0 & ...
                        event_norm_pos_run{ss}.A{Aonly_field_filtered{ss}(rr)} <= placeField_filtered_max_posnorm{ss}{tt}{rr}(2));
                    
                    events_in_field_temp_2 = find(event_norm_pos_run{ss}.A{Aonly_field_filtered{ss}(rr)} >= placeField_filtered_max_posnorm{ss}{tt}{rr}(1) & ...
                        event_norm_pos_run{ss}.A{Aonly_field_filtered{ss}(rr)} <= 1 );
                    
                    %merge and sort indices here
                    events_in_field{ss}{tt}{rr} = sort([events_in_field_temp_1; events_in_field_temp_2]);
                end
                
                %register the corresponding lap of in-field filtered event
                event_in_field_laps{ss}{tt}{rr} = event_lap_idx{ss}.A{Aonly_field_filtered{ss}(rr)}(events_in_field{ss}{tt}{rr});
                %get number of unique events (those occuring on each lap)
                event_in_field_nb{ss}{tt}{rr} = size(unique(event_in_field_laps{ss}{tt}{rr}),1);
                %get position of in-field events
                events_in_field_pos{ss}{tt}{rr} = event_norm_pos_run{ss}.A{Aonly_field_filtered{ss}(rr)}(events_in_field{ss}{tt}{rr});
                
            elseif tt == selectTrial(2) %correct B trials or all B
                %check if first edge is less than 2nd field - correct
                %detection                
                if placeField_filtered_max_posnorm{ss}{tt}{rr}(1) < placeField_filtered_max_posnorm{ss}{tt}{rr}(2)
                events_in_field{ss}{tt}{rr} = find(event_norm_pos_run{ss}.B{Bonly_field_filtered{ss}(rr)} >= placeField_filtered_max_posnorm{ss}{tt}{rr}(1) & ...
                    event_norm_pos_run{ss}.B{Bonly_field_filtered{ss}(rr)} <= placeField_filtered_max_posnorm{ss}{tt}{rr}(2));
                else
                    events_in_field_temp_1 = find(event_norm_pos_run{ss}.B{Bonly_field_filtered{ss}(rr)} >= 0 & ...
                        event_norm_pos_run{ss}.B{Bonly_field_filtered{ss}(rr)} <= placeField_filtered_max_posnorm{ss}{tt}{rr}(2));
                    
                    events_in_field_temp_2 = find(event_norm_pos_run{ss}.B{Bonly_field_filtered{ss}(rr)} >= placeField_filtered_max_posnorm{ss}{tt}{rr}(1) & ...
                        event_norm_pos_run{ss}.B{Bonly_field_filtered{ss}(rr)} <= 1 );
                    
                    %merge and sort indices here
                    events_in_field{ss}{tt}{rr} = sort([events_in_field_temp_1; events_in_field_temp_2]);                    
                end
                
                %register the corresponding lap of in-field filtered event
                event_in_field_laps{ss}{tt}{rr} = event_lap_idx{ss}.B{Bonly_field_filtered{ss}(rr)}(events_in_field{ss}{tt}{rr});
                %get number of unique events (those occuring on each lap)
                event_in_field_nb{ss}{tt}{rr} = size(unique(event_in_field_laps{ss}{tt}{rr}),1);
                %get position of in-field events
                events_in_field_pos{ss}{tt}{rr} = event_norm_pos_run{ss}.B{Bonly_field_filtered{ss}(rr)}(events_in_field{ss}{tt}{rr});
            end
        end
    end
    
    
    
    %check which task-selective ROIs have less than 5 events
    %correct A trials
    event_thres_exclude_log{ss}.A  = cell2mat(event_in_field_nb{ss}{selectTrial(1)}) < 5;
    %correct B trials
    event_thres_exclude_log{ss}.B  = cell2mat(event_in_field_nb{ss}{selectTrial(2)}) < 5;

    
    %update indices with event
    ROI_field_filtered_event{ss}.A = Aonly_field_filtered{ss}(~event_thres_exclude_log{ss}.A);
    ROI_field_filtered_event{ss}.B = Bonly_field_filtered{ss}(~event_thres_exclude_log{ss}.B);
    
    %update assn place fields
    placeField_eventFilt{ss}{selectTrial(1)} = placeField_filtered_max_posnorm{ss}{selectTrial(1)}(~event_thres_exclude_log{ss}.A);
    placeField_eventFilt{ss}{selectTrial(2)} = placeField_filtered_max_posnorm{ss}{selectTrial(2)}(~event_thres_exclude_log{ss}.B);
    
    %update event position (normalized)
    event_pos_inField{ss}{selectTrial(1)} = events_in_field_pos{ss}{selectTrial(1)}(~event_thres_exclude_log{ss}.A);
    event_pos_inField{ss}{selectTrial(2)} = events_in_field_pos{ss}{selectTrial(2)}(~event_thres_exclude_log{ss}.B);
end

%% Make sure the animal was in a run epoch in the min/max range of space on opposing laps (al least 80%) of space on at least 6 laps

% Check that in running epoch within 3 bins to the left or right of each
% ROI_field_filtered_event.A
% ROI_field_filtered_event.B

for ss=sessionSelect
    %take the median position of the events in field and min/max position
    for tt=selectTrial %for correct A and B trials
        for rr=1:size(event_pos_inField{ss}{tt},2)
            med_pos_event{ss}{tt}(rr) = median(event_pos_inField{ss}{tt}{rr});
            %into one matrix min and max of each event
            min_max_pos_event{ss}{tt}(rr,1) = min(event_pos_inField{ss}{tt}{rr});
            min_max_pos_event{ss}{tt}(rr,2) = max(event_pos_inField{ss}{tt}{rr});
        end
    end
    
    %get correct A and B laps idx
    corr_lap_idx{ss}{selectTrial(1)} = unique(Behavior_split{ss}{selectTrial(1)}.resampled.lapNb);
    corr_lap_idx{ss}{selectTrial(2)} = unique(Behavior_split{ss}{selectTrial(2)}.resampled.lapNb);
    
    
    
    %get indices across all laps of the ranges
    for tt=selectTrial %for correct A and B trials
        for rr=1:size(event_pos_inField{ss}{tt},2)
            %get indices that match the position range (ALL LAPS)
            pos_range_indices{ss}{tt}{rr} = find( Behavior_full{ss}.resampled.normalizedposition >= min_max_pos_event{ss}{tt}(rr,1) & ...
                Behavior_full{ss}.resampled.normalizedposition <= min_max_pos_event{ss}{tt}(rr,2));
            %get lap idx of corresponding idxs
            lap_idx_range{ss}{tt}{rr} = Behavior_full{ss}.resampled.lapNb(pos_range_indices{ss}{tt}{rr});
        end
    end
    
    
    %for correct A or B trials
    for tt=selectTrial
        %for each ROI in correct A trials
        for rr=1:size(lap_idx_range{ss}{tt},2)
            %extact logical with only laps correponding to opposing trial laps- dependent on B parameter
            if tt == selectTrial(1) %if looking on run status in B trials for correct A trials
                lap_opposed_idx{ss}{tt}{rr} = ismember(lap_idx_range{ss}{tt}{rr},corr_lap_idx{ss}{selectTrial(2)});
            elseif tt == selectTrial(2) %if looking on run status in A trials for correct B trials
                lap_opposed_idx{ss}{tt}{rr} = ismember(lap_idx_range{ss}{tt}{rr},corr_lap_idx{ss}{selectTrial(1)});
            end
            %get the lap number associated with each frame in the opposing trials
            lap_label_opposed{ss}{tt}{rr} = lap_idx_range{ss}{tt}{rr}(lap_opposed_idx{ss}{tt}{rr});
            %the binary indicating in animal in run epoch with that range
            lap_runEpoch_opposed{ss}{tt}{rr} = Behavior_full{ss}.run_ones(pos_range_indices{ss}{tt}{rr}(lap_opposed_idx{ss}{tt}{rr}));
            %extract the associated positions
            lap_pos_opposed{ss}{tt}{rr} = Behavior_full{ss}.resampled.normalizedposition(pos_range_indices{ss}{tt}{rr}(lap_opposed_idx{ss}{tt}{rr}));
        end
    end
    
    %split into individual laps for A events, look in B laps
    %for correct A or B trials
    for tt=selectTrial
        %for each ROI in correct A trials
        for rr=1:size(lap_idx_range{ss}{tt},2)
            %for each opposing lap
            if tt == selectTrial(1) %if looking on run status in B trials for correct A trials
                for ll=1:size(corr_lap_idx{ss}{selectTrial(2)},1)
                    % - depdendent on B parameter
                    split_lap_idxs{ss}{tt}{rr}{ll} = find(lap_label_opposed{ss}{tt}{rr} == corr_lap_idx{ss}{selectTrial(2)}(ll));
                    
                    split_lap_pos{ss}{tt}{rr}{ll} = lap_pos_opposed{ss}{tt}{rr}(split_lap_idxs{ss}{tt}{rr}{ll});
                    split_lap_runEpoch{ss}{tt}{rr}{ll} = lap_runEpoch_opposed{ss}{tt}{rr}(split_lap_idxs{ss}{tt}{rr}{ll});
                end
            elseif tt == selectTrial(2)
                for ll=1:size(corr_lap_idx{ss}{selectTrial(1)},1)
                    split_lap_idxs{ss}{tt}{rr}{ll} = find(lap_label_opposed{ss}{tt}{rr} == corr_lap_idx{ss}{selectTrial(1)}(ll));
                    
                    split_lap_pos{ss}{tt}{rr}{ll} = lap_pos_opposed{ss}{tt}{rr}(split_lap_idxs{ss}{tt}{rr}{ll});
                    split_lap_runEpoch{ss}{tt}{rr}{ll} = lap_runEpoch_opposed{ss}{tt}{rr}(split_lap_idxs{ss}{tt}{rr}{ll});
                end
            end
            
        end
    end



%for correct A or B trials
for tt=selectTrial
    %for each ROI in correct A trials
    for rr=1:size(lap_idx_range{ss}{tt},2)
        %for each lap
        for ll=1:size(split_lap_pos{ss}{tt}{rr},2)
            %check if empty on any laps
            if ~isempty(split_lap_pos{ss}{tt}{rr}{ll})
                unique_pos{ss}{tt}{rr}{ll} = unique(split_lap_pos{ss}{tt}{rr}{ll});
                %for each unique position, check if entirety in run epoch
                for pos_idx=1:size(unique_pos{ss}{tt}{rr}{ll},1)
                    %get idx associated with given unique pos in the lap
                    pos_idxs_each{ss}{tt}{rr}{ll}{pos_idx} = find(split_lap_pos{ss}{tt}{rr}{ll} == unique_pos{ss}{tt}{rr}{ll}(pos_idx));
                    run_Epoch_each{ss}{tt}{rr}{ll}{pos_idx} = split_lap_runEpoch{ss}{tt}{rr}{ll}(pos_idxs_each{ss}{tt}{rr}{ll}{pos_idx});
                    %check which position did animal spend all of it in run state
                    
                end
            else
                %can't be in run epoch in position not id'd (rare)
                run_Epoch_each{ss}{tt}{rr}{ll}{1} = 0;
            end
            %calculate logical of  run at each position and get fraction in run
            %state for each lap
            frac_run_lap{ss}{tt}{rr}(ll) = sum(cellfun(@prod,run_Epoch_each{ss}{tt}{rr}{ll}))/size(run_Epoch_each{ss}{tt}{rr}{ll},2);
        end
    end
end



%for all laps check how above 80% of space in run epoch and check if this
%occurs in at least 6 laps
%for correct A or B trials
for tt=selectTrial
    %for each ROI in correct A trials
    for rr=1:size(lap_idx_range{ss}{tt},2)
        if sum(frac_run_lap{ss}{tt}{rr} >= 0.8) >= 6
            %generate logical with 1's include ROI, and 0 exclude ROI
            run_epoch_filt_include_log{ss}{tt}(rr) = 1;
        else
            run_epoch_filt_include_log{ss}{tt}(rr) = 0;
        end
    end
end

%apply final filter to ROI indices
%select ROI indices (from original)
final_filtered_ROI{ss}.A = ROI_field_filtered_event{ss}.A(logical(run_epoch_filt_include_log{ss}{selectTrial(1)}));
final_filtered_ROI{ss}.B = ROI_field_filtered_event{ss}.B(logical(run_epoch_filt_include_log{ss}{selectTrial(2)}));


%place field width positions
placeField_final{ss}{selectTrial(1)} = placeField_eventFilt{ss}{selectTrial(1)}(logical(run_epoch_filt_include_log{ss}{selectTrial(1)}));
placeField_final{ss}{selectTrial(2)} = placeField_eventFilt{ss}{selectTrial(2)}(logical(run_epoch_filt_include_log{ss}{selectTrial(2)}));

%event position (normalized)
event_pos_inField_final{ss}{selectTrial(1)} = event_pos_inField{ss}{selectTrial(1)}(logical(run_epoch_filt_include_log{ss}{selectTrial(1)}));
event_pos_inField_final{ss}{selectTrial(2)} = event_pos_inField{ss}{selectTrial(2)}(logical(run_epoch_filt_include_log{ss}{selectTrial(2)}));
end

%% Plot as shaded area to verify correct id of place field onto normalized

if options.dispFigure ==1
    %plot normalized position
    %show A selective first
    figure('Position', [1930 130 1890 420])
    hold on
    title('A-selective filtered')
    for rr=1:size(final_filtered_ROI{ss}.A,2)%1:size(Aonly_notSIb_idx,2)
        %ROI = rr%AandB_tuned_idx(rr);
        ROI = final_filtered_ROI{ss}.A(rr); %Aonly_notSIb_idx(rr);
        hold on
        title(num2str(ROI))
        yticks([0 0.5 1])
        ylabel('Normalized position')
        xlabel('Time [min]');
        xticks(0:3:12);
        %ylim([0 1])
        set(gca,'FontSize',14)
        set(gca,'LineWidth',1)
        %A laps
        for ii=1:size(lap_idxs{ss}.A,1)
            plot(Imaging_split{ss}{selectTrial(1)}.time_restricted(lap_idxs{ss}.A(ii,1):lap_idxs{ss}.A(ii,2))/60,...
                Behavior_split{ss}{selectTrial(1)}.resampled.position_norm(lap_idxs{ss}.A(ii,1):lap_idxs{ss}.A(ii,2)),...
                'Color',[0 0 1 0.6],'LineWidth',1.5)
        end
        %B laps
        for ii=1:size(lap_idxs{ss}.B,1)
            plot(Imaging_split{ss}{selectTrial(2)}.time_restricted(lap_idxs{ss}.B(ii,1):lap_idxs{ss}.B(ii,2))/60,...
                Behavior_split{ss}{selectTrial(2)}.resampled.position_norm(lap_idxs{ss}.B(ii,1):lap_idxs{ss}.B(ii,2)),...
                'Color',[1 0 0 0.6],'LineWidth',1.5)
        end
        %overlay significant calcium run events
        %A
        scatter(event_norm_time{ss}.A{ROI},event_norm_pos_run{ss}.A{ROI},[],[0 0 1],'*')
        %B
        scatter(event_norm_time{ss}.B{ROI},event_norm_pos_run{ss}.B{ROI},[],[1 0 0],'*')
        
        %plot horz lines signifying start and end of place field
        %start
        %lineS = refline(0,placeField_final{ss}{rr}(selectTrial(1)))
        %lineS.Color = 'g';
        %end
        %lineE = refline(0,placeField_final{ss}{rr}(selectTrial(2)))
        %lineE.Color = 'g';
        
        pause
        clf
    end
    
%FIX THIS
    %show B selective second
    figure('Position', [1930 130 1890 420])
    hold on
    title('B-selective filtered')
    for rr=1:size(final_filtered_ROI{ss}.B,2)%1:size(Aonly_notSIb_idx,2)
        %ROI = rr%AandB_tuned_idx(rr);
        ROI = final_filtered_ROI{ss}.B(rr); %Aonly_notSIb_idx(rr);
        hold on
        title(num2str(ROI))
        yticks([0 0.5 1])
        ylabel('Normalized position')
        xlabel('Time [min]');
        xticks(0:3:12);
        %ylim([0 1])
        set(gca,'FontSize',14)
        set(gca,'LineWidth',1)
        %A laps
        for ii=1:size(lap_idxs{ss}.A,1)
            plot(Imaging_split{ss}{selectTrial(1)}.time_restricted(lap_idxs{ss}.A(ii,1):lap_idxs{ss}.A(ii,2))/60,...
                Behavior_split{ss}{selectTrial(1)}.resampled.position_norm(lap_idxs{ss}.A(ii,1):lap_idxs{ss}.A(ii,2)),...
                'Color',[0 0 1 0.6],'LineWidth',1.5)
        end
        %B laps
        for ii=1:size(lap_idxs{ss}.B,1)
            plot(Imaging_split{ss}{selectTrial(2)}.time_restricted(lap_idxs{ss}.B(ii,1):lap_idxs{ss}.B(ii,2))/60,...
                Behavior_split{ss}{selectTrial(2)}.resampled.position_norm(lap_idxs{ss}.B(ii,1):lap_idxs{ss}.B(ii,2)),...
                'Color',[1 0 0 0.6],'LineWidth',1.5)
        end
        %overlay significant calcium run events
        %A
        scatter(event_norm_time{ss}.A{ROI},event_norm_pos_run{ss}.A{ROI},[],[0 0 1],'*')
        %B
        scatter(event_norm_time{ss}.B{ROI},event_norm_pos_run{ss}.B{ROI},[],[1 0 0],'*')
        
        %plot horz lines signifying start and end of place field
        %start
        %lineS = refline(0,placeField_final{2}{rr}(1))
        lineS.Color = 'g';
        %end
        %lineE = refline(0,placeField_final{2}{rr}(2))
        lineE.Color = 'g';
        
        pause
        clf
    end
end

%% Export task-selective ROIs in struct
for ss=sessionSelect
    %export indices for task selective neurons
    task_selective_ROIs{ss}.A.idx = final_filtered_ROI{ss}.A;
    task_selective_ROIs{ss}.B.idx = final_filtered_ROI{ss}.B;
    
    %field margin for task selective neurons
    task_selective_ROIs{ss}.A.field_margin = placeField_final{ss}{selectTrial(1)};
    task_selective_ROIs{ss}.B.field_margin = placeField_final{ss}{selectTrial(2)};
    
    %event associated with selected place field
    task_selective_ROIs{ss}.A.fieldEvents = event_pos_inField_final{ss}{selectTrial(1)};
    task_selective_ROIs{ss}.B.fieldEvents = event_pos_inField_final{ss}{selectTrial(2)};
end


end

