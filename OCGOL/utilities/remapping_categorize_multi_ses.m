function [remapping_ROIs,partial_field_idx] = remapping_categorize_multi_ses(cent_diff, tuned_logical,pf_vector, session_vars,max_transient_peak, pf_count_filtered,select_fields, options)
%split mutually tuned neurons by remapping category: 
%common (less than certain centroid difference between max
%tuned
partial_field_idx = 0;

%% Define input variables

sessionSelect = options.sessionSelect;
selectTrial = options.selectTrial;

%% Set parameters

%number of degrees of centroid difference
deg_thres = options.deg_thres;
%degree range for splitting the global remapping neurons
deg_ranges = options.deg_ranges;
%significance of level for AUC comparison test
p_sig = options.p_sig;

%% Get ROI indices of A selective, B selective, A&B selective neurons

for ss=sessionSelect
    %find A&B tuned by either criterion (input option below)
    %unclassified neurons get moved into mixed category
    tuned_A_si_ts{ss} = tuned_logical(ss).si.Atuned  | tuned_logical(ss).ts.Atuned;
    tuned_B_si_ts{ss} = tuned_logical(ss).si.Btuned  | tuned_logical(ss).ts.Btuned;
    tuned_AB_si_ts{ss} = tuned_A_si_ts{ss} & tuned_B_si_ts{ss};
    
    %rate map, common, and global
    %only TS tuned neurons
    rate_remap_idx_start{ss} = find(tuned_logical(ss).ts.AandB_tuned ==1);
    %rate_remap_idx_start = find(tuned_AB_si_ts == 1);
    
    %partial remapping - both A and B tuned to SI
    %partial_remap_idx_start = find(tuned_logical.si.AandB_tuned == 1);
    partial_remap_idx_start{ss} = find(tuned_AB_si_ts{ss} == 1);
end

for ss=sessionSelect
    %choose if SI or TS tuned
    switch options.tuning_criterion
        case 'si' %spatial information
            AandB_tuned_idx{ss} = find(tuned_logical(ss).si.AandB_tuned == 1);
            Aonly_tuned_idx{ss} = find(tuned_logical(ss).si.onlyA_tuned == 1);
            Bonly_tuned_idx{ss} = find(tuned_logical(ss).si.onlyB_tuned == 1);
            
        case 'ts' %tuning specificity
            %both A and B tuned by TS
            AandB_tuned_idx{ss} = find(tuned_logical(ss).ts.AandB_tuned ==1);
            %only A tuned by TS
            Aonly_tuned_idx{ss} = find(tuned_logical(ss).ts.onlyA_tuned == 1);
            %only B tuned by tS
            Bonly_tuned_idx{ss} = find(tuned_logical(ss).ts.onlyB_tuned == 1);
            
            %all A tuned by SI
            A_tuned_si_idx{ss} = find(tuned_logical(ss).si.Atuned == 1);
            %all B tuned by SI
            B_tuned_si_idx{ss} = find(tuned_logical(ss).si.Btuned == 1);
            
            %only A tuned by TS, but also not SI B tuned
            Aonly_notSIb_idx{ss} = setdiff(Aonly_tuned_idx{ss},B_tuned_si_idx{ss});
            
            %only B tuned by TS, but also not SI A tuned
            Bonly_notSIa_idx{ss} = setdiff(Bonly_tuned_idx{ss},A_tuned_si_idx{ss});
            
    end
end

%% Filter make sure rate/global/common remapping category has single PF in each and meets min event if field criteriym

for ss=sessionSelect
    %use place field calculation and logical from place properties calculation
    single_pf_idxs{ss} = find(sum((pf_count_filtered{ss}(selectTrial,:) == 1),1) ==2);
    %make copy of original rate remap vector
    rate_remap_idx_orig{ss} = rate_remap_idx_start{ss};
    %extract only the indices
    rate_remap_idx_start{ss} = intersect(rate_remap_idx_orig{ss},single_pf_idxs{ss});
end

%% Parse starting partial remap neurons (2PF vs 1 PF)
for ss=sessionSelect
    %ROIs with siginificant 2PF vs. 1PFs
    single_double_pf_idxs{ss} = find(sum((pf_count_filtered{ss}(selectTrial,:) == 1) + 2*(pf_count_filtered{ss}(selectTrial,:) == 2),1) ==3);
    %make copy of original rate remap vector
    partial_remap_idx_orig{ss} = partial_remap_idx_start{ss};
    %extract only the indices that have 2PF vs. 1 PF place fields
    partial_remap_idx_start{ss} = intersect(partial_remap_idx_orig{ss},single_double_pf_idxs{ss});
end

%% Screen all both trial matching idxs to make at least 1 PF/min 5 event in field
for ss=sessionSelect
    at_least_one_pf_idxs{ss} = find(sum((pf_count_filtered{ss}(selectTrial,:) ~= 0),1) ==2);
    tuned_AB_si_ts_filt_idx{ss} = intersect(find(tuned_AB_si_ts{ss} ==1),at_least_one_pf_idxs{ss});
end

%% Truncate the centroid difference struct to only include the selected ROIs
for ss=sessionSelect
    cent_diff_AandB{ss}.angle_diff = cent_diff.angle_diff{ss}(rate_remap_idx_start{ss});
    cent_diff_AandB{ss}.max_bin = cent_diff.angle_diff{ss}(:,rate_remap_idx_start{ss});
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

%% Get lap indices for each lap in all A or B trials

%only correct
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
%trials chosen with seletTrial option

for ss=sessionSelect
    %for each ROI
    for rr=1:size(events{ss}{1},2)
        %time of significant run events in A
        event_norm_time{ss}.A{rr} = Imaging_split{ss}{selectTrial(1)}.time_restricted(find(Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr) == 1))/60;
        %normalizesd position of significant run events in A
        event_norm_pos_run{ss}.A{rr} = Behavior_split{ss}{selectTrial(1)}.resampled.position_norm(find(Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr) == 1));
        %lap assignment for A
        event_lap_idx{ss}.A{rr} = Behavior_split{ss}{selectTrial(1)}.resampled.lapNb(logical(Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr)));
        %get respective AUC values
        event_AUC{ss}.A{rr} = Event_split{ss}{selectTrial(1)}.Run.properties.AUC{rr};
        
        %time of significant run events in B
        event_norm_time{ss}.B{rr} = Imaging_split{ss}{selectTrial(2)}.time_restricted(find(Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr) == 1))/60;
        %normalizesd position of significant run events in B
        event_norm_pos_run{ss}.B{rr} = Behavior_split{ss}{selectTrial(2)}.resampled.position_norm(find(Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr) == 1));
        %lap assignment for B
        event_lap_idx{ss}.B{rr} = Behavior_split{ss}{selectTrial(2)}.resampled.lapNb(logical(Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr)));
        %get respective AUC values
        event_AUC{ss}.B{rr} = Event_split{ss}{selectTrial(2)}.Run.properties.AUC{rr};
        
    end
end

%% Plot the run epochs, corresponding position bin edges of place field

%overlay area of max place field bin

%TODO: MODIFY THIS TO TAKE IN THE CHOSEN TRIALS TO
%USE VS. FIXED AS IS NOW

%spatial bin assignment for each run-epoch frame (100 bins) 
% %A trials
% session_vars{1}.Place_cell{1}.Bin{8};
% %B trials
% session_vars{1}.Place_cell{2}.Bin{8};
% %get the run epoch binaries for each set of trials - A session - 5339 run
% %frames
% %plot the run epochs as patches
% session_vars{1}.Behavior_split{1}.run_ones;
% %get the run epoch binaries for each set of trials - B session - 4547 run
% %frames
% session_vars{1}.Behavior_split{2}.run_ones; 
% %run events - binary onset across entire run/no run interval
% Event_split{1}{1}.Run.run_onset_ones;

for ss=sessionSelect
    %get bins
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


%% Remove ROIs idx's without a id'd place field _ THIS SHOULD BE REDUNDANT GIVEN FILTERING DONE ABOVE

for ss=sessionSelect
    %edges of all identified common TS neurons
    %correct A
    placeFieldEdges{ss}{selectTrial(1)} = Place_cell{ss}{selectTrial(1)}.placeField.edge(rate_remap_idx_start{ss});
    %correct B
    placeFieldEdges{ss}{selectTrial(2)} = Place_cell{ss}{selectTrial(2)}.placeField.edge(rate_remap_idx_start{ss});
    
    %in A trials
    idx_wo_placeFields{ss}{selectTrial(1)} = rate_remap_idx_start{ss}(find(cellfun(@isempty,placeFieldEdges{ss}{selectTrial(1)}) ==1));
    %in B trials
    idx_wo_placeFields{ss}{selectTrial(2)} = rate_remap_idx_start{ss}(find(cellfun(@isempty,placeFieldEdges{ss}{selectTrial(2)}) ==1));
    
    %merge empty field indices into 1 vector for removal from list
    rm_idx_no_pf{ss} = unique(cell2mat(idx_wo_placeFields{ss}(:,selectTrial)));
    
    %remove ROIs for A and B that do not have place fields
    if ~isempty(rm_idx_no_pf{ss})
        rate_remap_idx_filtered{ss} = setdiff(rate_remap_idx_start{ss},rm_idx_no_pf{ss});
        
    else %keep the previous ROIs (copy only)
        rate_remap_idx_filtered{ss} = rate_remap_idx_start{ss};
    end
    
    %filter out centroid difference data based on no-id'd neuron tuned by TS
    %find indices of removed neurons from list
    select_pf_filtered_log{ss} = ismember(rate_remap_idx_start{ss},rate_remap_idx_filtered{ss});
    
    %update centroid input data with pf filtered idxs
    cent_diff_AandB_pf_filt{ss}.angle_diff = cent_diff_AandB{ss}.angle_diff(select_pf_filtered_log{ss});
    cent_diff_AandB_pf_filt{ss}.max_bin = cent_diff_AandB{ss}.max_bin(:,select_pf_filtered_log{ss});
    
    %get the ROI Idxs associated with place filtered ROIS
    remapping_pf_filtered{ss} = rate_remap_idx_start{ss}(select_pf_filtered_log{ss});
    
end

%% Determine the equivalent normalized position range of the max place field

for ss=sessionSelect
    %extract the max place field index for A and B trial using filtered A tuned
    %and B tuned ROIs
    max_field_idx{ss}{selectTrial(1)} = max_transient_peak{ss}{selectTrial(1)}(remapping_pf_filtered{ss});
    max_field_idx{ss}{selectTrial(2)} = max_transient_peak{ss}{selectTrial(2)}(remapping_pf_filtered{ss});
    
    %get the edges of the max transient place field for each set of idxs
    %edges of all identified
    %correct A
    placeField_filtered{ss}{selectTrial(1)} = Place_cell{ss}{selectTrial(1)}.placeField.edge(remapping_pf_filtered{ss});
    %correct B
    placeField_filtered{ss}{selectTrial(2)} = Place_cell{ss}{selectTrial(2)}.placeField.edge(remapping_pf_filtered{ss});
    
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
    median_track_len{ss} = median(Behavior_full{ss}.position_lap(:,2));
    
    %conversion factor (norm_pos/cm length) - 100 bins
    norm_conv_factor{ss} = median_track_len{ss}/100;
end

%% Filter out neurons that do not have at least 5 sig events in max place field in at least 5 distinct laps

%calcium event position and absolute restrict time
%all neuron idx space
%relevant variables
% event_norm_pos_run.A;
% event_norm_time.A;  
% event_lap_idx.A;
% 
% Aonly_field_filtered;
% Bonly_field_filtered;
for ss=sessionSelect
    %find events occuring within max place field for each ROI
    for tt=selectTrial
        for rr=1:size(placeField_filtered_max_posnorm{ss}{tt},2)
            %get idxs of events with max place field
            %need this selector because of struct assignment of A vs.B
            if tt == selectTrial(1) %correct A trials
                %OLDER VERSIONS - REPLACED BY CODE BELOW
                %events_in_field{tt}{rr} = find(event_norm_pos_run.A{remapping_pf_filtered(rr)} >= placeField_filtered_max_posnorm{tt}{rr}(1) & ...
                %   event_norm_pos_run.A{remapping_pf_filtered(rr)} <= placeField_filtered_max_posnorm{tt}{rr}(2));
                
                %check of the normalized field position edges cross the start of the track
                %(first edge) - i.e. if first edge is greater than
                %the second edge
                if placeField_filtered_max_posnorm{ss}{tt}{rr}(1) < placeField_filtered_max_posnorm{ss}{tt}{rr}(2)
                    events_in_field{ss}{tt}{rr} = find(event_norm_pos_run{ss}.A{remapping_pf_filtered{ss}(rr)} >= placeField_filtered_max_posnorm{ss}{tt}{rr}(1) & ...
                        event_norm_pos_run{ss}.A{remapping_pf_filtered{ss}(rr)} <= placeField_filtered_max_posnorm{ss}{tt}{rr}(2));
                else %do two event finds and merge (and sort) into one set of indices
                    events_in_field_temp_1 = find(event_norm_pos_run{ss}.A{remapping_pf_filtered{ss}(rr)} >= 0 & ...
                        event_norm_pos_run{ss}.A{remapping_pf_filtered{ss}(rr)} <= placeField_filtered_max_posnorm{ss}{tt}{rr}(2));
                    %typo was here replace with correct edge
                    events_in_field_temp_2 = find(event_norm_pos_run{ss}.A{remapping_pf_filtered{ss}(rr)} >= placeField_filtered_max_posnorm{ss}{tt}{rr}(1) & ...
                        event_norm_pos_run{ss}.A{remapping_pf_filtered{ss}(rr)} <= 1 );
                    
                    %merge and sort indices here
                    events_in_field{ss}{tt}{rr} = sort([events_in_field_temp_1; events_in_field_temp_2]);
                end
                
                
                %register the corresponding lap of in-field filtered event
                event_in_field_laps{ss}{tt}{rr} = event_lap_idx{ss}.A{remapping_pf_filtered{ss}(rr)}(events_in_field{ss}{tt}{rr});
                %get number of unique events (those occuring on each lap)
                event_in_field_nb{ss}{tt}{rr} = size(unique(event_in_field_laps{ss}{tt}{rr}),1);
                %get position of in-field events
                events_in_field_pos{ss}{tt}{rr} = event_norm_pos_run{ss}.A{remapping_pf_filtered{ss}(rr)}(events_in_field{ss}{tt}{rr});
                
            elseif tt == selectTrial(2) %correct B trials
                %events_in_field{tt}{rr} = find(event_norm_pos_run.B{remapping_pf_filtered(rr)} >= placeField_filtered_max_posnorm{tt}{rr}(1) & ...
                %    event_norm_pos_run.B{remapping_pf_filtered(rr)} <= placeField_filtered_max_posnorm{tt}{rr}(2));
                
                %check of the normalized field position edges cross the start of the track
                %(first edge) - i.e. if first edge is greater than
                %the second edge
                if placeField_filtered_max_posnorm{ss}{tt}{rr}(1) < placeField_filtered_max_posnorm{ss}{tt}{rr}(2)
                    events_in_field{ss}{tt}{rr} = find(event_norm_pos_run{ss}.B{remapping_pf_filtered{ss}(rr)} >= placeField_filtered_max_posnorm{ss}{tt}{rr}(1) & ...
                        event_norm_pos_run{ss}.B{remapping_pf_filtered{ss}(rr)} <= placeField_filtered_max_posnorm{ss}{tt}{rr}(2));
                else %do two event finds and merge (and sort) into one set of indices
                    events_in_field_temp_1 = find(event_norm_pos_run{ss}.B{remapping_pf_filtered{ss}(rr)} >= 0 & ...
                        event_norm_pos_run{ss}.B{remapping_pf_filtered{ss}(rr)} <= placeField_filtered_max_posnorm{ss}{tt}{rr}(2));
                    %typo was below - replaced with correct edge
                    events_in_field_temp_2 = find(event_norm_pos_run{ss}.B{remapping_pf_filtered{ss}(rr)} >= placeField_filtered_max_posnorm{ss}{tt}{rr}(1) & ...
                        event_norm_pos_run{ss}.B{remapping_pf_filtered{ss}(rr)} <= 1 );
                    
                    %merge and sort indices here
                    events_in_field{ss}{tt}{rr} = sort([events_in_field_temp_1; events_in_field_temp_2]);
                end
                
                %register the corresponding lap of in-field filtered event
                event_in_field_laps{ss}{tt}{rr} = event_lap_idx{ss}.B{remapping_pf_filtered{ss}(rr)}(events_in_field{ss}{tt}{rr});
                %get number of unique events (those occuring on each lap)
                event_in_field_nb{ss}{tt}{rr} = size(unique(event_in_field_laps{ss}{tt}{rr}),1);
                %get position of in-field events
                events_in_field_pos{ss}{tt}{rr} = event_norm_pos_run{ss}.B{remapping_pf_filtered{ss}(rr)}(events_in_field{ss}{tt}{rr});
            end
        end
    end
end

for ss=sessionSelect
    %check which task-selective ROIs have less than 5 events
    %correct A trials
    event_thres_exclude_log{ss}.A  = cell2mat(event_in_field_nb{ss}{selectTrial(1)}) < 5;
    %correct B trials
    event_thres_exclude_log{ss}.B  = cell2mat(event_in_field_nb{ss}{selectTrial(2)}) < 5;
    
    %merge the excludes neurosn idx
    event_thres_excludeAB_log{ss} =  event_thres_exclude_log{ss}.A | event_thres_exclude_log{ss}.B;
    
    %update indices with event - update absolute idxs
    ROI_field_filtered_event{ss}.A = remapping_pf_filtered{ss}(~event_thres_excludeAB_log{ss});
    ROI_field_filtered_event{ss}.B = remapping_pf_filtered{ss}(~event_thres_excludeAB_log{ss});
    
    %select the event indices for event/place filtered ROIs
    events_in_field_pf_filtered{ss}{selectTrial(1)} = events_in_field{ss}{selectTrial(1)}(~event_thres_excludeAB_log{ss});
    events_in_field_pf_filtered{ss}{selectTrial(2)} = events_in_field{ss}{selectTrial(2)}(~event_thres_excludeAB_log{ss});
    
    %update assn place fields
    placeField_eventFilt{ss}{selectTrial(1)} = placeField_filtered_max_posnorm{ss}{selectTrial(1)}(~event_thres_excludeAB_log{ss});
    placeField_eventFilt{ss}{selectTrial(2)} = placeField_filtered_max_posnorm{ss}{selectTrial(2)}(~event_thres_excludeAB_log{ss});
    
    %update event position (normalized)
    event_pos_inField{ss}{selectTrial(1)} = events_in_field_pos{ss}{selectTrial(1)}(~event_thres_excludeAB_log{ss});
    event_pos_inField{ss}{selectTrial(2)} = events_in_field_pos{ss}{selectTrial(2)}(~event_thres_excludeAB_log{ss});
    
    %update centroid input data with pf filtered idxs
    cent_diff_AandB_eventCount{ss}.angle_diff = cent_diff_AandB_pf_filt{ss}.angle_diff(~event_thres_excludeAB_log{ss});
    cent_diff_AandB_eventCount{ss}.max_bin = cent_diff_AandB_pf_filt{ss}.max_bin(:,~event_thres_excludeAB_log{ss});
end

%% Make sure the animal was in a run epoch in the min/max range of space on opposing laps (al least 80%) of space on at least 6 laps

% Check that in running epoch within 3 bins to the left or right of each
%ROI_field_filtered_event.A
%ROI_field_filtered_event.B
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
end

for ss=sessionSelect
%for correct A or B trials
for tt=selectTrial
    %for each ROI in correct A trials 
    for rr=1:size(lap_idx_range{ss}{tt},2)
        %extact logical with only laps correponding to opposing trial laps- dependent on B parameter
        if tt ==selectTrial(1) %if looking on run status in B trials for correct A trials
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
        if tt==selectTrial(1)
            for ll=1:size(corr_lap_idx{ss}{selectTrial(2)},1)
                % - depdendent on B parameter
                %if tt == 1 %if looking on run status in B trials for correct A trials
                split_lap_idxs{ss}{tt}{rr}{ll} = find(lap_label_opposed{ss}{tt}{rr} == corr_lap_idx{ss}{selectTrial(2)}(ll));
                %elseif tt == 2
                %split_lap_idxs{tt}{rr}{ll} = find(lap_label_opposed{tt}{rr} == corr_lap_idx{1}(ll));
                %end
                split_lap_pos{ss}{tt}{rr}{ll} = lap_pos_opposed{ss}{tt}{rr}(split_lap_idxs{ss}{tt}{rr}{ll});
                split_lap_runEpoch{ss}{tt}{rr}{ll} = lap_runEpoch_opposed{ss}{tt}{rr}(split_lap_idxs{ss}{tt}{rr}{ll});
            end
        elseif tt ==selectTrial(2)
            for ll=1:size(corr_lap_idx{ss}{selectTrial(1)},1)
                % - depdendent on B parameter
                %if tt == 1 %if looking on run status in B trials for correct A trials
                %split_lap_idxs{tt}{rr}{ll} = find(lap_label_opposed{tt}{rr} == corr_lap_idx{2}(ll));
                %elseif tt == 2
                split_lap_idxs{ss}{tt}{rr}{ll} = find(lap_label_opposed{ss}{tt}{rr} == corr_lap_idx{ss}{selectTrial(1)}(ll));
                %end
                split_lap_pos{ss}{tt}{rr}{ll} = lap_pos_opposed{ss}{tt}{rr}(split_lap_idxs{ss}{tt}{rr}{ll});
                split_lap_runEpoch{ss}{tt}{rr}{ll} = lap_runEpoch_opposed{ss}{tt}{rr}(split_lap_idxs{ss}{tt}{rr}{ll});
            end
        end
        
    end
end
end

for ss=sessionSelect
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
    disp(ss)
end

for ss=sessionSelect
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
end


for ss=sessionSelect
    %combine filters for each epoch into 1 filter
    run_epoch_filt_both{ss} = logical(run_epoch_filt_include_log{ss}{selectTrial(1)}) | logical(run_epoch_filt_include_log{ss}{selectTrial(2)});
    
    %apply final filter to ROI indices
    
    %select ROI indices (from original)
    final_filtered_ROI{ss}.A = ROI_field_filtered_event{ss}.A(run_epoch_filt_both{ss});
    final_filtered_ROI{ss}.B = ROI_field_filtered_event{ss}.B(run_epoch_filt_both{ss});
    
    %place field width positions
    placeField_final{ss}{selectTrial(1)} = placeField_eventFilt{ss}{selectTrial(1)}(run_epoch_filt_both{ss});
    placeField_final{ss}{selectTrial(2)} = placeField_eventFilt{ss}{selectTrial(2)}(run_epoch_filt_both{ss});
    
    %event position (normalized)
    event_pos_inField_final{ss}{selectTrial(1)} = event_pos_inField{ss}{selectTrial(1)}(run_epoch_filt_both{ss});
    event_pos_inField_final{ss}{selectTrial(2)} = event_pos_inField{ss}{selectTrial(2)}(run_epoch_filt_both{ss});
    
    %select the event indices for event/place filtered ROIs
    events_in_field_final_filtered{ss}{selectTrial(1)} = events_in_field_pf_filtered{ss}{selectTrial(1)}(run_epoch_filt_both{ss});
    events_in_field_final_filtered{ss}{selectTrial(2)} = events_in_field_pf_filtered{ss}{selectTrial(2)}(run_epoch_filt_both{ss});
    
    %update centroid input data with pf filtered idxs
    cent_diff_AandB_runFilt{ss}.angle_diff = cent_diff_AandB_eventCount{ss}.angle_diff(run_epoch_filt_both{ss});
    cent_diff_AandB_runFilt{ss}.max_bin = cent_diff_AandB_eventCount{ss}.max_bin(:,run_epoch_filt_both{ss});
end


%% Centroid difference filter (select those with cent diff less than 10cm ~ 18 deg) - point of split between rate and global remap
%centroid difference in degrees (36 deg ~ 20 cm) (27 deg ~ 15 cm)
%deg_thres = 45;
for ss=sessionSelect
    
    %get the vector idx's of the neurons with near and far centroid differences
    near_field_idx{ss} = find((cent_diff_AandB_runFilt{ss}.angle_diff < deg2rad(deg_ranges(2))) == 1);
    %mid_field_idx = find((cent_diff_AandB_runFilt.angle_diff <= deg2rad(deg_thres)) == 1)
    mid_field_idx{ss} =find(((cent_diff_AandB_runFilt{ss}.angle_diff >= deg2rad(deg_ranges(2))) &...
        (cent_diff_AandB_runFilt{ss}.angle_diff < deg2rad(deg_ranges(3)))) ==1);
    %far centroid split ROIs
    far_field_idx{ss} = find((cent_diff_AandB_runFilt{ss}.angle_diff >= deg2rad(deg_ranges(3))) == 1);
    
    %translate the idxs to absolute ROIs idxs (A and B are the same)
    near_idx_centroid{ss} = final_filtered_ROI{ss}.A(near_field_idx{ss});
    %second cateogory
    global_remap_ROI{ss}{2} = final_filtered_ROI{ss}.A(mid_field_idx{ss});
    %third category
    global_remap_ROI{ss}{3} = final_filtered_ROI{ss}.A(far_field_idx{ss});
    
end

for ss=sessionSelect
    %place field width positions
    %potential rate remapping
    placeField_near_centroid{ss}{selectTrial(1)} = placeField_final{ss}{selectTrial(1)}(near_field_idx{ss});
    placeField_near_centroid{ss}{selectTrial(2)} = placeField_final{ss}{selectTrial(2)}(near_field_idx{ss});
    
    %global remapping (mid and far field)
    placeField_global_centroid{ss}{2}{selectTrial(1)} = placeField_final{ss}{selectTrial(1)}(mid_field_idx{ss});
    placeField_global_centroid{ss}{2}{selectTrial(2)} = placeField_final{ss}{selectTrial(2)}(mid_field_idx{ss});
    
    placeField_global_centroid{ss}{3}{selectTrial(1)} = placeField_final{ss}{selectTrial(1)}(far_field_idx{ss});
    placeField_global_centroid{ss}{3}{selectTrial(2)} = placeField_final{ss}{selectTrial(2)}(far_field_idx{ss});
end

for ss=sessionSelect
    %event position (normalized)
    %potential rate remapping
    event_pos_inField_near_centroid{ss}{selectTrial(1)} = event_pos_inField_final{ss}{selectTrial(1)}(near_field_idx{ss});
    event_pos_inField_near_centroid{ss}{selectTrial(2)} = event_pos_inField_final{ss}{selectTrial(2)}(near_field_idx{ss});
    
    %global remapping
    event_pos_inField_global_centroid{ss}{2}{selectTrial(1)} = event_pos_inField_final{ss}{selectTrial(1)}(mid_field_idx{ss});
    event_pos_inField_global_centroid{ss}{2}{selectTrial(2)} = event_pos_inField_final{ss}{selectTrial(2)}(mid_field_idx{ss});
    
    event_pos_inField_global_centroid{ss}{3}{selectTrial(1)} = event_pos_inField_final{ss}{selectTrial(1)}(far_field_idx{ss});
    event_pos_inField_global_centroid{ss}{3}{selectTrial(2)} = event_pos_inField_final{ss}{selectTrial(2)}(far_field_idx{ss});
    
    %select the event indices for event/place filtered ROIs
    events_in_field_near_centroid{ss}{selectTrial(1)} = events_in_field_final_filtered{ss}{selectTrial(1)}(near_field_idx{ss});
    events_in_field_near_centroid{ss}{selectTrial(2)} = events_in_field_final_filtered{ss}{selectTrial(2)}(near_field_idx{ss});
    
    %global remapping
    events_in_field_global_centroid{ss}{2}{selectTrial(1)} = events_in_field_final_filtered{ss}{selectTrial(1)}(mid_field_idx{ss});
    events_in_field_global_centroid{ss}{2}{selectTrial(2)} = events_in_field_final_filtered{ss}{selectTrial(2)}(mid_field_idx{ss});
    
    events_in_field_global_centroid{ss}{3}{selectTrial(1)} = events_in_field_final_filtered{ss}{selectTrial(1)}(far_field_idx{ss});
    events_in_field_global_centroid{ss}{3}{selectTrial(2)} = events_in_field_final_filtered{ss}{selectTrial(2)}(far_field_idx{ss});
    
end

%% Do Mann Whitney U test for AUC of events in nearby place field (rate remapping ROIs) 

for ss=sessionSelect
    %extract AUC values for only filtered neurons
    event_AUC_filtered{ss}.A = event_AUC{ss}.A(near_idx_centroid{ss});
    event_AUC_filtered{ss}.B = event_AUC{ss}.B(near_idx_centroid{ss});
    
    %extract AUC values for in field events
    for rr =1:size(event_AUC_filtered{ss}.A,2)
        %for A trials
        event_AUC_event_filtered{ss}.A{rr} = event_AUC_filtered{ss}.A{rr}(events_in_field_near_centroid{ss}{selectTrial(1)}{rr});
        %for B trials
        event_AUC_event_filtered{ss}.B{rr} = event_AUC_filtered{ss}.B{rr}(events_in_field_near_centroid{ss}{selectTrial(2)}{rr});
    end
    
    %for each ROI, run Mann-Whitney U comparing in field AUC values
    for rr =1:size(event_AUC_filtered{ss}.A,2)
        %return p value for each comparison
        switch options.AUC_test
            case 'ranksum'
                [p_val{ss}(rr),~] = ranksum(event_AUC_event_filtered{ss}.A{rr},event_AUC_event_filtered{ss}.B{rr});
            case 'ks'
                [~,p_val{ss}(rr),~] = kstest2(event_AUC_event_filtered{ss}.A{rr},event_AUC_event_filtered{ss}.B{rr});
        end
    end
end

for ss=sessionSelect
    %logical to select out ROIs from final filter
    rate_remap_pos_log{ss} = p_val{ss} < p_sig;
    common_log{ss} = p_val{ss} >= p_sig;
    
    %rate remapping ROIs + filtered out params
    %select ROI indices (from original)
    rate_remapping_ROI{ss} = near_idx_centroid{ss}(rate_remap_pos_log{ss});
    common_ROI{ss} = near_idx_centroid{ss}(common_log{ss});
    
    %place field width positions
    %rate remapping neurons
    placeField_rate_remap{ss}{selectTrial(1)} = placeField_near_centroid{ss}{selectTrial(1)}(rate_remap_pos_log{ss});
    placeField_rate_remap{ss}{selectTrial(2)} = placeField_near_centroid{ss}{selectTrial(2)}(rate_remap_pos_log{ss});
    
    %common neurons
    placeField_common{ss}{selectTrial(1)} = placeField_near_centroid{ss}{selectTrial(1)}(common_log{ss});
    placeField_common{ss}{selectTrial(2)} = placeField_near_centroid{ss}{selectTrial(2)}(common_log{ss});
    
    %event position (normalized)
    %rate remapping
    event_pos_inField_rate_remap{ss}{selectTrial(1)} = event_pos_inField_near_centroid{ss}{selectTrial(1)}(rate_remap_pos_log{ss});
    event_pos_inField_rate_remap{ss}{selectTrial(2)} = event_pos_inField_near_centroid{ss}{selectTrial(2)}(rate_remap_pos_log{ss});
    
    %common
    event_pos_inField_common{ss}{selectTrial(1)} = event_pos_inField_near_centroid{ss}{selectTrial(1)}(common_log{ss});
    event_pos_inField_common{ss}{selectTrial(2)} = event_pos_inField_near_centroid{ss}{selectTrial(2)}(common_log{ss});
end


%% Partial remapping ROIs

for ss=sessionSelect
    %remove other i'd'd indices from this set
    global_overlap_near_log{ss} = ismember(partial_remap_idx_start{ss},global_remap_ROI{ss}{2});
    global_overlap_far_log{ss} = ismember(partial_remap_idx_start{ss},global_remap_ROI{ss}{3});
    rate_overlap_log{ss} = ismember(partial_remap_idx_start{ss},rate_remapping_ROI{ss});
    common_overlap_log{ss} = ismember(partial_remap_idx_start{ss},common_ROI{ss});
    
    %combine in matrix and do or
    %check for overlaps with partial remapping neurons
    %should be 0 idxs removed b/c only use neurons with single fields
    overlap_mat{ss} = [global_overlap_near_log{ss}; global_overlap_far_log{ss}; rate_overlap_log{ss}; common_overlap_log{ss}];
    overlap_log{ss} = sum(overlap_mat{ss},1);
    
    %make shared logical of all overlapping and remove from
    %partial_remap_idx_start_ROI
    remove_previous_ROI_log{ss} = overlap_log{ss};
    
    %filter out previously categorized neurons
    partial_idx_previous_removed{ss} = partial_remap_idx_start{ss}(~remove_previous_ROI_log{ss});
end

%% Partial ROI filtering here

%centroid diff - common (less than 10 cm) and one more than 30
%make sure than animal was running both fields on either trials
%regardless of remap in each zone
%runs through every session inside of the function
[partial_remap_filtered,partial_field_idx] = filter_partial_remappers_multi_ses(partial_remap_idx_start,cent_diff,select_fields,Place_cell,edges,pf_vector,...
                            Behavior_full,Behavior_split,event_norm_pos_run, event_lap_idx,options);


%% Create mixed category - all neurons that do not fit any category by are tuned in both trials by either criteria

%all neurons tuned in both trial by either category (filtered by 1PF and
%min 5 events in field)

for ss=sessionSelect
    %convert to logical
    tuned_AB_si_ts_filt_log{ss} = logical(zeros(1,size(tuned_AB_si_ts{ss},2)));
    tuned_AB_si_ts_filt_log{ss}(tuned_AB_si_ts_filt_idx{ss}) = 1;
    
    %convert idxs of all selected neurons to logical vectors equal to selected
    %ROIs id'd for the session - create common logical matrix
    blank_logical{ss} = zeros(6,size(tuned_AB_si_ts{ss},2));
    %Assign true to eahc row in the following order
    %1 - common; 2 - global near; 3 -global far
    %4 - rate; 5 - partial; 6 - mixed
    blank_logical{ss}(1,common_ROI{ss}) = 1;
    blank_logical{ss}(2,global_remap_ROI{ss}{2}) = 1;
    blank_logical{ss}(3,global_remap_ROI{ss}{3}) = 1;
    blank_logical{ss}(4,rate_remapping_ROI{ss}) = 1;
    blank_logical{ss}(5,partial_remap_filtered{ss}) = 1;
    %sum the first 5 logicals - should not sum to more than 1
    %create mixed/unclassified category for neurons that do not match the
    %classification above
    blank_logical{ss}(6,:) = (~logical(sum(blank_logical{ss}(1:5,:))) & tuned_AB_si_ts_filt_log{ss});
    %get indices of mixed ROIs
    mixed_ROI{ss} = find(blank_logical{ss}(6,:) == 1);
    
    %Check if each cateory split count adds up to whole
    if(sum(sum(blank_logical{ss},2)) == length(find(tuned_AB_si_ts_filt_log{ss} == 1)))
        disp('Category split adds up to whole of A&B neurons.')
        disp(sum(blank_logical{ss},2))
    else
        disp('Category split does NOT add up to whole of A&B neurons!')
    end
end


%% Export indices of neuron in each category as in struct
%SAVE this such that each session is saved to save each
for ss=sessionSelect
    remapping_ROIs{ss}.global_near = global_remap_ROI{ss}{2};
    remapping_ROIs{ss}.global_far = global_remap_ROI{ss}{3};
    
    remapping_ROIs{ss}.rate = rate_remapping_ROI{ss};
    remapping_ROIs{ss}.common = common_ROI{ss};

    remapping_ROIs{ss}.partial = partial_remap_filtered{ss};
    remapping_ROIs{ss}.mixed = mixed_ROI{ss};
end



%% Plot as shaded area to verify correct id of place field onto normalized
if 0
    %rate remapping
    if options.dispFigure ==1
        %plot normalized position
        %show A selective first
        figure('Position', [1930 130 1890 420])
        hold on
        title('A-selective filtered')
        for rr=562%:size(rate_remapping_ROI,2)%1:size(Aonly_notSIb_idx,2)
            %ROI = rr%AandB_tuned_idx(rr);
            ROI = 562 %rate_remapping_ROI(rr); %Aonly_notSIb_idx(rr);
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
            for ii=1:size(lap_idxs.A,1)
                plot(Imaging_split{1}{1}.time_restricted(lap_idxs.A(ii,1):lap_idxs.A(ii,2))/60,...
                    Behavior_split{1}{1}.resampled.position_norm(lap_idxs.A(ii,1):lap_idxs.A(ii,2)),...
                    'Color',[0 0 1 0.6],'LineWidth',1.5)
            end
            %B laps
            for ii=1:size(lap_idxs.B,1)
                plot(Imaging_split{1}{2}.time_restricted(lap_idxs.B(ii,1):lap_idxs.B(ii,2))/60,...
                    Behavior_split{1}{2}.resampled.position_norm(lap_idxs.B(ii,1):lap_idxs.B(ii,2)),...
                    'Color',[1 0 0 0.6],'LineWidth',1.5)
            end
            %overlay significant calcium run events
            %A
            scatter(event_norm_time.A{ROI},event_norm_pos_run.A{ROI},[],[0 0 1],'*')
            %B
            scatter(event_norm_time.B{ROI},event_norm_pos_run.B{ROI},[],[1 0 0],'*')
            
            %plot horz lines signifying start and end of place field
            %start
            lineS = refline(0,placeField_rate_remap{1}{rr}(1))
            lineS.Color = 'g';
            %end
            lineE = refline(0,placeField_rate_remap{1}{rr}(2))
            lineE.Color = 'g';
            
            %pause
            %clf
        end
        
        
    end
    
    %global
    if options.dispFigure ==1
        %plot normalized position
        %show A selective first
        figure('Position', [1930 130 1890 420])
        hold on
        title('Global remapping neurons')
        for rr=1:size(global_remap_ROI,2)%1:size(Aonly_notSIb_idx,2)
            %ROI = rr%AandB_tuned_idx(rr);
            ROI = global_remap_ROI(rr); %Aonly_notSIb_idx(rr);
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
            for ii=1:size(lap_idxs.A,1)
                plot(Imaging_split{1}{1}.time_restricted(lap_idxs.A(ii,1):lap_idxs.A(ii,2))/60,...
                    Behavior_split{1}{1}.resampled.position_norm(lap_idxs.A(ii,1):lap_idxs.A(ii,2)),...
                    'Color',[0 0 1 0.6],'LineWidth',1.5)
            end
            %B laps
            for ii=1:size(lap_idxs.B,1)
                plot(Imaging_split{1}{2}.time_restricted(lap_idxs.B(ii,1):lap_idxs.B(ii,2))/60,...
                    Behavior_split{1}{2}.resampled.position_norm(lap_idxs.B(ii,1):lap_idxs.B(ii,2)),...
                    'Color',[1 0 0 0.6],'LineWidth',1.5)
            end
            %overlay significant calcium run events
            %A
            scatter(event_norm_time.A{ROI},event_norm_pos_run.A{ROI},[],[0 0 1],'*')
            %B
            scatter(event_norm_time.B{ROI},event_norm_pos_run.B{ROI},[],[1 0 0],'*')
            
            %plot horz lines signifying start and end of place field
            %start
            lineS = refline(0,placeField_global_centroid{1}{rr}(1))
            lineS.Color = 'g';
            %end
            lineE = refline(0,placeField_global_centroid{1}{rr}(2))
            lineE.Color = 'g';
            
            pause
            clf
        end
        
    end
    
    %partial
    if options.dispFigure ==1
        %plot normalized position
        %show A selective first
        figure('Position', [1930 130 1890 420])
        hold on
        title('Partial remapping neurons')
        for rr=1:size(partial_idx_previous_removed,2)%1:size(Aonly_notSIb_idx,2)
            %ROI = rr%AandB_tuned_idx(rr);
            ROI = partial_idx_previous_removed(rr); %Aonly_notSIb_idx(rr);
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
            for ii=1:size(lap_idxs.A,1)
                plot(Imaging_split{1}{1}.time_restricted(lap_idxs.A(ii,1):lap_idxs.A(ii,2))/60,...
                    Behavior_split{1}{1}.resampled.position_norm(lap_idxs.A(ii,1):lap_idxs.A(ii,2)),...
                    'Color',[0 0 1 0.6],'LineWidth',1.5)
            end
            %B laps
            for ii=1:size(lap_idxs.B,1)
                plot(Imaging_split{1}{2}.time_restricted(lap_idxs.B(ii,1):lap_idxs.B(ii,2))/60,...
                    Behavior_split{1}{2}.resampled.position_norm(lap_idxs.B(ii,1):lap_idxs.B(ii,2)),...
                    'Color',[1 0 0 0.6],'LineWidth',1.5)
            end
            %overlay significant calcium run events
            %A
            scatter(event_norm_time.A{ROI},event_norm_pos_run.A{ROI},[],[0 0 1],'*')
            %B
            scatter(event_norm_time.B{ROI},event_norm_pos_run.B{ROI},[],[1 0 0],'*')
            
            %plot horz lines signifying start and end of place field
            %start
            %lineS = refline(0,placeField_global_centroid{1}{rr}(1))
            %lineS.Color = 'g';
            %end
            %lineE = refline(0,placeField_global_centroid{1}{rr}(2))
            %lineE.Color = 'g';
            
            pause
            clf
        end
        
    end
end

%% Export task-selective ROIs in struct

%export indices for task selective neurons
% task_selective_ROIs.A.idx = final_filtered_ROI.A;
% task_selective_ROIs.B.idx = final_filtered_ROI.B;
% 
% %field margin for task selective neurons
% task_selective_ROIs.A.field_margin = placeField_final{1};
% task_selective_ROIs.B.field_margin = placeField_final{2};
% 
% %event associated with selected place field
% task_selective_ROIs.A.fieldEvents = event_pos_inField_final{1};
% task_selective_ROIs.B.fieldEvents = event_pos_inField_final{2};




end

