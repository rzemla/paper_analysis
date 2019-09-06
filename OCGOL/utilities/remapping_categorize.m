function [remapping_ROIs] = remapping_categorize(cent_diff, tuned_logical, pf_vector_max, session_vars,max_transient_peak,pf_count_filtered_log, pf_count_filtered, options)
%split mutually tuned neurons by remapping category: 
%common (less than certain centroid difference between max
%tuned_log = tunedLogical.ts.AandB_tuned;

%% Set parameters

%number of degrees of centroid difference
deg_thres = options.deg_thres;
%degree range for splitting the global remapping neurons
deg_ranges = options.deg_ranges;

%% Get ROI indices of A selective, B selective, A&B selective neurons

%find A&B tuned by either criterion (input option below)
%unclassified neurons get moved into mixed category
tuned_A_si_ts = tuned_logical.si.Atuned  | tuned_logical.ts.Atuned;
tuned_B_si_ts = tuned_logical.si.Btuned  | tuned_logical.ts.Btuned;
tuned_AB_si_ts = tuned_A_si_ts & tuned_B_si_ts;

%rate map, common, and global
%only TS tuned neurons
rate_remap_idx_start = find(tuned_logical.ts.AandB_tuned ==1);
%rate_remap_idx_start = find(tuned_AB_si_ts == 1); 

%partial remapping - both A and B tuned to SI
%partial_remap_idx_start = find(tuned_logical.si.AandB_tuned == 1);
partial_remap_idx_start = find(tuned_AB_si_ts == 1);

%choose if SI or TS tuned
switch options.tuning_criterion
    case 'si' %spatial information
        AandB_tuned_idx = find(tuned_logical.si.AandB_tuned == 1);
        Aonly_tuned_idx = find(tuned_logical.si.onlyA_tuned == 1);
        Bonly_tuned_idx = find(tuned_logical.si.onlyB_tuned == 1);
        
    case 'ts' %tuning specificity
        %both A and B tuned by TS
        AandB_tuned_idx = find(tuned_logical.ts.AandB_tuned ==1);
        %only A tuned by TS
        Aonly_tuned_idx = find(tuned_logical.ts.onlyA_tuned == 1);
        %only B tuned by tS
        Bonly_tuned_idx = find(tuned_logical.ts.onlyB_tuned == 1);
        
        %all A tuned by SI
        A_tuned_si_idx = find(tuned_logical.si.Atuned == 1);
        %all B tuned by SI
        B_tuned_si_idx = find(tuned_logical.si.Btuned == 1);
        
        %only A tuned by TS, but also not SI B tuned
        Aonly_notSIb_idx = setdiff(Aonly_tuned_idx,B_tuned_si_idx);
        
        %only B tuned by TS, but also not SI A tuned
        Bonly_notSIa_idx = setdiff(Bonly_tuned_idx,A_tuned_si_idx);

end

%% Filter make sure rate/global/common remapping category has single PF in each and meets min event if field criteriym
%use place field calculation and logical from place properties calculation
single_pf_idxs = find(sum((pf_count_filtered == 1),1) ==2);
%make copy of original rate remap vector
rate_remap_idx_orig = rate_remap_idx_start;
%extract only the indices
rate_remap_idx_start = intersect(rate_remap_idx_orig,single_pf_idxs);

%% Parse starting partial remap neurons (2PF vs 1 PF)
%ROIs with siginificant 2PF vs. 1PFs
single_double_pf_idxs = find(sum((pf_count_filtered == 1) + 2*(pf_count_filtered == 2),1) ==3);
%make copy of original rate remap vector
partial_remap_idx_orig = partial_remap_idx_start;
%extract only the indices that have 2PF vs. 1 PF place fields
partial_remap_idx_start = intersect(partial_remap_idx_orig,single_double_pf_idxs);

%% Truncate the centroid difference struct to only include the selected ROIs
cent_diff_AandB.angle_diff = cent_diff.angle_diff(rate_remap_idx_start);
cent_diff_AandB.max_bin = cent_diff.angle_diff(:,rate_remap_idx_start);

%% Define/load variables for each session

%for each session
for ii = 1:size(session_vars,2)
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
for ss=1:size(session_vars,2)
    %for each lap
    for ii=1:size(run_intervals{ss},2)
        events{ss}{ii} = events_full{ss}{ii}(logical(run_intervals{ss}{ii}),:);
    end
end

%% Get lap indices for each lap in all A or B trials

%only correct

%get unique lap indices
lapA_idxs = unique(Behavior_split{1}{1}.resampled.lapNb);
lapB_idxs = unique(Behavior_split{1}{2}.resampled.lapNb);

%get lap start and end indices for all A or B trials
%all A
for ll=1:size(lapA_idxs,1)
    lap_idxs.A(ll,1) = find(Behavior_split{1}{1}.resampled.lapNb == lapA_idxs(ll),1,'first');
    lap_idxs.A(ll,2) = find(Behavior_split{1}{1}.resampled.lapNb == lapA_idxs(ll),1,'last');
end

%all B
for ll=1:size(lapB_idxs,1)
    lap_idxs.B(ll,1) = find(Behavior_split{1}{2}.resampled.lapNb == lapB_idxs(ll),1,'first');
    lap_idxs.B(ll,2) = find(Behavior_split{1}{2}.resampled.lapNb == lapB_idxs(ll),1,'last');
end

%% Event onsets in run interval
%only correct trials (1,2) - 
%TODO: MODIFY THIS TO TAKE IN THE CHOSEN TRIALS TO
%USE VS. FIXED AS IS NOW
%for each ROI
for rr=1:size(events{ss}{1},2)
    %time of significant run events in A
    event_norm_time.A{rr} = Imaging_split{1}{1}.time_restricted(find(Event_split{1}{1}.Run.run_onset_binary(:,rr) == 1))/60;
    %normalizesd position of significant run events in A
    event_norm_pos_run.A{rr} = Behavior_split{1}{1}.resampled.position_norm(find(Event_split{1}{1}.Run.run_onset_binary(:,rr) == 1));
    %lap assignment for A
    event_lap_idx.A{rr} = Behavior_split{1}{1}.resampled.lapNb(logical(Event_split{1, 1}{1, 1}.Run.run_onset_binary(:,rr)));
    %get respective AUC values
    event_AUC.A{rr} = Event_split{1}{1}.Run.properties.AUC{rr};
    
    %time of significant run events in B
    event_norm_time.B{rr} = Imaging_split{1}{2}.time_restricted(find(Event_split{1}{2}.Run.run_onset_binary(:,rr) == 1))/60;
    %normalizesd position of significant run events in B
    event_norm_pos_run.B{rr} = Behavior_split{1}{2}.resampled.position_norm(find(Event_split{1}{2}.Run.run_onset_binary(:,rr) == 1));
    %lap assignment for B
    event_lap_idx.B{rr} = Behavior_split{1}{2}.resampled.lapNb(logical(Event_split{1}{2}.Run.run_onset_binary(:,rr)));
    %get respective AUC values
    event_AUC.B{rr} = Event_split{1}{2}.Run.properties.AUC{rr};
    
end


%% Plot the run epochs, corresponding position bin edges of place field

%overlay area of max place field bin

%TODO: MODIFY THIS TO TAKE IN THE CHOSEN TRIALS TO
%USE VS. FIXED AS IS NOW

%spatial bin assignment for each run-epoch frame (100 bins) 
%A trials
session_vars{1}.Place_cell{1}.Bin{8};
%B trials
session_vars{1}.Place_cell{2}.Bin{8};
%get the run epoch binaries for each set of trials - A session - 5339 run
%frames
%plot the run epochs as patches
session_vars{1}.Behavior_split{1}.run_ones;
%get the run epoch binaries for each set of trials - B session - 4547 run
%frames
session_vars{1}.Behavior_split{2}.run_ones; 
%run events - binary onset across entire run/no run interval
Event_split{1}{1}.Run.run_onset_ones;
%get bins
%get edges for corresponding bins!! - find place in spatial info where 
%correct A trials
run_position_norm{1} = Behavior_split{1}{1}.resampled.run_position_norm;
%correct B trials
run_position_norm{2} = Behavior_split{1}{2}.resampled.run_position_norm;
%Bin running position in 100 bins and get edges for each set of laps:
%for each number of bins, bin the normalized position during run epochs
%for correct A trials
[count_bin{1},edges{1},bin{1}] = histcounts(run_position_norm{1}, 100);
%for correct B trials
[count_bin{2},edges{2},bin{2}] = histcounts(run_position_norm{2}, 100);


%% Remove ROIs idx's without a id'd place field _ THIS SHOULD BE REDUNDANT GIVEN FILTERING DONE ABOVE

%edges of all identified common TS neurons
%correct A
placeFieldEdges{1} = Place_cell{1}{1}.placeField.edge(rate_remap_idx_start);
%correct B
placeFieldEdges{2} = Place_cell{1}{2}.placeField.edge(rate_remap_idx_start);

%in A trials
idx_wo_placeFields{1} = rate_remap_idx_start(find(cellfun(@isempty,placeFieldEdges{1}) ==1));
%in B trials
idx_wo_placeFields{2} = rate_remap_idx_start(find(cellfun(@isempty,placeFieldEdges{2}) ==1));

%merge empty field indices into 1 vector for removal from list
rm_idx_no_pf = unique(cell2mat(idx_wo_placeFields));

%remove ROIs for A and B that do not have place fields
if ~isempty(rm_idx_no_pf)
    rate_remap_idx_filtered = setdiff(rate_remap_idx_start,rm_idx_no_pf);

else %keep the previous ROIs (copy only)
    rate_remap_idx_filtered = rate_remap_idx_start;
end

%filter out centroid difference data based on no-id'd neuron tuned by TS
%find indices of removed neurons from list
select_pf_filtered_log = ismember(rate_remap_idx_start,rate_remap_idx_filtered);

%update centroid input data with pf filtered idxs
cent_diff_AandB_pf_filt.angle_diff = cent_diff_AandB.angle_diff(select_pf_filtered_log);
cent_diff_AandB_pf_filt.max_bin = cent_diff_AandB.max_bin(:,select_pf_filtered_log);

%get the ROI Idxs associated with place filtered ROIS
remapping_pf_filtered = rate_remap_idx_start(select_pf_filtered_log);

%% Determine the equivalent normalized position range of the max place field

%extract the max place field index for A and B trial using filtered A tuned
%and B tuned ROIs
max_field_idx{1} = max_transient_peak{1}{1}(remapping_pf_filtered);
max_field_idx{2} = max_transient_peak{1}{2}(remapping_pf_filtered);

%get the edges of the max transient place field for each set of idxs
%edges of all identified
%correct A
placeField_filtered{1} = Place_cell{1}{1}.placeField.edge(remapping_pf_filtered);
%correct B
placeField_filtered{2} = Place_cell{1}{2}.placeField.edge(remapping_pf_filtered);

%check which field has more than 1 field and select edges of the one with
%higher transient rate

%for each trial (A and B )
for tt=1:2
    for rr=1:size(placeField_filtered{tt},2)
        if size(placeField_filtered{tt}{rr},1) > 1
            placeField_filtered_max{tt}{rr} = placeField_filtered{tt}{rr}(max_field_idx{tt}(rr),:);
        else
            placeField_filtered_max{tt}{rr} = placeField_filtered{tt}{rr};
        end
    end
end

%convert edges from relevant place field to normalized postion edges
for tt=1:2
    for rr=1:size(placeField_filtered_max{tt},2)
        %start position of PF
        placeField_filtered_max_posnorm{tt}{rr}(1) = edges{1}(placeField_filtered_max{tt}{rr}(1));
        %end position of PF
        placeField_filtered_max_posnorm{tt}{rr}(2) = edges{1}(placeField_filtered_max{tt}{rr}(2)+1)-0.01;
    end
end


%% Get normalized position distance converstion factor

%get median lap length based on the registered length of each lap
median_track_len = median(Behavior_full{1}.position_lap(:,2));

%conversion factor (norm_pos/cm length) - 100 bins
norm_conv_factor = median_track_len/100;
%median_track_len/1;

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

%find events occuring within max place field for each ROI
for tt=1:2
    for rr=1:size(placeField_filtered_max_posnorm{tt},2)
        %get idxs of events with max place field
        if tt == 1 %correct A trials
        %OLDER VERSIONS - REPLACED BY CODE BELOW
            %events_in_field{tt}{rr} = find(event_norm_pos_run.A{remapping_pf_filtered(rr)} >= placeField_filtered_max_posnorm{tt}{rr}(1) & ...
             %   event_norm_pos_run.A{remapping_pf_filtered(rr)} <= placeField_filtered_max_posnorm{tt}{rr}(2));
            
            %check of the normalized field position edges cross the start of the track
            %(first edge) - i.e. if first edge is greater than
            %the second edge
            if placeField_filtered_max_posnorm{tt}{rr}(1) < placeField_filtered_max_posnorm{tt}{rr}(2)
                events_in_field{tt}{rr} = find(event_norm_pos_run.A{remapping_pf_filtered(rr)} >= placeField_filtered_max_posnorm{tt}{rr}(1) & ...
                    event_norm_pos_run.A{remapping_pf_filtered(rr)} <= placeField_filtered_max_posnorm{tt}{rr}(2));
            else %do two event finds and merge (and sort) into one set of indices
                events_in_field_temp_1 = find(event_norm_pos_run.A{remapping_pf_filtered(rr)} >= 0 & ...
                    event_norm_pos_run.A{remapping_pf_filtered(rr)} <= placeField_filtered_max_posnorm{tt}{rr}(2));
                events_in_field_temp_2 = find(event_norm_pos_run.A{remapping_pf_filtered(rr)} >= placeField_filtered_max_posnorm{tt}{rr}(2) & ...
                    event_norm_pos_run.A{remapping_pf_filtered(rr)} <= 1 );
                
                %merge and sort indices here
                events_in_field{tt}{rr} = sort([events_in_field_temp_1; events_in_field_temp_2]);
            end


            %register the corresponding lap of in-field filtered event
            event_in_field_laps{tt}{rr} = event_lap_idx.A{remapping_pf_filtered(rr)}(events_in_field{tt}{rr});
            %get number of unique events (those occuring on each lap)
            event_in_field_nb{tt}{rr} = size(unique(event_in_field_laps{tt}{rr}),1);
            %get position of in-field events
            events_in_field_pos{tt}{rr} = event_norm_pos_run.A{remapping_pf_filtered(rr)}(events_in_field{tt}{rr});
            
        elseif tt == 2 %correct B trials
            %events_in_field{tt}{rr} = find(event_norm_pos_run.B{remapping_pf_filtered(rr)} >= placeField_filtered_max_posnorm{tt}{rr}(1) & ...
            %    event_norm_pos_run.B{remapping_pf_filtered(rr)} <= placeField_filtered_max_posnorm{tt}{rr}(2));

            %check of the normalized field position edges cross the start of the track
            %(first edge) - i.e. if first edge is greater than
            %the second edge
            if placeField_filtered_max_posnorm{tt}{rr}(1) < placeField_filtered_max_posnorm{tt}{rr}(2)
                events_in_field{tt}{rr} = find(event_norm_pos_run.B{remapping_pf_filtered(rr)} >= placeField_filtered_max_posnorm{tt}{rr}(1) & ...
                    event_norm_pos_run.B{remapping_pf_filtered(rr)} <= placeField_filtered_max_posnorm{tt}{rr}(2));
            else %do two event finds and merge (and sort) into one set of indices
                events_in_field_temp_1 = find(event_norm_pos_run.B{remapping_pf_filtered(rr)} >= 0 & ...
                    event_norm_pos_run.B{remapping_pf_filtered(rr)} <= placeField_filtered_max_posnorm{tt}{rr}(2));
                events_in_field_temp_2 = find(event_norm_pos_run.B{remapping_pf_filtered(rr)} >= placeField_filtered_max_posnorm{tt}{rr}(2) & ...
                    event_norm_pos_run.B{remapping_pf_filtered(rr)} <= 1 );
                
                %merge and sort indices here
                events_in_field{tt}{rr} = sort([events_in_field_temp_1; events_in_field_temp_2]);
            end

            %register the corresponding lap of in-field filtered event
            event_in_field_laps{tt}{rr} = event_lap_idx.B{remapping_pf_filtered(rr)}(events_in_field{tt}{rr});
            %get number of unique events (those occuring on each lap)
            event_in_field_nb{tt}{rr} = size(unique(event_in_field_laps{tt}{rr}),1);
            %get position of in-field events
            events_in_field_pos{tt}{rr} = event_norm_pos_run.B{remapping_pf_filtered(rr)}(events_in_field{tt}{rr});
        end
    end
end

%check which task-selective ROIs have less than 5 events
%correct A trials
event_thres_exclude_log.A  = cell2mat(event_in_field_nb{1}) < 5;
%correct B trials
event_thres_exclude_log.B  = cell2mat(event_in_field_nb{2}) < 5;

%merge the excludes neurosn idx
event_thres_excludeAB_log =  event_thres_exclude_log.A | event_thres_exclude_log.B;

%update indices with event - update absolute idxs
ROI_field_filtered_event.A = remapping_pf_filtered(~event_thres_excludeAB_log);
ROI_field_filtered_event.B = remapping_pf_filtered(~event_thres_excludeAB_log);

%select the event indices for event/place filtered ROIs
events_in_field_pf_filtered{1} = events_in_field{1}(~event_thres_excludeAB_log);
events_in_field_pf_filtered{2} = events_in_field{2}(~event_thres_excludeAB_log);

%update assn place fields
placeField_eventFilt{1} = placeField_filtered_max_posnorm{1}(~event_thres_excludeAB_log);
placeField_eventFilt{2} = placeField_filtered_max_posnorm{2}(~event_thres_excludeAB_log);

%update event position (normalized)
event_pos_inField{1} = events_in_field_pos{1}(~event_thres_excludeAB_log);
event_pos_inField{2} = events_in_field_pos{2}(~event_thres_excludeAB_log);

%update centroid input data with pf filtered idxs
cent_diff_AandB_eventCount.angle_diff = cent_diff_AandB_pf_filt.angle_diff(~event_thres_excludeAB_log);
cent_diff_AandB_eventCount.max_bin = cent_diff_AandB_pf_filt.max_bin(:,~event_thres_excludeAB_log);

%% Make sure the animal was in a run epoch in the min/max range of space on opposing laps (al least 80%) of space on at least 6 laps

% Check that in running epoch within 3 bins to the left or right of each
ROI_field_filtered_event.A
ROI_field_filtered_event.B

%take the median position of the events in field and min/max position
for tt=1:2 %for correct A and B trials
    for rr=1:size(event_pos_inField{tt},2)
        med_pos_event{tt}(rr) = median(event_pos_inField{tt}{rr}); 
        %into one matrix min and max of each event
        min_max_pos_event{tt}(rr,1) = min(event_pos_inField{tt}{rr}); 
        min_max_pos_event{tt}(rr,2) = max(event_pos_inField{tt}{rr}); 
    end
end

%get correct A and B laps idx
corr_lap_idx{1} = unique(Behavior_split{1}{1}.resampled.lapNb);
corr_lap_idx{2} = unique(Behavior_split{1}{2}.resampled.lapNb);

%get indices across all laps of the ranges
for tt=1:2 %for correct A and B trials
    for rr=1:size(event_pos_inField{tt},2)
        %get indices that match the position range (ALL LAPS)
        pos_range_indices{tt}{rr} = find( Behavior_full{1}.resampled.normalizedposition >= min_max_pos_event{tt}(rr,1) & ...
            Behavior_full{1}.resampled.normalizedposition <= min_max_pos_event{tt}(rr,2));
        %get lap idx of corresponding idxs
        lap_idx_range{tt}{rr} = Behavior_full{1}.resampled.lapNb(pos_range_indices{tt}{rr});
    end
end

%for correct A or B trials
for tt=1:2
    %for each ROI in correct A trials 
    for rr=1:size(lap_idx_range{tt},2)
        %extact logical with only laps correponding to opposing trial laps- dependent on B parameter
        if tt == 1 %if looking on run status in B trials for correct A trials
            lap_opposed_idx{tt}{rr} = ismember(lap_idx_range{tt}{rr},corr_lap_idx{2});
        elseif tt == 2 %if looking on run status in A trials for correct B trials
            lap_opposed_idx{tt}{rr} = ismember(lap_idx_range{tt}{rr},corr_lap_idx{1});
        end
        %get the lap number associated with each frame in the opposing trials
        lap_label_opposed{tt}{rr} = lap_idx_range{tt}{rr}(lap_opposed_idx{tt}{rr});
        %the binary indicating in animal in run epoch with that range
        lap_runEpoch_opposed{tt}{rr} = Behavior_full{1}.run_ones(pos_range_indices{tt}{rr}(lap_opposed_idx{tt}{rr}));
        %extract the associated positions
        lap_pos_opposed{tt}{rr} = Behavior_full{1}.resampled.normalizedposition(pos_range_indices{tt}{rr}(lap_opposed_idx{tt}{rr}));
    end
end

%split into individual laps for A events, look in B laps
%for correct A or B trials
for tt=1:2
    %for each ROI in correct A trials
    for rr=1:size(lap_idx_range{tt},2)
        %for each opposing lap
        for ll=1:size(corr_lap_idx{2},1)
            % - depdendent on B parameter
            if tt == 1 %if looking on run status in B trials for correct A trials
                split_lap_idxs{tt}{rr}{ll} = find(lap_label_opposed{tt}{rr} == corr_lap_idx{2}(ll));
            elseif tt == 2
                split_lap_idxs{tt}{rr}{ll} = find(lap_label_opposed{tt}{rr} == corr_lap_idx{1}(ll));
            end
            
            split_lap_pos{tt}{rr}{ll} = lap_pos_opposed{tt}{rr}(split_lap_idxs{tt}{rr}{ll});
            split_lap_runEpoch{tt}{rr}{ll} = lap_runEpoch_opposed{tt}{rr}(split_lap_idxs{tt}{rr}{ll});
        end
    end
end

%for correct A or B trials
for tt=1:2
    %for each ROI in correct A trials
    for rr=1:size(lap_idx_range{tt},2)
        %for each lap
        for ll=1:size(split_lap_pos{tt}{rr},2)
            unique_pos{tt}{rr}{ll} = unique(split_lap_pos{tt}{rr}{ll});
            %for each unique position, check if entirety in run epoch
            for pos_idx=1:size(unique_pos{tt}{rr}{ll},1)
                %get idx associated with given unique pos in the lap
                pos_idxs_each{tt}{rr}{ll}{pos_idx} = find(split_lap_pos{tt}{rr}{ll} == unique_pos{tt}{rr}{ll}(pos_idx));
                run_Epoch_each{tt}{rr}{ll}{pos_idx} = split_lap_runEpoch{tt}{rr}{ll}(pos_idxs_each{tt}{rr}{ll}{pos_idx});
                %check which position did animal spend all of it in run state
                
            end
            %calculate logical of  run at each position and get fraction in run
            %state for each lap
            frac_run_lap{tt}{rr}(ll) = sum(cellfun(@prod,run_Epoch_each{tt}{rr}{ll}))/size(run_Epoch_each{tt}{rr}{ll},2);
        end
    end
end

%for all laps check how above 80% of space in run epoch and check if this
%occurs in at least 6 laps
%for correct A or B trials
for tt=1:2
    %for each ROI in correct A trials
    for rr=1:size(lap_idx_range{tt},2)
        if sum(frac_run_lap{tt}{rr} >= 0.8) >= 6
            %generate logical with 1's include ROI, and 0 exclude ROI
            run_epoch_filt_include_log{tt}(rr) = 1;
        else
            run_epoch_filt_include_log{tt}(rr) = 0;
        end
    end
end

%combine filters for each epoch into 1 filter
run_epoch_filt_both = logical(run_epoch_filt_include_log{1}) | logical(run_epoch_filt_include_log{1});

%apply final filter to ROI indices

%select ROI indices (from original)
final_filtered_ROI.A = ROI_field_filtered_event.A(run_epoch_filt_both);
final_filtered_ROI.B = ROI_field_filtered_event.B(run_epoch_filt_both);

%place field width positions
placeField_final{1} = placeField_eventFilt{1}(run_epoch_filt_both);
placeField_final{2} = placeField_eventFilt{2}(run_epoch_filt_both);

%event position (normalized)
event_pos_inField_final{1} = event_pos_inField{1}(run_epoch_filt_both);
event_pos_inField_final{2} = event_pos_inField{2}(run_epoch_filt_both);

%select the event indices for event/place filtered ROIs
events_in_field_final_filtered{1} = events_in_field_pf_filtered{1}(run_epoch_filt_both);
events_in_field_final_filtered{2} = events_in_field_pf_filtered{2}(run_epoch_filt_both);

%update centroid input data with pf filtered idxs
cent_diff_AandB_runFilt.angle_diff = cent_diff_AandB_eventCount.angle_diff(run_epoch_filt_both);
cent_diff_AandB_runFilt.max_bin = cent_diff_AandB_eventCount.max_bin(:,run_epoch_filt_both);


%% Centroid difference filter (select those with cent diff less than 10cm ~ 18 deg) - point of split between rate and global remap
%centroid difference in degrees (36 deg ~ 20 cm) (27 deg ~ 15 cm)
%deg_thres = 45;

%get the vector idx's of the neurons with near and far centroid differences
near_field_idx = find((cent_diff_AandB_runFilt.angle_diff < deg2rad(deg_ranges(2))) == 1);
%mid_field_idx = find((cent_diff_AandB_runFilt.angle_diff <= deg2rad(deg_thres)) == 1)
mid_field_idx =find(((cent_diff_AandB_runFilt.angle_diff >= deg2rad(deg_ranges(2))) &...
 (cent_diff_AandB_runFilt.angle_diff < deg2rad(deg_ranges(3)))) ==1);
%far centroid split ROIs
far_field_idx = find((cent_diff_AandB_runFilt.angle_diff >= deg2rad(deg_ranges(3))) == 1);



%translate the idxs to absolute ROIs idxs (A and B are the same)
near_idx_centroid = final_filtered_ROI.A(near_field_idx);
%second cateogory
global_remap_ROI{2} = final_filtered_ROI.A(mid_field_idx);
%third category
global_remap_ROI{3} = final_filtered_ROI.A(far_field_idx);

%place field width positions
%potential rate remapping
placeField_near_centroid{1} = placeField_final{1}(near_field_idx);
placeField_near_centroid{2} = placeField_final{2}(near_field_idx);
 %RESUME HERE
%global remapping (mid and far field)
placeField_global_centroid{3}{1} = placeField_final{1}(far_field_idx);
placeField_global_centroid{3}{2} = placeField_final{2}(far_field_idx);

%event position (normalized)
%potential rate remapping
event_pos_inField_near_centroid{1} = event_pos_inField_final{1}(near_field_idx);
event_pos_inField_near_centroid{2} = event_pos_inField_final{2}(near_field_idx);

%global remapping
event_pos_inField_global_centroid{1} = event_pos_inField_final{1}(far_field_idx);
event_pos_inField_global_centroid{2} = event_pos_inField_final{2}(far_field_idx);

%select the event indices for event/place filtered ROIs
events_in_field_near_centroid{1} = events_in_field_final_filtered{1}(near_field_idx);
events_in_field_near_centroid{2} = events_in_field_final_filtered{2}(near_field_idx);

%global remapping
events_in_field_global_centroid{1} = events_in_field_final_filtered{1}(far_field_idx);
events_in_field_global_centroid{2} = events_in_field_final_filtered{2}(far_field_idx);

%% Do Mann Whitney U test for AUC of events in nearby place field (rate remapping ROIs) 

%extract AUC values for only filtered neurons
event_AUC_filtered.A = event_AUC.A(near_idx_centroid);
event_AUC_filtered.B = event_AUC.B(near_idx_centroid);

%extract AUC values for in field events
for rr =1:size(event_AUC_filtered.A,2)
    %for A trials
    event_AUC_event_filtered.A{rr} = event_AUC_filtered.A{rr}(events_in_field_near_centroid{1}{rr});
    %for B trials
    event_AUC_event_filtered.B{rr} = event_AUC_filtered.B{rr}(events_in_field_near_centroid{2}{rr});
end

%for each ROI, run Mann-Whitney U comparing in field AUC values
for rr =1:size(event_AUC_filtered.A,2)
    %return p value for each comparison
    switch options.AUC_test
        case 'ranksum'
            [p_val(rr),~] = ranksum(event_AUC_event_filtered.A{rr},event_AUC_event_filtered.B{rr});
        case 'ks'
            [~,p_val(rr),~] = kstest2(event_AUC_event_filtered.A{rr},event_AUC_event_filtered.B{rr});
    end
end
%logical to select out ROIs from final filter
rate_remap_pos_log = p_val <0.01;
common_log = p_val >= 0.01;

%rate remapping ROIs + filtered out params
%select ROI indices (from original)
rate_remapping_ROI = near_idx_centroid(rate_remap_pos_log);
common_ROI = near_idx_centroid(common_log);

%place field width positions
%rate remapping neurons
placeField_rate_remap{1} = placeField_near_centroid{1}(rate_remap_pos_log);
placeField_rate_remap{2} = placeField_near_centroid{2}(rate_remap_pos_log);

%common neurons
placeField_common{1} = placeField_near_centroid{1}(common_log);
placeField_common{2} = placeField_near_centroid{2}(common_log);

%event position (normalized)
%rate remapping
event_pos_inField_rate_remap{1} = event_pos_inField_near_centroid{1}(rate_remap_pos_log);
event_pos_inField_rate_remap{2} = event_pos_inField_near_centroid{2}(rate_remap_pos_log);

%common
event_pos_inField_common{1} = event_pos_inField_near_centroid{1}(common_log);
event_pos_inField_common{2} = event_pos_inField_near_centroid{2}(common_log);

%select the event indices for event/place filtered ROIs
%events_in_field_final_filtered{1} = events_in_field_pf_filtered{1}(run_epoch_filt_both);
%events_in_field_final_filtered{2} = events_in_field_pf_filtered{2}(run_epoch_filt_both);

%% Partial remapping ROIs

%remove other i'd'd indices from this set
global_overlap_log = ismember(partial_remap_idx_start,global_remap_ROI);
rate_overlap_log = ismember(partial_remap_idx_start,rate_remapping_ROI);
common_overlap_log = ismember(partial_remap_idx_start,common_ROI);

%make shared logical of all overlapping and remove from
%partial_remap_idx_start_ROI
remove_previous_ROI_log =  (common_overlap_log | (global_overlap_log | rate_overlap_log));

%filter out previously categorized neurons
partial_idx_previous_removed = partial_remap_idx_start(~remove_previous_ROI_log);

%% Export indices of neuron in each category as in struct

remapping_ROIs.global = global_remap_ROI;
remapping_ROIs.rate = rate_remapping_ROI;
remapping_ROIs.common = common_ROI;

%% Plot as shaded area to verify correct id of place field onto normalized

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
%% Export task-selective ROIs in struct

%export indices for task selective neurons
task_selective_ROIs.A.idx = final_filtered_ROI.A;
task_selective_ROIs.B.idx = final_filtered_ROI.B;

%field margin for task selective neurons
task_selective_ROIs.A.field_margin = placeField_final{1};
task_selective_ROIs.B.field_margin = placeField_final{2};

%event associated with selected place field
task_selective_ROIs.A.fieldEvents = event_pos_inField_final{1};
task_selective_ROIs.B.fieldEvents = event_pos_inField_final{2};

%% Patch generator for run epochs
%creates x range of patches
if 0
figure;
hold on
xRange = {5:10; 20; 40:42; 50; 60:75; 90:95};
for k1 = 1:size(xRange,1)
    q = xRange{k1};
    %if very brief period
    if length(q) == 1
        q = [q q+0.1];
    end
    %start x end x; end x start x
    qx = [min(q) max(q)  max(q)  min(q)];
    yl = ylim;
    %start y x2; end y x2
    qy = [[1 1]*yl(1) [1 1]*yl(2)];
    %plot the patches
    %green
    %patch(qx, qy, [0 1 0],'EdgeColor', 'none','FaceAlpha', 0.3)
    %red
    patch(qx, qy, [1 0 0],'EdgeColor', 'none','FaceAlpha', 0.3)
end
end

%% SECOND PLOTTER - discard later
    %show B selective second
%     figure('Position', [1930 130 1890 420])
%     hold on
%     title('B-selective filtered')
%     for rr=1:size(rate_remapping_ROI,2)%1:size(Aonly_notSIb_idx,2)
%         %ROI = rr%AandB_tuned_idx(rr);
%         ROI = rate_remapping_ROI(rr); %Aonly_notSIb_idx(rr);
%         hold on
%         title(num2str(ROI))
%         yticks([0 0.5 1])
%         ylabel('Normalized position')
%         xlabel('Time [min]');
%         xticks(0:3:12);
%         %ylim([0 1])
%         set(gca,'FontSize',14)
%         set(gca,'LineWidth',1)
%         %A laps
%         for ii=1:size(lap_idxs.A,1)
%             plot(Imaging_split{1}{1}.time_restricted(lap_idxs.A(ii,1):lap_idxs.A(ii,2))/60,...
%                 Behavior_split{1}{1}.resampled.position_norm(lap_idxs.A(ii,1):lap_idxs.A(ii,2)),...
%                 'Color',[0 0 1 0.6],'LineWidth',1.5)
%         end
%         %B laps
%         for ii=1:size(lap_idxs.B,1)
%             plot(Imaging_split{1}{2}.time_restricted(lap_idxs.B(ii,1):lap_idxs.B(ii,2))/60,...
%                 Behavior_split{1}{2}.resampled.position_norm(lap_idxs.B(ii,1):lap_idxs.B(ii,2)),...
%                 'Color',[1 0 0 0.6],'LineWidth',1.5)
%         end
%         %overlay significant calcium run events
%         %A
%         scatter(event_norm_time.A{ROI},event_norm_pos_run.A{ROI},[],[0 0 1],'*')
%         %B
%         scatter(event_norm_time.B{ROI},event_norm_pos_run.B{ROI},[],[1 0 0],'*')
%         
%         %plot horz lines signifying start and end of place field
%         %start
%         lineS = refline(0,placeField_final{2}{rr}(1))
%         lineS.Color = 'g';
%         %end
%         lineE = refline(0,placeField_final{2}{rr}(2))
%         lineE.Color = 'g';
%         
%         pause
%         clf
%     end

%% OLD PLOTTING CODE
% %plot normalized position
% figure('Position', [1930 130 1890 420])
% for rr=1:size(Bonly_field_filtered,2)%1:size(Aonly_notSIb_idx,2)
%     %ROI = rr%AandB_tuned_idx(rr);
%     ROI = Bonly_field_filtered(rr); %Aonly_notSIb_idx(rr);
%     hold on
%     title(num2str(ROI))
%     yticks([0 0.5 1])
%     ylabel('Normalized position')
%     xlabel('Time [min]');
%     xticks(0:3:12);
%     %ylim([0 1])
%     set(gca,'FontSize',14)
%     set(gca,'LineWidth',1)
%     %A laps
%     for ii=1:size(lap_idxs.A,1)
%         plot(Imaging_split{1}{1}.time_restricted(lap_idxs.A(ii,1):lap_idxs.A(ii,2))/60,...
%             Behavior_split{1}{1}.resampled.position(lap_idxs.A(ii,1):lap_idxs.A(ii,2)),...
%             'Color',[0 0 1 0.6],'LineWidth',1.5)
%     end
%     %B laps
%     for ii=1:size(lap_idxs.B,1)
%         plot(Imaging_split{1}{2}.time_restricted(lap_idxs.B(ii,1):lap_idxs.B(ii,2))/60,...
%             Behavior_split{1}{2}.resampled.position(lap_idxs.B(ii,1):lap_idxs.B(ii,2)),...
%             'Color',[1 0 0 0.6],'LineWidth',1.5)
%     end
%     %overlay significant calcium run events
%     %A
%     scatter(event_norm_time.A{ROI},event_norm_pos_run.A{ROI},[],[0 0 1],'*')
%     %B
%     scatter(event_norm_time.B{ROI},event_norm_pos_run.B{ROI},[],[1 0 0],'*')
%     
%     %generate run and non-run patched across track
%     if 0
%     xRange = {5:10; 20; 40:42; 50; 60:75; 90:95};
%     
%     for pr =1:size(Behavior_split{1, 1}{1, 1}.run_on_off_idx)
%     
%        xRange{pr,1} =  Imaging_split{1}{1}.time_restricted(Behavior_split{1}{1}.run_on_off_idx(pr,1): Behavior_split{1}{1}.run_on_off_idx(pr,2))./60;
%        
%     end
%     for k1 = 1:size(xRange,1)
%         q = xRange{k1};
%         %if very brief period
%         if length(q) == 1
%             q = [q q+0.1];
%         end
%         %start x end x; end x start x
%         qx = [min(q) max(q)  max(q)  min(q)];
%         yl = ylim;
%         %start y x2; end y x2
%         qy = [[1 1]*yl(1) [1 1]*yl(2)];
%         %plot the patches
%         %green
%         %patch(qx, qy, [0 1 0],'EdgeColor', 'none','FaceAlpha', 0.3)
%         %red
%         patch(qx, qy, [1 0 0],'EdgeColor', 'none','FaceAlpha', 0.3)
%     end
%     end
%     pause
%     clf
% end

end

