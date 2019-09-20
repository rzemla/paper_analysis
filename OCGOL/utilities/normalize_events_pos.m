function [norm_events] = normalize_events_pos(session_vars,options)

%% Define session/trial inputs
selectTrial = options.selectTrial;
sessionSelect = options.sessionSelect;

%% Load behavioral data for each session

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


%% Event onsets in run interval
%use trial definition set by trialSelect
for ss=sessionSelect
    %for each ROI
    for rr=1:size(events{ss}{1},2) %get the ROI size
        %time of significant run events in A
        event_norm_time{ss}.A{rr} = Imaging_split{ss}{selectTrial(1)}.time_restricted(find(Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr) == 1))/60;
        %normalizesd position of significant run events in A
        event_norm_pos_run{ss}.A{rr} = Behavior_split{ss}{selectTrial(1)}.resampled.position_norm(find(Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr) == 1));
        %lap assignment for A
        event_lap_idx{ss}.A{rr} = Behavior_split{ss}{selectTrial(1)}.resampled.lapNb(logical(Event_split{ss}{selectTrial(1)}.Run.run_onset_binary(:,rr)));
        
        %time of significant run events in B
        event_norm_time{ss}.B{rr} = Imaging_split{ss}{selectTrial(2)}.time_restricted(find(Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr) == 1))/60;
        %normalizesd position of significant run events in B
        event_norm_pos_run{ss}.B{rr} = Behavior_split{ss}{selectTrial(2)}.resampled.position_norm(find(Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr) == 1));
        %lap assignment for B
        event_lap_idx{ss}.B{rr} = Behavior_split{ss}{selectTrial(2)}.resampled.lapNb(logical(Event_split{ss}{selectTrial(2)}.Run.run_onset_binary(:,rr)));
    end
end

%% Export
norm_events.event_norm_time = event_norm_time;
norm_events.event_norm_pos_run = event_norm_pos_run;
norm_events.event_lap_idx =  event_lap_idx;

end

