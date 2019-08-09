function [field_event_rates,pf_vector,field_total_events, select_fields] = transient_rate_in_field_multi_ses(session_vars,options)
%input - cell of structs with animal data

%% Define input variables 


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


%try for 1 session first

%for each session
for ss=1:size(session_vars,2)
    %for each trial (A or B) regardless if correct
    for tt=options.selectSes
        %Place field edge data
        placeField_edges{ss}{tt} = session_vars{ss}.Place_cell{tt}.placeField.edge;
        
        %Event rate for 100 bins
        event_rate{ss}{tt} = session_vars{ss}.Place_cell{tt}.Spatial_Info.rate_map{8};
        
        %event map (not rate)% 100 bins
        event_map{ss}{tt} = session_vars{ss}.Place_cell{tt}.Spatial_Info.event_map{8};
        
        %Occupancy for 100 bins (seconds, not normalized (0-1))
        occupancy{ss}{tt} = session_vars{ss}.Place_cell{tt}.Spatial_Info.occupancy_map{8};
    end
end

%% Get bin edges from normalized position

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

%% Event onsets in run interval
%only correct trials (1,2)
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

%should equal event rate (it does)
%rate_map = event_map./occupancy';
%isequal(rate_map,event_rate)

%% Find and store transient rate and total events in each place field

for ss=1:size(session_vars,2)
    %for each trial (A or B) regardless if correct (4,5) or only correct
    %(1,2)
    for tt=options.selectSes
        %find field event rate for each session and trial types
        [field_event_rates{ss}{tt}, field_total_events{ss}{tt}] = field_rate(event_map{ss}{tt},occupancy{ss}{tt},placeField_edges{ss}{tt});
    end
end


%% Make at at least 5 sig events on distinct laps
%get the edges of all neurons (all fields, not just max rate field)
%correct A
placeField_edge{1} = Place_cell{1}{options.selectSes(1)}.placeField.edge;
%correct B
placeField_edge{2} = Place_cell{1}{options.selectSes(2)}.placeField.edge;

%check which field has more than 1 field and select edges of the one with
%higher transient rate

%shouldn't need this code - just extracts the place fields with max
%transient rate
%for each trial (A and B )
% for tt=1:2
%     for rr=1:size(placeField_edge{tt},2)
%         %if more than 2 place fields
%         if size(placeField_edge{tt}{rr},1) > 1
%             placeField_filtered_max{tt}{rr} = placeField_edge{tt}{rr}(max_field_idx{tt}(rr),:);
%         else %if 1 place field
%             placeField_filtered_max{tt}{rr} = placeField_edge{tt}{rr};
%         end
%     end
% end

%convert edges from relevant place field to normalized postion edges
%edges are the normalized position equivalents of the bin edges identified
%in bin space (using 100 bins)

for tt=options.selectSes
    %for each ROI
    for rr=1:size(placeField_edge{tt},2)
        %start position of PF %edges 1 trial for both since the binning is
        %the same
        %for each place field
        if ~isempty(placeField_edge{tt}{rr})
            for pp=1:size(placeField_edge{tt}{rr},1)
                placeField_posnorm{tt}{rr}(pp,1) = edges{1}(placeField_edge{tt}{rr}(pp,1));
                %end position of PF
                placeField_posnorm{tt}{rr}(pp,2) = edges{1}(placeField_edge{tt}{rr}(pp,2)+1)-0.01;
            end
        else
            placeField_posnorm{tt}{rr} = [];
        end
        
    end
end

%% Filter out neurons that do not have at least 5 sig events in max place field in at least 5 distinct laps

%find events occuring within each place field for each ROI
for tt=options.selectSes
    for rr=1:size(placeField_posnorm{tt},2)
        if tt ==options.selectSes(1)%correct A trials (or all A trials)
            if ~isempty(placeField_posnorm{tt}{rr})
                %for each id'd place field
                for pp=1:size(placeField_posnorm{tt}{rr},1)
                    events_in_field{tt}{rr}{pp} = find(event_norm_pos_run.A{rr} >= placeField_posnorm{tt}{rr}(pp,1) & ...
                        event_norm_pos_run.A{rr} <= placeField_posnorm{tt}{rr}(pp,2));
                    %register the corresponding lap of in-field filtered event
                    event_in_field_laps{tt}{rr}{pp} = event_lap_idx.A{rr}(events_in_field{tt}{rr}{pp});
                    %get number of unique events (those occuring on each lap)
                    event_in_field_nb{tt}{rr}(pp) = size(unique(event_in_field_laps{tt}{rr}{pp}),1);
                    %get position of in-field events
                    events_in_field_pos{tt}{rr}{pp} = event_norm_pos_run.A{rr}(events_in_field{tt}{rr}{pp});
                end
            else
                events_in_field{tt}{rr} = [];
                event_in_field_laps{tt}{rr} = [];
                event_in_field_nb{tt}{rr} = [];
                events_in_field_pos{tt}{rr} = [];
            end
            
        elseif tt == options.selectSes(2) %correct B trials (or all B trials)
            if ~isempty(placeField_posnorm{tt}{rr})
                %for each id'd place field
                for pp=1:size(placeField_posnorm{tt}{rr},1)
                    events_in_field{tt}{rr}{pp} = find(event_norm_pos_run.B{rr} >= placeField_posnorm{tt}{rr}(pp,1) & ...
                        event_norm_pos_run.B{rr} <= placeField_posnorm{tt}{rr}(pp,2));
                    %register the corresponding lap of in-field filtered event
                    event_in_field_laps{tt}{rr}{pp} = event_lap_idx.B{rr}(events_in_field{tt}{rr}{pp});
                    %get number of unique events (those occuring on each lap)
                    event_in_field_nb{tt}{rr}(pp) = size(unique(event_in_field_laps{tt}{rr}{pp}),1);
                    %get position of in-field events
                    events_in_field_pos{tt}{rr}{pp} = event_norm_pos_run.B{rr}(events_in_field{tt}{rr}{pp});
                end
            else
                events_in_field{tt}{rr} = [];
                event_in_field_laps{tt}{rr} = [];
                event_in_field_nb{tt}{rr} = [];
                events_in_field_pos{tt}{rr} = [];
            end
        end
    end
end

%% Create logical selection vectors for neurons with place fields that contain at least 5 events on distinct laps

%create logical for each ROI on each set of trial types
for tt=options.selectSes
    for rr=1:size(event_in_field_nb{tt},2)
            %get logical of place fields with more than 5 distinct calcium
            %events in field
            select_fields{tt}{rr} = event_in_field_nb{tt}{rr} >= 5; 
    end
end

%% Plot
%which session
session_nb = 1;

%histogram of rates in id'd fields
figure;
subplot(2,1,1)
hold on;
xlim([0 1])
title('In-field transient rates for all neurons - A trials');
histogram(cell2mat(field_event_rates{session_nb}{options.selectSes(1)}))
subplot(2,1,2)
hold on
xlim([0 1])
title('In-field transient rates for all neurons - B trials');
histogram(cell2mat(field_event_rates{session_nb}{options.selectSes(2)}))

%% Recalculate centroid based on peak with highest transient rate

for ss=1:size(session_vars,2)
    %for each trial (A or B) regardless if correct
    for tt=options.selectSes
        %tuning vectors for each ROI
        tun_vectors{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity.tuning_vector;
        %original sum vector
        tun_vector{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity.tuning_vector_specificity;
    end
end

% Turn into function
%input edges and tuning vector
%output - adjusted tuning vector

%make separate output that calculates the tuning vectors for all fields

options.pf.skipDisplay = 0;

for ss=1:size(session_vars,2)
    %for each trial (A or B) regardless if correct
    for tt=options.selectSes
        [pf_vector{ss}{tt}] = adjust_tuning_vector(tun_vectors{ss}{tt},tun_vector{ss}{tt},placeField_edges{ss}{tt},options);
    end
end

%% Plot rate map, place field edges - use for debug, final check

% figure;
% for rr=21:22
%     hold on;
%     %plot unsmoothed/non-normalized event rate
%     plot(event_rate(:,rr), 'k');
%     %plot identified place fields
%     for pp=1:size(placeField_edges{rr},1)
%         stem(placeField_edges{rr}(pp,:), [1,1], 'r')
%     end
%     
%     pause
%     clf
% end



%%

end

