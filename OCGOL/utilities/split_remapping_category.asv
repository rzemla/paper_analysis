function [outputArg1,outputArg2] = split_remapping_category(cent_diff_AandB, tuned_log, pf_vector_max, session_vars,max_transient_peak)
%split mutually tuned neurons by remapping category: 
%common (less than certain centroid difference between max
%tuned_log = tunedLogical.ts.AandB_tuned;

%% Get ROI indices
AandB_tuned_idx = find(tuned_log ==1);

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

%% Extract place fields corresponding to max transient rate


max_transient_peak

Place_cell{1}{1}.placeField.edge{1, 1}

max_placeField_edge
max_placeField_center

%% Get lap indices for each lap in all B or B trials

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
%only correct trials (1,2)
%for each ROI
for rr=1:size(events{ss}{1},2)
    %time of significant run events in A
    event_norm_time.A{rr} = Imaging_split{1}{1}.time_restricted(find(Event_split{1}{1}.Run.run_onset_binary(:,rr) == 1))/60;
    %normalizesd position of significant run events in A
    event_norm_pos_run.A{rr} = Behavior_split{1}{1}.resampled.position_norm(find(Event_split{1}{1}.Run.run_onset_binary(:,rr) == 1));
    %time of significant run events in B
    event_norm_time.B{rr} = Imaging_split{1}{2}.time_restricted(find(Event_split{1}{2}.Run.run_onset_binary(:,rr) == 1))/60;
    %normalizesd position of significant run events in B
    event_norm_pos_run.B{rr} = Behavior_split{1}{2}.resampled.position_norm(find(Event_split{1}{2}.Run.run_onset_binary(:,rr) == 1));
end


%% Plot the run epochs, corresponding position bin edges of place field

%plot normalized position
figure('Position', [1930 130 1890 420])
for rr=1:size(AandB_tuned_idx,2)
    ROI = AandB_tuned_idx(rr)
    hold on
    yticks([0 0.5 1])
    ylabel('Normalized position')
    xlabel('Time [min]');
    xticks(0:3:12);
    ylim([0 1])
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
    
    pause
    clf
end

%spatial bin assignment for each run-epoch frame (100 bins) 
session_vars{1}.Place_cell{1}.Bin{8}
%get edges for corresponding bins





%% Patch generator for run epochs
%creates x range of patches
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

 
 
%% Check that at least 5 calcium events (on distinct laps?) occured in field with max rate


%% Check that in running epoch within 3 bins to the left or right of each


end

