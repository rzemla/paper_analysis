function [session_vars] = extract_speed_for_each_trial_set(session_vars)

%across entire restricted session
run_epoch_all_laps = session_vars{1}.Behavior.run_ones;

%across every lap
run_epoch_each_lap = session_vars{1}.Behavior_split_lap.run_ones;

%speed across all complete laps (downsampled to match frames)
speed_all_laps = session_vars{1}.Behavior.speed;

%normalized position across all laps
norm_pos_all_laps = session_vars{1}.Behavior.resampled.normalizedposition;

%trial order of each lap
%2 - A lap
%3 - B lap
%added 0 - (i.e. - 20 or 30 is wrong lap is wrong A or B lap,respectively)
trialOrder = session_vars{1}.Behavior.performance.trialOrder;

%lap number of each frame index
lapNb_fr = session_vars{1}.Behavior.resampled.lapNb;

%% Split by lap
max_lap = max(lapNb_fr);

%get start and end lap idxs across (in frame space) for each lap
for ll=1:max_lap
    lap_idxs(ll,1) = find(lapNb_fr == ll,1,'first');
    lap_idxs(ll,2) = find(lapNb_fr == ll,1,'last');
end

%normalized position and speed on each lap
for ll=1:max_lap
   norm_pos_each_lap{ll} = norm_pos_all_laps(lap_idxs(ll,1):lap_idxs(ll,2));  
   speed_each_lap{ll} = speed_all_laps(lap_idxs(ll,1):lap_idxs(ll,2));
end


%% 

%update Behavior.split_lap.position (lap by lap)
%correct A laps
session_vars{1}.Behavior_split{1}.speed = cell2mat(speed_each_lap(find(trialOrder == 2))');
%get run only speed
session_vars{1}.Behavior_split{1}.speed_runEpoch = session_vars{1}.Behavior_split{1}.speed(logical(session_vars{1}.Behavior_split{1}.run_ones));

%correct B laps
session_vars{1}.Behavior_split{2}.speed = cell2mat(speed_each_lap(find(trialOrder == 3))');
%get run only speed
session_vars{1}.Behavior_split{2}.speed_runEpoch = session_vars{1}.Behavior_split{2}.speed(logical(session_vars{1}.Behavior_split{2}.run_ones));

%all laps
session_vars{1}.Behavior_split{3}.speed = cell2mat(speed_each_lap');
%get run only speed
session_vars{1}.Behavior_split{3}.speed_runEpoch = session_vars{1}.Behavior_split{3}.speed(logical(session_vars{1}.Behavior_split{3}.run_ones));

%all A laps
session_vars{1}.Behavior_split{4}.speed = cell2mat(speed_each_lap(find(trialOrder == 2 | trialOrder == 20))');
%get run only speed
session_vars{1}.Behavior_split{4}.speed_runEpoch = session_vars{1}.Behavior_split{4}.speed(logical(session_vars{1}.Behavior_split{4}.run_ones));

%all B laps
session_vars{1}.Behavior_split{5}.speed = cell2mat(speed_each_lap(find(trialOrder == 3 | trialOrder == 30))');
%get run only speed
session_vars{1}.Behavior_split{5}.speed_runEpoch = session_vars{1}.Behavior_split{5}.speed(logical(session_vars{1}.Behavior_split{5}.run_ones));


end

