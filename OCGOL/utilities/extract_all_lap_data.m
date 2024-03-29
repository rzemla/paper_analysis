function [all_data] = extract_all_lap_data(session_vars,registered)

%complete  lap data only
%get all of the data below for each session 
%time vector - ok
%run state vector - ok
%lap idx vector - ok
%position vector  in cm - ok
%position vector normalized - ok
%speed vector - ok
%lick signal (binary) vector - ok
%trial type vector (2  = A), (3 = B), -10 - punish, -1 = technical - ok
%trial correct vector 1 = correct; 0 - incorrect, nan = punish or technical
%- ok
%transients (binary) matrix ROI idx x time(frames) - these need to be
%extracted from the match indexes -ok

%run a QC to pull random ROIs from matching and transients and see if the
%vectors pulled out are identical


%% Extract all available data vectors that are already constructed
for ss=1:size(session_vars,2)
    %time resampled
    all_data{ss}.time = session_vars{ss}.Behavior.resampled.time;
    %run state
    all_data{ss}.run_state = session_vars{ss}.Behavior.run_ones;
    %position (cm)
    all_data{ss}.position_cm = session_vars{ss}.Behavior.resampled.position;
    %position norm
    all_data{ss}.position_norm = session_vars{ss}.Behavior.resampled.normalizedposition;
    %speed
    all_data{ss}.speed = session_vars{ss}.Behavior.speed;
    %lap number (complete) of each frame
    all_data{ss}.lap_nb = session_vars{ss}.Behavior.resampled.lapNb;
end

%% Construct binary lick vector for each session

%run for each session - see what data correspondence is from other dataset
%(rate map export)
for ss=1:size(session_vars,2)
    [all_data{ss}.lick] = return_lick_vec(session_vars,ss);
end

%% Extract trial order data

%legend for trialOrder and trialCorrect vector

% session_vars{}.Behavior.performance.trialOrder 
% 2 � correct A trial
% 3 � correct B trial
% 20 � incorrect A trial
% 30 � incorrect B trial
% -10 � punish trial
% -1 - techical exclusion

% session_vars{}.Behavior.performance.trialCorrect 
% 0 means incorrect trial
% 1 mean correct trial
% -2 means punish lap
%-1 means technical exclusion

for ss=1:size(session_vars,2)
    %find A or B trials (20 and 2) - A incorr/corr and (30 and 3) - B incorr/corr
    t_order = session_vars{ss}.Behavior.performance.trialOrder;
    %extract the trials for each session
    [all_data{ss}.trialType] = return_trial_order(t_order,all_data{ss}.lap_nb);
end

%% Vector for each frames indicating whether frame in in correct state

for ss=1:size(session_vars,2)
    %correct vs. incorrect vector
    %1 = correct lap; 0 = incorrect lap; nan = technical or punish lap
    t_corr = session_vars{ss}.Behavior.performance.trialCorrect;
    %extract the trials for each session
    [all_data{ss}.trialCorr] = return_trial_corr(t_corr,all_data{ss}.lap_nb);
end

%% Make calcium transient matrix (binary - all events run and no run)

%match list to all sessions (ROI match list 
match_list = registered.multi.assigned_filtered;

for ss=1:size(session_vars,2)
    %all the binary events corresponding to the imaging frames
    event_onset_all = session_vars{ss}.Events.onset_binary';
    [all_data{ss}.transient_mat] = return_transient_matrix(event_onset_all,match_list,all_data{ss}.lap_nb,ss);
end

%% QC check - plot and see how the data aligns
%{
%take random ROIs from match list, get transient vector and compare
ROI_nb_match = 150;
ses_nb = 2;
ROI_test = match_list(ROI_nb_match,ses_nb);
event_test_mat = session_vars{ses_nb}.Events.onset_binary';

%check if the match indexes correspond to the the correct ROI on that day
isequal(event_test_mat(ROI_test,:),all_data{ses_nb}.transient_mat(ROI_nb_match,:));
%}

%% Export additional lap descriptive data

for ss=1:size(session_vars,2)
    %lap time onset/offset - complete data
    all_data{ss}.lap_data.lap_onset_offset = session_vars{ss}.Behavior.lap;
    %lap/trial type based on reward signal and trial signal
    all_data{ss}.lap_data.id = session_vars{ss}.Behavior.lap_id;
    %performance data merge with lap
    all_data{ss}.lap_data.performance = session_vars{ss}.Behavior.performance;
end


%% Plot all and sub-section of data for any given dataset

%choose random ROI
ROI_nb = randi(size(match_list,1));

%for each session
for ss=1:size(session_vars,2)
    figure
    subplot(5,1,1)
    hold on
    title('Normalized postion and Ca2+ transients (binary)')
    xlabel('Time [s]')
    plot(all_data{ss}.time,all_data{ss}.position_norm)
    plot(all_data{ss}.time,all_data{ss}.transient_mat(ROI_nb,:))
    
    subplot(5,1,2)
    hold on
    title('Normalized postion and lick')
    xlabel('Time [s]')
    plot(all_data{ss}.time,all_data{ss}.position_norm)
    plot(all_data{ss}.time,all_data{ss}.lick)
    
    subplot(5,1,3)
    hold on
    title('Normalized postion and run state')
    xlabel('Time [s]')
    plot(all_data{ss}.time,all_data{ss}.position_norm)
    plot(all_data{ss}.time,all_data{ss}.run_state)
    
    subplot(5,1,4)
    hold on
    title('Speed')
    ylabel('cm/s')
    xlabel('Time [s]')
    plot(all_data{ss}.time,all_data{ss}.speed)
    
    subplot(5,1,5)
    hold on
    ylim([-10, 3])
    title('Trial type and trial correct')
    xlabel('Time [s]')
    plot(all_data{ss}.time,all_data{ss}.trialType)
    %plot animal position between trial correct and trial type
    plot(all_data{ss}.time,all_data{ss}.position_norm+1)
    %correct or incorrect trial
    plot(all_data{ss}.time,all_data{ss}.trialCorr)
    
end


end

