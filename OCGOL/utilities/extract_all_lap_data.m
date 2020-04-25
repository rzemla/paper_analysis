function [all_lap_data_transients] = extract_all_lap_data(session_vars)

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
% 2 – correct A trial
% 3 – correct B trial
% 20 – incorrect A trial
% 30 – incorrect B trial
% -10 – punish trial
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
    t_corr = session_vars{1, 1}.Behavior.performance.trialCorrect;
    %extract the trials for each session
    [all_data{ss}.trialCorr] = return_trial_corr(t_corr,all_data{ss}.lap_nb);
end

%% Make calcium transient matrix (binary - all events run and no run)
%RESUME HERE - last thing is to return the transients and QC them
%ADD OTHER VECTORs FOR EASE OF BROWSING IN THE FUTURE


%match list to all sessions
match_list = registered.multi.assigned_filtered;

%make blank 2D matrix for each each assigned and not assigned neuron for
%give session
transient_mat_all = nan(size(match_list,1),size(lap_nb,1));
%fill indexes
neuron_occup_log = ~isnan(match_list(:,1));
%ROIs corresponding to the matching indices
neuron_idx_ordered = match_list(neuron_occup_log,1);
%all the binary events corresponding to the imaging frames
event_onset_all = session_vars{1, 1}.Events.onset_binary';
%assign these event maps to the matching transient matrix
transient_mat_all(neuron_occup_log,:) = event_onset_all(neuron_idx_ordered,:);




end

