function OCGOL_recall_remapping_Jason_export_STC_shortened_V0(path_dir,crossdir)

%% Import variables and define options

%run componenet registration across sessions
options.register = 0;

%whether to load place field data processed below
options.loadPlaceField_data = 1;

%load extracted ROI zooms/outlines
options.load_ROI_zooms_outlines = 1;

%visualize ROI outlines of matches across sessions
options.visualize_match = 0;

%load SCE data shuffled n=50/100 (re-shuffle later on cluster with n =1000)
options.loadSCE = 0;

options.selectTrial = [1 2];

%moved below to calculate based on number of input sessions
%which session to include in calculation
%options.sessionSelect = [1 2 3 4 5 6];

%for use in workspace
%selectTrial = options.selectTrial;
%sessionSelect = options.sessionSelect;

%ANIMAL #1
%I42R 1
% path_dir = {'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d1_032118_1',...
%     'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d2_032218_2',...
%     'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d3_032318_3',...
%     'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d6_032618_4_2',...
%     'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d7_032718_5',...
%     'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d8_032818_6',...
%     'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d9_032918_7'};
% %cross session directory
% crossdir = 'G:\OCGOL_stability_recall\I42R_1\crossSession';

%ANIMAL #2
% I46
%  path_dir = {'G:\OCGOL_stability_recall\I46\I46_AB_d1_062018_1',...
%      'G:\OCGOL_stability_recall\I46\I46_AB_d2_062118_2',...
%      'G:\OCGOL_stability_recall\I46\I46_AB_d3_062218_3',...
%      'G:\OCGOL_stability_recall\I46\I46_AB_d6_062518_4',...
%      'G:\OCGOL_stability_recall\I46\I46_AB_d7_062618_5',...
%      'G:\OCGOL_stability_recall\I46\I46_AB_d8_062718_6',...
%      'G:\OCGOL_stability_recall\I46\I46_AB_d9_062818_7'};
% %cross session directory
% crossdir = 'G:\OCGOL_stability_recall\I46\crossSession';

%ANIMAL #3
%I45
% path_dir = {'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d1_062018_1',...
%     'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d2_062118_2',...
%     'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d3_062218_3',...
%     'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d6_062518_4',...
%     'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d7_062618_5',...
%     'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d8_062718_6',...
%     'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d9_062818_7'};
% %cross session directory
% crossdir = 'G:\OCGOL_stability_recall\I45_RT\crossSession';

%ANIMAL #4
%I42L 1
%  path_dir = {'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d1_032118_1',...
%      'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d2_032218_2',...
%      'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d3_032318_3',...
%      'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d6_032618_4',...
%      'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d7_032718_5',...
%      'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d8_032818_6',...
%      'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d9_032918_7'};
% %cross session directory
% crossdir = 'G:\OCGOL_stability_recall\I42L_1\crossSession';

%ANIMAL #5
%I47 LP
% path_dir = {'G:\OCGOL_stability_recall\I47_LP\I47_LP_AB_d1_062018_1',...
%      'G:\OCGOL_stability_recall\I47_LP\I47_LP_AB_d2_062118_2',...
%      'G:\OCGOL_stability_recall\I47_LP\I47_LP_AB_d3_062218_3',...
%      'G:\OCGOL_stability_recall\I47_LP\I47_LP_AB_d6_062518_4',...
%      'G:\OCGOL_stability_recall\I47_LP\I47_LP_AB_d7_062618_5',...
%      'G:\OCGOL_stability_recall\I47_LP\I47_LP_AB_d8_062718_6',...
%      'G:\OCGOL_stability_recall\I47_LP\I47_LP_AB_d9_062818_7'};
% %cross session directory
% crossdir = 'G:\OCGOL_stability_recall\I47_LP\crossSession';


                            %%%%% LONG-TERM (30 day sessions) %%%%
%I47_LP
% path_dir = {'D:\OCGOL_learning_long_term\I47_LP\I47_LP_AB_d1_062018_1',...
%      'D:\OCGOL_learning_long_term\I47_LP\I47_LP_AB_d6_062518_2',...
%      'D:\OCGOL_learning_long_term\I47_LP\I47_LP_AB_d16_070518_3',...
%      'D:\OCGOL_learning_long_term\I47_LP\I47_LP_AB_d20_070918_4',...
%      'D:\OCGOL_learning_long_term\I47_LP\I47_LP_AB_d25_071318_5',...
%      'D:\OCGOL_learning_long_term\I47_LP\I47_LP_AB_d30_071818_6'};
% %cross session directory
% crossdir = 'D:\OCGOL_learning_long_term\I47_LP\crossSession';

%I45 RT
% path_dir = {'D:\OCGOL_learning_long_term\I45_RT\I45_RT_AB_d1_062018_1',...
%      'D:\OCGOL_learning_long_term\I45_RT\I45_RT_AB_d6_062518_2',...
%      'D:\OCGOL_learning_long_term\I45_RT\I45_RT_AB_d16_070518_3',...
%      'D:\OCGOL_learning_long_term\I45_RT\I45_RT_AB_d20_070918_4',...
%      'D:\OCGOL_learning_long_term\I45_RT\I45_RT_AB_d25_071318_5',...
%      'D:\OCGOL_learning_long_term\I45_RT\I45_RT_AB_d30_071818_6'};
% %cross session directory
% crossdir = 'D:\OCGOL_learning_long_term\I45_RT\crossSession';

%I46
% path_dir = {'D:\OCGOL_learning_long_term\I46\I46_AB_d1_062018_1',...
%      'D:\OCGOL_learning_long_term\I46\I46_AB_d6_062518_2',...
%      'D:\OCGOL_learning_long_term\I46\I46_AB_d16_070518_3',...
%      'D:\OCGOL_learning_long_term\I46\I46_AB_d20_070918_4',...
%      'D:\OCGOL_learning_long_term\I46\I46_AB_d25_071318_5',...
%      'D:\OCGOL_learning_long_term\I46\I46_AB_d30_071818_6'};
% %cross session directory
% crossdir = 'D:\OCGOL_learning_long_term\I46\crossSession';

%% Determine number of sessions to analyze for animal

%which session to include in calculation
options.sessionSelect = 1:size(path_dir,2);

%for use as var in global workspace
sessionSelect = options.sessionSelect;

%% Override to selection (only 1st for animal 1)
% options.sessionSelect = 1;
% sessionSelect = 1;

%% Load place cell variables for each session
%get mat directories in each output folder
for ii=options.sessionSelect
    %get matfile names for each session
    matfiles{ii} = dir([path_dir{ii},'\output','\*.mat']);
    
    %for first day, read in the cleaned ROIs for recurrence analysis
    if ii==1
        rm_d1_files = dir([path_dir{ii},'\removedROI_clean','\*.mat']);
    end
end

%load in place cell variables (and others later)
for ii = options.sessionSelect%1:size(path_dir,2)
    %add event variables
    disp(ii)
    %decide which variables here do not need to be loaded
    session_vars{ii} = load(fullfile(matfiles{ii}.folder,matfiles{ii}.name),'Place_cell', 'Behavior',...
        'Behavior_split_lap','Behavior_split','Events_split','Events_split_lap', 'Imaging_split');
    
    %read in cleaned soma ROIs from day 1 (only first session)
    if ii == 1
        load(fullfile(rm_d1_files.folder,rm_d1_files.name),'removedROI_clean')
    end
end

%load additional data from struct
for ii = options.sessionSelect
    %add event variables
    disp(ii)
    session_vars_append{ii} = load(fullfile(matfiles{ii}.folder,matfiles{ii}.name),'Imaging','updated_dff','Events');
    %session_vars_append{ii} = load(fullfile(matfiles{ii}.folder,matfiles{ii}.name),'Events');   
end

%assign to main session variable struct (additional variables)
for ii = options.sessionSelect
    session_vars{ii}.Imaging = session_vars_append{ii}.Imaging;
    session_vars{ii}.updated_dff = session_vars_append{ii}.updated_dff;
    session_vars{ii}.Events = session_vars_append{ii}.Events;   
end

%% Match ROIs from across OCGOL days

if options.register == 1
    %run cross registration
    disp('Running registration of components');
    [registered] = match_ROIs_V2(path_dir,crossdir);
    
    %save registration variable in crosssession
    disp('Saving registered component matchings');
    save(fullfile(crossdir,'registered.mat'),'registered');
    
elseif options.register == 0
    %load the registered struct
    disp('Loading registered component matchings...');
    load(fullfile(crossdir,'registered.mat'),'registered');
    disp('Loaded.');
    
end

%% Load filtered ROI matches into registered struct

%get dir path with wildcard match to .mat files
filtered_ROI_dir_path = subdir(fullfile(crossdir,'filtered_match_ROI','*.mat'));
%load in temp var
match_var = load(filtered_ROI_dir_path.name);
%load in registered struct
registered.multi.assigned_filtered = match_var.ROI_assign_multi_filtered;

%% Get ROI_zooms and ROI_outlines for each neuron on each day
%number of sessions (runs even if not all session vars are loaded)
%already soma parsed
nbSes = size(session_vars,2);

%extract and save ROI zooms/outlines for all neurons
if options.load_ROI_zooms_outlines == 0
    %calculate ROI zooms/outlines
    [ROI_zooms, ROI_outlines] = defineOutlines_eachSes(nbSes,session_vars, path_dir);
    %save to directory
    save(fullfile(crossdir,'ROI_zooms_outlines.mat'),'ROI_zooms','ROI_outlines');
else %load from save files
    load(fullfile(crossdir,'ROI_zooms_outlines.mat'),'ROI_zooms','ROI_outlines');
    disp('Loaded ROI zooms and outlines.')
end

%% Visualize the matching ROIs that were matched above (match on every session only!)
%number of ROIs (rows) by sessions (cols)
rows = 20;
cols = 6; %take # of sessions as input

%number of sessions to look at
nb_ses = cols;

if options.visualize_match ==1
    visualize_matches_filtered(rows,cols,registered,ROI_zooms,ROI_outlines,nb_ses,crossdir);
end


%% Calculate relevant place fields

if options.loadPlaceField_data == 0
%use rate map - number of event onsets/ occupancy across all laps
options.gSigma = 3;
%which place cell struct to do placefield extraction on
%iterate through place_cell cells of interest
%4 - all A regardless if correct
%5 - all B regardless if correct
%I57 RTLS - problem with 4,4 - fixed
%I57 LT - problem with ses 4, trial 5 adjust (set to -2) - narrow as opposed to
%extend field - apply to rest of animals

for ss =options.sessionSelect%1:size(session_vars,2) %1,2,3,4,5,6 OK
    %for ss= [4]
    disp(['Running session: ', num2str(ss)]);
    for ii = options.selectTrial
        options.place_struct_nb = ii;
        disp(['Running trial type: ', num2str(ii)]);
        [session_vars{ss}.Place_cell] = place_field_finder_gaussian(session_vars{ss}.Place_cell,options);
    end
end

%save whole place cell struct and load in and replace for each session in
%the future
%make post-processing directory (postProcess)
mkdir(crossdir,'postProcess')

%for each Place_cell session extract placeField struct
%use trial types here
for ss = options.sessionSelect
    for tt=options.selectTrial
        session_pf{ss}(tt).placeField = session_vars{ss}.Place_cell{tt}.placeField;
    end
end

%save Place_cell struct in that directory
save(fullfile(crossdir,'postProcess','placeField_upd_struct.mat'),'session_pf')

else
    tic;
    disp('Loading place field data')
    load(fullfile(crossdir,'postProcess','placeField_upd_struct.mat'));
    toc
    %replace the Place_cell struct in the session_vars cell
    for ss = options.sessionSelect
        %all A and all B
        for tt=options.selectTrial
            session_vars{ss}.Place_cell{tt}.placeField = session_pf{ss}(tt).placeField;
        end
    end
end

%% Define tuned logical vectors

%flag to all A or B trial or only correct A or B trials
%all correct = 0 ==> uses trials 4,5
%all correct = 1 ==> uses trials 1,2
options.allCorrect = 1;

%which sessions to run
%options.sessionSelect = [1 2 3 4 5 6];
%returns struct of structs
[tunedLogical] = defineTunedLogicals(session_vars,options);

%% Calculate the transient rates in each of the place fields (integrate later) and recalculate centroids based on highest transient rate field

%funtion to calculate transient rate in field
%take raw event rate and divide by occupancy (s) transients/s
% TODO: see ifrecalculate to see if dividing by normalized occupancy (fractional 0-1)
%yields different result
%[field_event_rates,pf_vector] = transient_rate_in_field(session_vars);

%which trials to use to calculate the in field transient rate
%[1 2] - only correct A B trials
%[4 5] - all A B trials
%A correct/B correct or all
%options.selectTrial = [4 5];
%which sessions to run
%options.sessionSelect = [1 2 3 4 5 6];
%continue to modify 
[field_event_rates,pf_vector,field_total_events, select_fields] = transient_rate_in_field_multi_ses(session_vars,options);


%% Get max transient peak here
%get field event rates of max peak

%select whether to use TS vector of adjusted vector for cells with single fields
%originally used 0, but make switch to 1
options.select_adj_vec = 0;

[max_bin_rate,max_transient_peak] = max_transient_rate_multi_ses(session_vars,field_event_rates,pf_vector,options);


%% Filter filtered matching components for SI or TS tuning for at least on id'd place field and 5 events in firld

%which trials to use to calculate the in field transient rate
%options.selectTrial = [4 5];
%which session to include in calculation
%options.sessionSelect = [1 2 3 4 5 6];
%select fields has logical 1 for whichever neurons has a place field at at
%least 5 events on distinct laps within that PF - otherwise not PF
[registered] = filter_matching_components(registered,tunedLogical,select_fields,options);

%% Recurrence and fraction active analysis

[recurr,frac_active,recurr_ex,frac_active_ex] = recurrence_analysis(registered,removedROI_clean,session_vars,tunedLogical,select_fields,options);

%save recurrence and fraction active of neurons

save(fullfile(crossdir,'recurrence.mat'),'recurr','frac_active','recurr_ex','frac_active_ex');

%% Centroid difference (max transient rate) - also returns vector of max angle (of place field) - pf_vector_max
%MODIFY TO BEHAVE LIKE FOR SINGLE SESSIONS

options.tuning_criterion = 'ts';
%which trials to use 
%options.selectTrial = [4 5];
%which session to include in calculation
%options.sessionSelect = [1 2 3 4 5 6];
[cent_diff,pf_vector_max] = centroid_diff_multi_ses(session_vars,tunedLogical, pf_vector,field_event_rates,select_fields,registered,options);

%save pf_vector_max and cent_diff for angle difference analysis

save(fullfile(crossdir,'pf_vector_max.mat'),'pf_vector_max');

%% Export matching STCs for Jason for detecting splitter cells across time
%QC checked - works
%all neurons first with minimum number of events

[matching_tun_curves] = export_day_matched_STCs_recall(session_vars,session_vars_append,registered);

%export the matching STCs
save(fullfile(crossdir,'matching_tun_curves.mat'),'matching_tun_curves');


%% Extract performance fractions across sessions (respective laps)
%check if agree with manual analysis
%turn into table with future code upgrade

%which sessions to use
%options.sessionSelect = [1 2 3 4 5 6];

%performance and total laps for each session
[ses_perf,ses_lap_ct] = session_performance(session_vars,options);

%modify into table format with descriptors
%1st row - performance on all laps
%2nd row - A trial performance
%3rd row - B trial performance

%equivalent session days
day_recall_number = [1 2 3 6 7 8 9];

for ii=1:size(session_vars,2)
    var_names{ii} = ['Session_',num2str(ii),'_', 'Day_', num2str(day_recall_number(ii))];
end

%convert into table
ses_perf_table = array2table(ses_perf,'RowNames',{'All trials','A trials','B trials'},'VariableNames',var_names);
ses_lap_ct_table = array2table(ses_lap_ct,'RowNames',{'All trials','A trials','B trials'},'VariableNames',var_names);

%export performance tables for Jason
save(fullfile(crossdir,'perf_lap_tables.mat'),'ses_perf_table','ses_lap_ct_table');

%export session performance data
save(fullfile(crossdir,'ses_perf.mat'),'ses_perf','ses_lap_ct');


end
