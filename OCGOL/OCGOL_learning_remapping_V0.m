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

%all A and B trials used for learning (sessions determined below)
options.selectTrial = [4 5];

%for use in global workspace
selectTrial = options.selectTrial;

%I56_RTLS
%input directories to matching function
%  path_dir = {'G:\OCGOL_learning_short_term\I56_RTLS\I56_RLTS_5AB_041019_1',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_5AB_041119_2',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_3A3B_041219_3',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_3A3B_041319_4',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_no_punish_041519_5',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_no_punish_041619_6',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_punish_041719_7'};
% % %cross session directory
%  crossdir = 'G:\OCGOL_learning_short_term\I56_RTLS\crossSession';

%I57_RTLS
 path_dir = {'G:\OCGOL_learning_short_term\I57_RTLS\I57_RLTS_5AB_041019_1',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_5AB_041119_2',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_3A3B_041219_3',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_1A1B_041319_4',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_no_punish_041519_5',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_no_punish_041619_6',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_punish_041719_7'};

%cross session directory
crossdir = 'G:\OCGOL_learning_short_term\I57_RTLS\crossSession';

%I57_LT
%  path_dir = {'G:\OCGOL_learning_short_term\I57_LT\I57_LT_5A5B_041619_1',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_5A5B_041719_2',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_3A3B_041819_3',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_3A3B_041919_4',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_no_punish_042019_5',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_no_punish_042119_6',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_punish_042219_7',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_punish_042319_8'};
%  
%  crossdir = 'G:\OCGOL_learning_short_term\I57_LT\crossSession';

%I58 RT
%input directories to matching function
%  path_dir = {'E:\OCGOL_learning_short_term\I58_RT\I58_RT_5A5B_073019_1',...
%      'E:\OCGOL_learning_short_term\I58_RT\I58_RT_5A5B_073119_2',...
%      'E:\OCGOL_learning_short_term\I58_RT\I58_RT_3A3B_080119_3',...
%      'E:\OCGOL_learning_short_term\I58_RT\I58_RT_3A3B_080219_4',...
%      'E:\OCGOL_learning_short_term\I58_RT\I58_RT_randAB_no_punish_080319_5',...
%      'E:\OCGOL_learning_short_term\I58_RT\I58_RT_randAB_no_punish_080419_6',...
%      'E:\OCGOL_learning_short_term\I58_RT\I58_RT_ABrand_punish_080519_7',...
%      'E:\OCGOL_learning_short_term\I58_RT\I58_RT_randAB_punish_080619_8',...
%      'E:\OCGOL_learning_short_term\I58_RT\I58_RT_randAB_punish_080719_9'};
% %cross session directory
% crossdir = 'E:\OCGOL_learning_short_term\I58_RT\crossSession';

%I58 LT
% path_dir = {'E:\OCGOL_learning_short_term\I58_LT\I58_LT_5A5B_080419_1',...
%      'E:\OCGOL_learning_short_term\I58_LT\I58_LT_5A5B_080519_2',...
%      'E:\OCGOL_learning_short_term\I58_LT\I58_LT_3A3B_080619_3',...
%      'E:\OCGOL_learning_short_term\I58_LT\I58_LT_3A3B_080719_4',...
%      'E:\OCGOL_learning_short_term\I58_LT\I58_LT_randAB_no_punish_080819_5',...
%      'E:\OCGOL_learning_short_term\I58_LT\I58_LT_randAB_no_punish_080919_6',...
%      'E:\OCGOL_learning_short_term\I58_LT\I58_LT_randAB_punish_081119_7'};
% %cross session directory
% crossdir = 'E:\OCGOL_learning_short_term\I58_LT\crossSession';

%I58 RTLP
%  path_dir = {'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_5A5B_080419_1',...
%      'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_5A5B_080519_2',...
%      'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_3A3B_080619_3',...
%      'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_3A3B_080719_4',...
%      'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_randAB_no_punish_080819_5',...
%      'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_randAB_no_punish_080919_6'};
% %cross session directory
% crossdir = 'E:\OCGOL_learning_short_term\I58_RTLP\crossSession';

%% Determine number of sessions to analyze for animal

%which session to include in calculation
options.sessionSelect = 1:size(path_dir,2);

%for use as var in global workspace
sessionSelect = options.sessionSelect;


%% Load place cell variables for each session
%get mat directories in each output folder
for ii=options.sessionSelect
    %get matfile names for each session
    matfiles{ii} = dir([path_dir{ii},'\output','\*.mat']);
end

%load in place cell variables (and others later)
for ii = options.sessionSelect%1:size(path_dir,2)
    %add event variables
    disp(ii)
    %decide which variables here do not need to be loaded
    session_vars{ii} = load(fullfile(matfiles{ii}.folder,matfiles{ii}.name),'Place_cell', 'Behavior',...
        'Behavior_split_lap','Behavior_split','Events_split','Events_split_lap', 'Imaging_split');
end

%load additional data from struct
for ii = options.sessionSelect
    %add event variables
    disp(ii)
    session_vars_append{ii} = load(fullfile(matfiles{ii}.folder,matfiles{ii}.name),'Imaging','updated_dff','Events');
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


%% Insert save checkpoint here to avoid re-preprocessing above data
if 0
%create cross session processed/loaded data directory
mkdir(fullfile(crossdir,'cross_data'))

%save the loaded and processed componenet/place cell data
disp('Saving place field, session variables, and component matching data');
save(fullfile(crossdir,'cross_data','cross_loaded.mat'),'session_vars','ROI_outlines','ROI_zooms','registered','-v7.3');
end

%% Define tuned logical vectors

%flag to all A or B trial or only correct A or B trials
%all correct = 0 ==> uses trials 4,5
%all correct = 1 ==> uses trials 1,2
options.allCorrect = 0;
%select which session to use
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
%adjust center vector to point to the field
options.select_adj_vec = 1;
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

%set to blank for learning analysis
removedROI_clean = [];

[recurr,frac_active,recurr_ex,frac_active_ex] = recurrence_analysis(registered,removedROI_clean,session_vars,tunedLogical,select_fields,options);

%save recurrence and fraction active of neurons

save(fullfile(crossdir,'recurrence.mat'),'recurr','frac_active','recurr_ex','frac_active_ex');


%% Centroid difference (max transient rate)
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
%all neurons first

export_day_matched_STCs(session_vars,session_vars_append,registered)

%% Split neurons by A or B task selective category - A or B selective (exclusive)
%which criterion to use for task-selective ROIs
%ts or both - ts selects only selective neurons based on TS tuning
%criterion
%both - uses both SI and TS criterion to select selectiven neurons
options.tuning_criterion = 'both';
%display events vs position for each task selective neuron in A or B
options.dispFigure = 0;
[task_selective_ROIs] = task_selective_categorize_multi_ses(tunedLogical,session_vars, max_transient_peak,options);


% %quick look at fractions
% for ss=1:6
% counts_task_sel_each(1,ss) = length(task_selective_ROIs{ss}.A.idx)/size(session_vars{ss}.Place_cell{4}.Spatial_Info.ROI_pvalue,2)
% counts_task_sel_each(2,ss) = length(task_selective_ROIs{ss}.B.idx)/size(session_vars{ss}.Place_cell{4}.Spatial_Info.ROI_pvalue,2)
% end
% 
% figure;
% subplot(1,2,1)
% hold on
% ylim([0 0.4])
% % plot(both_count')
% % plot(sum(both_count,1))
% subplot(1,2,2)
% hold on
% ylim([0 0.4])
% plot(counts_task_sel_each')
% plot(sum(counts_task_sel_each,1))


%% Number of place fields and widths for each sub-class of neurons
%add filter for classfing whether each field is significant (min 5 events)

%doesn;t matter
options.tuning_criterion = 'si'; %si or ts
%A correct/B correct or all
%options.selectTrial = [4 5];
%options.sessionSelect = [1 2 3 4 5 6];
[placeField_dist, pf_count_filtered_log, pf_count_filtered] = placeField_properties_multi_ses(session_vars,tunedLogical,select_fields,task_selective_ROIs,options);
%save the place field distributions output data
%save(fullfile(path_dir{1},'cumul_analysis','placeField_dist.mat'),'placeField_dist');

%% Calculate fraction tuned S.I. vs T.S for every session

[tuned_fractions,tuned_logicals] = fractionTuned_multi_ses(tunedLogical,pf_count_filtered_log,options);

%save fractional count
save(fullfile(crossdir,'tuned_fractions.mat'),'tuned_fractions');
%save tuned logicals
save(fullfile(crossdir,'tuned_logicals.mat'),'tuned_logicals');

%% Look at spatial information scores in matching neurons between days in A trials and B trials
%put into separate script
 si_score_comparison(session_vars, registered)

%% Task remapping filter - split into remapping categories
%which criterion to use for task-selective ROIs
options.tuning_criterion = 'ts';
%display events vs position for each task selective neuron in A or B
options.dispFigure = 0;
%number of degrees of centroid difference
%45 deg ~25 cm; 
%36 deg ~20 cm;
%27 deg ~15 cm;
%18 dege ~10 cm
options.deg_thres = 18;
%ranges for splitting the global remappers
%0-10 cm; 10 - 30cm; 30+ cm
options.deg_ranges = [0 18 54];
%degree threshold for partial remappers
options.partial_deg_thres = [18 36];
%choice between KS test of unpaired Mann Whitney U (later)
%either 'ranksum' or ks
options.AUC_test = 'ranksum';
%significance level of test
options.p_sig = 0.05;
%make sure that this function does not overwrite the the previous
%task_selective_ROIs structure

[task_remapping_ROIs,partial_field_idx] = remapping_categorize_multi_ses(cent_diff, tunedLogical ,pf_vector, session_vars,...
                        max_transient_peak, pf_count_filtered,select_fields,options);

%calculate total A&B neurons for each session
for ss=options.sessionSelect
    remap_cat_count(ss,1) = length(task_remapping_ROIs{ss}.global_near);
    remap_cat_count(ss,2) = length(task_remapping_ROIs{ss}.global_far);
    remap_cat_count(ss,3) = length(task_remapping_ROIs{ss}.rate);
    remap_cat_count(ss,4) = length(task_remapping_ROIs{ss}.common);
    remap_cat_count(ss,5) = length(task_remapping_ROIs{ss}.partial);
    remap_cat_count(ss,6) = length(task_remapping_ROIs{ss}.mixed);
    task_remapping_ROIs{ss}.nbROI = sum(remap_cat_count(ss,:));
end

fraction_rel_AB = remap_cat_count./repmat(sum(remap_cat_count,2),1,6);

%convert to fraction
figure
hold on
plot(fraction_rel_AB(:,6))

%% Save task-selective and task remapping neurons into struct neurons 

%get # of ROIs in each session
for ss=sessionSelect
    ses_nbROI(ss) = size(session_vars{ss}.Place_cell{selectTrial(1)}.Tuned_ROI_mask,2);
end

save(fullfile(crossdir,'task_neurons.mat'),'task_selective_ROIs','task_remapping_ROIs','ses_nbROI');

%% Load selective neurons
load(fullfile(crossdir,'task_neurons.mat'),'task_selective_ROIs','task_remapping_ROIs','ses_nbROI');

%% Assign each matching neuron to remapping category

match_mat = registered.multi.assigned_filtered;
%make cell with categorical values
cat_registered_cell = cell(size(match_mat,1),size(match_mat,2));

%for each session determine which category each cell belongs to
for ss=options.sessionSelect
    [~,idx_match] = intersect(match_mat(:,ss),task_selective_ROIs{ss}.A.idx);
    cat_registered_cell(idx_match,ss) = {'A-selective'};
    [~,idx_match] = intersect(match_mat(:,ss),task_selective_ROIs{ss}.B.idx);
    cat_registered_cell(idx_match,ss) = {'B-selective'};
    [~,idx_match] = intersect(match_mat(:,ss),task_remapping_ROIs{ss}.common);
    cat_registered_cell(idx_match,ss) = {'common'};
    [~,idx_match] = intersect(match_mat(:,ss),task_remapping_ROIs{ss}.rate);
    cat_registered_cell(idx_match,ss) = {'rate'};
    [~,idx_match] = intersect(match_mat(:,ss),task_remapping_ROIs{ss}.global_near);
    cat_registered_cell(idx_match,ss) = {'near'};
    [~,idx_match] = intersect(match_mat(:,ss),task_remapping_ROIs{ss}.global_far);
    cat_registered_cell(idx_match,ss) = {'far'};
    [~,idx_match] = intersect(match_mat(:,ss),task_remapping_ROIs{ss}.partial);
    cat_registered_cell(idx_match,ss) = {'partial'};
    [~,idx_match] = intersect(match_mat(:,ss),task_remapping_ROIs{ss}.mixed);
    cat_registered_cell(idx_match,ss) = {'mixed'};
end

%extract 2 session index
ses_comp = [4,5];
selMatchIdxs = find(sum(~isnan(match_mat(:,ses_comp)),2)==2);

cat_registered_cell(selMatchIdxs,ses_comp);

%% Extract normalized events
tic;
[norm_events] = normalize_events_pos(session_vars,options);
toc;

%% Get AUC/min calculation for each ROI in A vs. B
%AUC/min and frequency of events events/min in run epoch for A or B trials
%all added no run epochs

[session_vars] = AUC_rate(session_vars,options);

%% Raster spiral - prepare inputs for multi session display
% save this in the future as well and load
options.spiral_width = 0.1;
[plot_raster_vars] = prepare_inputs_raster_spiral_multi_ses(session_vars,options);

%% Raster, STC, ROI image across days
%insert dummy SCE data

plot_raster_spiral_STC_ROI_multi_ses(plot_raster_vars,session_vars,registered,cat_registered_cell,ROI_zooms,ROI_outlines,options)

%% Decicated ALL session spiral plotter with categorical type display


plot_raster_spiral_multi_ses_label_check(plot_raster_vars,session_vars,registered,cat_registered_cell,multi_ses_SCE_data,options)


%% Decicated two session spiral plotter with categorical type display

%insert dummy SCE data
multi_ses_SCE_data = [];
plot_raster_spiral_multi_ses_label_check(plot_raster_vars,session_vars,registered,cat_registered_cell,multi_ses_SCE_data,options)


%% Decicated two session spiral plotter with categorical type display
%keep working on this split out other functionalities in other functions
plot_spiral_raster_SCE(plot_raster_vars,session_vars,registered,cat_registered_cell,SCE,options)



%% Visualize place fields and events for each neurons
%use this to check code detection accuracy

visualize_neuron_characteristics(plot_raster_vars,norm_events,registered,session_vars,task_selective_ROIs,cat_registered_cell,select_fields,options)



%% Plot spiral raster using inputs
%ploy 
ROI_categories.task_selective_ROIs = task_selective_ROIs;

%works with inputs
plot_raster_spiral_multi_ses(plot_raster_vars,session_vars,registered,ROI_zooms, ROI_outlines,ROI_categories,options)

%% Task-selective neurons - AUC/min and event freq distributions

%code for rose plot like event density is here
activity_distributions(session_vars,task_selective_ROIs,options)

%% S.I. score for matching neurons across days
%resume later
track_si_score_change(session_vars,registered,options)


%% TC correlation for matching (+/-) tuned ROIs

[tc_corr_match] = tc_corr_matching_neurons(session_vars,registered,options);

%save to output file for cumulative analysis
save(fullfile(crossdir,'tc_corr_match.mat'),'tc_corr_match')

%% Extract performance fractions across sessions (respective laps)
%check if agree with manual analysis
%turn into table with future code upgrade

%which sessions to use
%options.sessionSelect = [1 2 3 4 5 6];

[ses_perf,ses_lap_ct] = session_performance(session_vars,options);

%export session performance data
save(fullfile(crossdir,'ses_perf.mat'),'ses_perf','ses_lap_ct');

%% PV and TC correlations for all matching neurons (PV) in A and B trials across days (line plot); TC corr (for A tuned or B tuned on both days)

%set option as to how to select neurons for plots
options.tuning_criterion = 'si'; %si or ts
%options.sessionSelect = [1 2 3 4 5 6 ];
options.selectSes = [4 5];
%learning or recall datasets
options.learning_data = 1;
[PV_TC_corr] = PV_TC_corr_across_days(session_vars,tunedLogical,registered,options);

%save to output file for cumulative analysis
save(fullfile(crossdir,'PV_TC_corr.mat'),'PV_TC_corr')

%% Raster plot across days (FIGURE 4) (at least 2 match between sessions) raster (non_norm)
%this generates the figure rasters in 4

%set option as to how to select neurons for plots
options.tuning_criterion = 'si'; %si or ts
%options.sessionSelect = [1 2 3 4 5 6];
%chose all A/B (learning) vs. only correct A/B (recall)
options.selectTrial = [4,5];
%is it a learning set (for plot/raster annotation)
options.learning_data = 1;
non_norm_matching_STC_rasters(session_vars,tunedLogical,registered,options,crossdir)


%% PV and TC correlations for all matching neurons (PV) in A and B trials across days (line plot); TC corr (for A tuned or B tuned on both days)

%set option as to how to select neurons for plots
options.tuning_criterion = 'si'; %si or ts
%options.sessionSelect = [1 2 3 4 5 6 ];
options.selectSes = [4 5];
%learning or recall datasets (used for tick labels)
options.learning_data = 0;
[PV_TC_corr] = PV_TC_corr_across_days(session_vars,tunedLogical,registered,options);

%save to output file for cumulative analysis
save(fullfile(crossdir,'PV_TC_corr.mat'),'PV_TC_corr')

%% Plot meanTC histogram for all SCEs across days of learning

%match matrix
matching_ROI_matrix = registered.multi.assigned_filtered(:,1:6);

%create matching A neuron SCE participation matrix
SCE_A_ROI_engage = zeros(size(matching_ROI_matrix,1), size(matching_ROI_matrix,2));

%create matching B neuron SCE participation matrix
SCE_B_ROI_engage = zeros(size(matching_ROI_matrix,1), size(matching_ROI_matrix,2));

%ROI participation in A or B trials
for ss=1:6
    SCE_part_A{ss} = sum(SCE{ss}.sce_activity.A,2);
    SCE_part_B{ss} = sum(SCE{ss}.sce_activity.B,2);
    SCE_part_all{ss} = sum(SCE{ss}.sce_activity_matrix  ,2);
end

%A 
for ss=1:6
    assign_counts = SCE_part_A{ss}(matching_ROI_matrix(~isnan(matching_ROI_matrix(:,ss)),ss));
    SCE_A_ROI_engage(~isnan(matching_ROI_matrix(:,ss)),ss) = assign_counts
    SCE_A_ROI_engage(isnan(matching_ROI_matrix(:,ss)),ss) = nan;
end

%B
for ss=1:6
    assign_counts = SCE_part_B{ss}(matching_ROI_matrix(~isnan(matching_ROI_matrix(:,ss)),ss));
    SCE_B_ROI_engage(~isnan(matching_ROI_matrix(:,ss)),ss) = assign_counts
    SCE_B_ROI_engage(isnan(matching_ROI_matrix(:,ss)),ss) = nan;
end

%all
for ss=1:6
    assign_counts = SCE_part_all{ss}(matching_ROI_matrix(~isnan(matching_ROI_matrix(:,ss)),ss));
    SCE_all_ROI_engage(~isnan(matching_ROI_matrix(:,ss)),ss) = assign_counts
    SCE_all_ROI_engage(isnan(matching_ROI_matrix(:,ss)),ss) = nan;
end

multi_ses_SCE_data.SCE_A_ROI_engage = SCE_A_ROI_engage;
multi_ses_SCE_data.SCE_B_ROI_engage = SCE_B_ROI_engage;
multi_ses_SCE_data.SCE_all_ROI_engage = SCE_all_ROI_engage;

figure
hold on
plot(SCE_A_ROI_engage')

%% Detect SCEs and measure number of SCE in each session A or B

%load or run shuffle to determine SCEs
if options.loadSCE == 0
    %how many shuffles to perform
    options.shuffle_nb =50;
    [SCE] = detect_SCE(session_vars,options);
    
    %extract number of SCEs on each day (total)
    for ss=sessionSelect
        SCE_total_count(ss) = SCE{ss}.nbSCE;
    end
    
elseif options.loadSCE == 1
    %export session performance data
    load(fullfile(crossdir,'SCE.mat'),'SCE');
end

%% Number of neurons in each SCE
for ss=sessionSelect
    %for each SCE,
    for cc =1:size(SCE{ss}.sync_range,1)
        %all unique ROIs in each SCE
        SCE{ss}.SCE_unique_ROIs{cc} = unique(cell2mat(SCE{ss}.SCE_ROIs(SCE{ss}.sync_range(cc,1):SCE{ss}.sync_range(cc,2))));
        %number of ROIs in each SCE
        SCE{ss}.SCE_nbROI(cc) = size(SCE{ss}.SCE_unique_ROIs{cc},2);
    end
end

%generate x tick labels
for ss=sessionSelect
    name_cell{ss} = ['',num2str(ss)];
end

%get 
figure
hold on
for ss=sessionSelect
    cdfplot(SCE{ss}.SCE_nbROI);
end
legend(name_cell);

%% SCE onset order

[SCE] = sce_onset_order(session_vars,SCE,options);

%% Assign SCEs by trial type
%MODIFY TO TAKE IN wrong trials on learning sessions
[SCE] = assign_SCE_trials(session_vars,SCE,options);


%% Make SCE ROI participation matrix
%get SCE onset normalized position here as well
[SCE] = SCE_participation(session_vars,SCE,PV_TC_corr,options);

%% Save SCE struct for all sessions
%export session performance data
save(fullfile(crossdir,'SCE.mat'),'SCE');

%% Plot meanTC histogram for all SCEs across days of learning

%norm histogams
figure
for ss=1:6
    subplot(1,6,ss)
    hold on
    ylim([0 0.25])
    histogram(SCE{ss}.meanTC,[0:0.05:1],'Normalization','probability')
end

%plot cdf
figure
hold on
for ss=1:6
    ecdf(SCE{ss}.meanTC)
pause
end

%match matrix
matching_ROI_matrix = registered.multi.assigned_filtered(:,1:6);

%create matching A neuron SCE participation matrix
SCE_A_ROI_engage = zeros(size(matching_ROI_matrix,1), size(matching_ROI_matrix,2));

%create matching B neuron SCE participation matrix
SCE_B_ROI_engage = zeros(size(matching_ROI_matrix,1), size(matching_ROI_matrix,2));

%ROI participation in A or B trials
for ss=1:6
    SCE_part_A{ss} = sum(SCE{ss}.sce_activity.A,2);
    SCE_part_B{ss} = sum(SCE{ss}.sce_activity.B,2);
    SCE_part_all{ss} = sum(SCE{ss}.sce_activity_matrix  ,2);
end

%A 
for ss=1:6
    assign_counts = SCE_part_A{ss}(matching_ROI_matrix(~isnan(matching_ROI_matrix(:,ss)),ss));
    SCE_A_ROI_engage(~isnan(matching_ROI_matrix(:,ss)),ss) = assign_counts
    SCE_A_ROI_engage(isnan(matching_ROI_matrix(:,ss)),ss) = nan;
end

%B
for ss=1:6
    assign_counts = SCE_part_B{ss}(matching_ROI_matrix(~isnan(matching_ROI_matrix(:,ss)),ss));
    SCE_B_ROI_engage(~isnan(matching_ROI_matrix(:,ss)),ss) = assign_counts
    SCE_B_ROI_engage(isnan(matching_ROI_matrix(:,ss)),ss) = nan;
end

%all
for ss=1:6
    assign_counts = SCE_part_all{ss}(matching_ROI_matrix(~isnan(matching_ROI_matrix(:,ss)),ss));
    SCE_all_ROI_engage(~isnan(matching_ROI_matrix(:,ss)),ss) = assign_counts
    SCE_all_ROI_engage(isnan(matching_ROI_matrix(:,ss)),ss) = nan;
end

multi_ses_SCE_data.SCE_A_ROI_engage = SCE_A_ROI_engage;
multi_ses_SCE_data.SCE_B_ROI_engage = SCE_B_ROI_engage;
multi_ses_SCE_data.SCE_all_ROI_engage = SCE_all_ROI_engage;

%% SCE plots against session performance

%color maps -total, A, B
color_map = [139, 0, 139;  65,105,225; 220,20,60]/255;

figure
subplot(2,1,1)
hold on
ylabel('SCE count')
plot(SCE_total_count,'LineWidth',2)
xticks([1:6])
xticks([1:6])
xticklabels(name_cell);
xtickangle(45);
xlabel('Session');
set(gca,'FontSize',16)
set(gca,'LineWidth',2)

subplot(2,1,2)
hold on
ylabel('Performance')
xlabel('Session')

for ii=1:3
    p(ii) = plot(ses_perf(ii,:),'Color',color_map(ii,:),'LineWidth',2);
end
xticks([1:6])
xticklabels(name_cell);
xtickangle(45);
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
legend([p(1) p(2) p(3)],'All','A','B','Location','southeast')


%% Make SCE matrix
%row - neuron on/off
%column - each SCE
for ss=1:6
    ses_nbROI(ss) = size(session_vars{ss}.Place_cell{4}.Spatial_tuning_curve,2);
end

%create SCE matrix
sce_activity_matrix = zeros(ses_nbROI(ss),SCE{ss}.nbSCE);

%fill in each matrix
for cc=1:SCE{ss}.nbSCE
    sce_activity_matrix(SCE{ss}.SCE_unique_ROIs{cc},cc) = 1;
end

%% Detect SCE assemblies 
%max number of clusters in k-means
options.clust_max = 20;
detect_SCE_assembly(sce_activity_matrix,options)

%% Tuning specificity differences pre-learning vs. post-learning (matching neurons)

TS_score_diff(session_vars,tunedLogical,registered)

%% Generate dF/F maps of neurons tuned in either session and plot side by side
%customize to add options
%tuned in both sessions by SI score
%sorted by A trials

%set option as to how to select neurons for plots
%use SI or TS
options.tuning_criterion = 'ts'; %si or ts

%select which session to use
options.sessionSelect = [1 2 3 6];

plot_dFF_OCGOL_training(session_vars,tunedLogical,registered,options)

%% Generate STC maps of neurons tuned in either session and plot side by side
%customize to add options
%tuned in both sessions by SI score
%sorted by A trials

%set option as to how to select neurons for plots
options.tuning_criterion = 'ts'; %si or ts
options.sessionSelect = [1 2 3 6];
plot_STC_OCGOL_training(session_vars,tunedLogical,registered,options)

%% All neuron (at least 2 match between sessions) raster (non_norm)

%set option as to how to select neurons for plots
options.tuning_criterion = 'si'; %si or ts
options.sessionSelect = [1 2 3 4 5 6];
non_norm_matching_STC_rasters_learning(session_vars,tunedLogical,registered,options)


%% Measure PV and TC correlation between A/B trial on first training day and once learned
%expect greater dissimilarity once learned

%% Centroid distance between A and B trials on early vs. late training for matching ROIs

%% Event property difference (duration of sig events in fields vs

%% Place field width for A or B trials early vs late (all neurons vs matching neurons)


%% Plot spatial tuning curves according to transient rate 

% options.tuning_criterion = 'ts'; %si or ts
% plot_STC_transient_rate(session_vars,tunedLogical,registered,field_event_rates, pf_vector,options)


%% Calculate the transient rates in each of the place fields (integrate later) and recalculate centroids based on highest transient rate field

%funtion to calculate transient rate in field
%take raw event rate and divide by occupancy (s) transients/s
% TODO: see ifrecalculate to see if dividing by normalized occupancy (fractional 0-1)
%yields different result
%[field_event_rates,pf_vector] = transient_rate_in_field(session_vars);

% options.selectSes = [4 5];
% %continue to modify 
% [field_event_rates,pf_vector,field_total_events, select_fields] = transient_rate_in_field_multi_ses(session_vars,options);



%% Learning pre and post event spiral
%display event spiral
%dF/F across laps between unlearned and learned day
%ROI outlines that show matching between days - need day from match ROI
%script


%%% CODE FROM GOL TASK BELOW %%%%

%% Calculate the distance from centroid to reward (function)

%compare to which reward (obsolete)
%if 1 --> GOL Block 1 location; if 2 --> GOL Block 2
%options.rewardBlock = 1;

%[centroids] = centroid_diff_reward(session_vars, options);


%Danielson
%(for cells with multiple fields
%the field with the highest in-field transient rate was used to calculate the centroid)
% The place f e fraction of complete forward passes through the
% place field associated with a significant Ca2+ transient) was significantly higher in deep than in superficial (n = 14
% mice, p = 0.001, paired T-Test), but the (ii) place cell specificity (defined as the fraction of running-related transients
% occurring within the place field)

%% Plot mean centroid diff to reward for each session against performance (D vs. P) 

% %fraction of licks in reward zone for each session
% for ss=1:size(session_vars,2)
%     frac_rew(ss) = session_vars{ss}.Behavior.performance.frac_rew;
% end
% 
% %mean of centroid diff to reward of each session for TS sig neurons
% for ss=1:size(session_vars,2)
%     mean_Distance_TS_sig(ss) = nanmean(centroids.TS_sig.angleDiffReward{ss});
% end
% 
% %make scatter plot
% figure;
% hold on
% %ylim([0 0.8]);
% %xlim([0.5 2.25]);
% ylabel('Fraction of licks in reward zone')
% xticks([0.7854, 1.1781,1.5708, 1.9635])
% xticklabels({'\pi/4','3\pi/8','\pi/2','5\pi/8'})
% xlabel('Mean distance to reward');
% scatter(mean_Distance_TS_sig,frac_rew,'b');

%% Bar graph of distance to reward on 1st vs last day
%only tuned cells according to TS tuning criteria

% figure;
% hold on
% ylabel('Distance to reward');
% xticks([1 2 3 4 5 6 7]);
% xticklabels({'RF Day 0','GOL 1 Day 1','GOL 1 Day 2','GOL 1 Day 3',...
%     'GOL 2 Day 4','GOL 2 Day 5','GOL 2 Day 6'})
% xtickangle(45)
% yticks([0.7854,1.5708])
% ylim([0.7854,2])
% yticklabels({'\pi/4','\pi/2'})
% bar(mean_Distance_TS_sig,'b')
% refline(0,1.5708)


%% Raster plotting ROI notes

%nice examples from I56_RLTS training set
%27/50
%30/460
%34/55
%41/546
%54/76
%56/314
%58/475
%79/111
%82/115
%85/117
%93/130
%97/132
%99/135
%122/481
%135/167
%138/175
%143/186
%146/188
%152/192
%161/201
%166/210 - common to split cell!
%168/215 - shift to earlier point, followed by narrowing
%175/223 - task selective reward cell
%185/238
%190/241
%202/393
%216/268
%260/294 - precise narrowing
%263/287 
%264/290 - reward zone cells --> odor B responsive cell
%265/288 - stable field in one B trials; establishing additional field for A trials
%266/301- narrowing of fields
%267/295 - narrowing of fields
%268/332 - common to split cell!
%271/291 - precise narrowing
%272...
%274/311 - splitting of cell 
%279/286 - development of non-selective reward cell
%280/438 - narrowing of cell (very nice example!)
%284/309 - stable reward cell
%285/316 - nice examples of narrowing
%286/285 - emergence of B selective from no selective
%288/349 - narrowing
%293/284 - emergence of A selective peri-reward cell
%294/593 - emergence of A selective peri-odor cell
%296/282 - very nice narrowing example
%297/310 - multi-peaked place cells narrowing
%298/435 - flip cell (from B to A)
%300/344 - nice narrowing and remapping of A cell
%303/450 - tuning
%306/308 - reward cell tuning
%309/463 - tuning
%314/493 - reward tuning
%320/312 - tuning to A - scattered in B
%322/353 - flip from reward tuning to place cell tuned in A
%327/361 - dispering/broadening of activity
%330/392 - nice narrowing example 
%336/368 - narrowing
%337/206 - nice non-specific narrowing
%339/420 - nice convergence from separate to common activity 
%346/383 - nice narrowing of activity
%348/365 - nice split
%350/92 - nice development of place fields
%356/480 - split
%357/429 - nice A and B convergence
%359/379 - nice development of A field separate from B
%360/443 - remapping and divergence
%362/401 - convergence in 1 field; split in another
%374/371 - nice split
%412/452 - nice split
%413/511 - nice example
%421/484 - nice cleanup and split
%426/21 - non-selective peri-reward to selective perireward
%452/504 - nice split and remap
%465/303 - dual field non split
%474/534 - convergence

%cells that split 374, 474, 412, 348
%reward cell integration - 264


%separate processing part from display portion (speed reason)

%break down in simpler code in the future;
%interesting ROI I56_RTLS: 
%s5 - split over time
%s6 - start A sel and end A sel
%s14 - separate place field formation over time
%s20 - task selective place field formation over time

%session 1: 27 - selective, split, converge, diverge
%s1_31 - common and then diverge
%s1_41 - starts A and B and ends A&B
%45 -split
%47  - B sel deve
%54 - split
%58 common - split
%63 - split to common
%74 - flip contexts
%79 - narrow of place field (dF/F) and common -> split
%93 - common to selective
%95 - selective, split, common
%97 - common - split - B selective - split
%122 - common -> split
%128 - selective - split - back to selective
%135 - scatter - common 2 field
%138 - common - divergence - convergence
%143 - evolution of partial remapping neuron
%146 - selective - split - common
%152 - common to selective
%161 - task sel stable
%166 - common - split
%168 - split - common
%175 - split to task sel (B)
%181 - commom to common - same place
%185 - spread - converge to common
%186 - flip contexts
%190 - scatter to common
%207 - A sel - B sel - common
%216 - A sel B sel A sel (multiple selective flips)
%218 - common - B sel - split - maintain
%243 - scatter - A sel - split commong
%245 flip context
%246 split - maintain split + common
%263 - partial to common
%265 - stable split
%266 - stable
%268 common to split
%271 - scatter, split common
%279 - common, reward A , reward A in B, split
%280 - scatter, common,selective 
%285 - scatter -> common
%286 - scatter - narrow
%290 - selective common split partial
%293 - common sel A - split - sel A
%294 - non common split, conjunctive sel
%296 - scatter split comoon

%% Plot dF/F rasters side by side by laps for matching ROIs
%compare_sessions_raster_spiral(session_vars,registered,ROI_outlines,ROI_zooms, options);

