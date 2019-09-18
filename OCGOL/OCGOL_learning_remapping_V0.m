%% Import variables and define options

%run componenet registration across sessions
options.register = 0;

%whether to load place field data processed below
options.loadPlaceField_data = 1;

%load extracted ROI zooms/outlines
options.load_ROI_zooms_outlines = 1;

%visualize ROI outlines of matches across sessions
options.visualize_match = 0;

%lab workstation
%input directories to matching function
%  path_dir = {'G:\OCGOL_learning_short_term\I56_RTLS\I56_RLTS_5AB_041019_1',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_5AB_041119_2',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_3A3B_041219_3',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_3A3B_041319_4',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_no_punish_041519_5',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_no_punish_041619_6'};
% %cross session directory
% crossdir = 'G:\OCGOL_learning_short_term\I56_RTLS\crossSession';

% %I57_RTLS
 path_dir = {'G:\OCGOL_learning_short_term\I57_RTLS\I57_RLTS_5AB_041019_1',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_5AB_041119_2',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_3A3B_041219_3',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_1A1B_041319_4',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_no_punish_041519_5',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_no_punish_041619_6'};
     %'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_punish_041719_7'};
%cross session directory
crossdir = 'G:\OCGOL_learning_short_term\I57_RTLS\crossSession';

%I57_LT
%  path_dir = {'G:\OCGOL_learning_short_term\I57_LT\I57_LT_5A5B_041619_1',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_5A5B_041719_2',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_3A3B_041819_3',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_3A3B_041919_4',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_no_punish_042019_5',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_no_punish_042119_6'};
%  %,...
%      %'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_punish_042219_7',...
%      %'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_punish_042319_8'};
% %cross session directory
% crossdir = 'G:\OCGOL_learning_short_term\I57_LT\crossSession';

%% Load place cell variables for each session
%get mat directories in each output folder
for ii=1:size(path_dir,2)
    %get matfile names for each session
    matfiles{ii} = dir([path_dir{ii},'\output','\*.mat']);
end

%load in place cell variables (and others later)
for ii = [1 2 3 4 5 6]%1:size(path_dir,2)
    %add event variables
    disp(ii)
    %decide which variables here do not need to be loaded
    session_vars{ii} = load(fullfile(matfiles{ii}.folder,matfiles{ii}.name),'Place_cell', 'Behavior',...
        'Behavior_split_lap','Behavior_split','Events_split','Events_split_lap', 'Imaging_split');
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

for ss = [1 2 3 4 5 6]%1:size(session_vars,2) %1,2,3,4,5,6 OK
    %for ss= [4]
    disp(['Running session: ', num2str(ss)]);
    for ii = [4,5]
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
for ss = [1 2 3 4 5 6]
    for tt=4:5
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
    for ss = [1 2 3 4 5 6]
        %all A and all B
        for tt=4:5
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
options.sessionSelect = [1 2 3 4 5 6];
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
options.selectTrial = [4 5];
%which sessions to run
options.sessionSelect = [1 2 3 4 5 6];
%continue to modify 
[field_event_rates,pf_vector,field_total_events, select_fields] = transient_rate_in_field_multi_ses(session_vars,options);

%% Get max transient peak here
%get field event rates of max peak

[max_bin_rate,max_transient_peak] = max_transient_rate_multi_ses(session_vars,field_event_rates,pf_vector,options);

%% Filter filtered matching components for SI or TS tuning for at least on id'd place field and 5 events in firld

%which trials to use to calculate the in field transient rate
options.selectTrial = [4 5];
%which session to include in calculation
options.sessionSelect = [1 2 3 4 5 6];
%select fields has logical 1 for whichever neurons has a place field at at
%least 5 events on distinct laps within that PF - otherwise not PF
[registered] = filter_matching_components(registered,tunedLogical,select_fields,options);

%% Centroid difference (max transient rate)
%MODIFY TO BEHAVE LIKE FOR SINGLE SESSIONS

options.tuning_criterion = 'ts';
%which trials to use 
options.selectTrial = [4 5];
%which session to include in calculation
options.sessionSelect = [1 2 3 4 5 6];
centroid_diff_multi_ses(session_vars,tunedLogical, pf_vector,field_event_rates,select_fields,registered,options)


%% Split neurons by A or B task selective category - A or B selective (exclusive)
%which criterion to use for task-selective ROIs
%ts or both - ts selects only selective neurons based on TS tuning
%criterion
%both - uses both SI and TS criterion to select selectiven neurons
options.tuning_criterion = 'both';
%display events vs position for each task selective neuron in A or B
options.dispFigure = 0;
[task_selective_ROIs] = task_selective_categorize_multi_ses(tunedLogical,session_vars, max_transient_peak,options);


%quick look at fractions
for ss=1:6
counts_task_sel_each(1,ss) = length(task_selective_ROIs{ss}.A.idx)/size(session_vars{ss}.Place_cell{4}.Spatial_Info.ROI_pvalue,2)
counts_task_sel_each(2,ss) = length(task_selective_ROIs{ss}.B.idx)/size(session_vars{ss}.Place_cell{4}.Spatial_Info.ROI_pvalue,2)
end

%divide by total neurons each day

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
[task_remapping_ROIs,partial_field_idx] = remapping_categorize(cent_diff, tunedLogical, pf_vector_max ,pf_vector, session_vars,...
                        max_transient_peak,pf_count_filtered_log, pf_count_filtered,select_fields,options);

%% PV and TC correlations for all matching neurons (PV) in A and B trials across days (line plot); TC corr (for A tuned or B tuned on both days)

%set option as to how to select neurons for plots
options.tuning_criterion = 'si'; %si or ts
options.sessionSelect = [1 2 3 4 5 6 ];
options.selectSes = [4 5];
%learning or recall datasets
options.learning_data = 1;
[PV_TC_corr] = PV_TC_corr_across_days(session_vars,tunedLogical,registered,options);

%save to output file for cumulative analysis
save(fullfile(crossdir,'PV_TC_corr.mat'),'PV_TC_corr')

%% All neuron (at least 2 match between sessions) raster (non_norm)

%set option as to how to select neurons for plots
options.tuning_criterion = 'si'; %si or ts
options.sessionSelect = [1 2 3 4 5 6];
%chose all A/B (learning) vs. only correct A/B (recall)
options.selectTrial = [4,5];
%is it a learning set (for plot/raster annotation)
options.learning_data = 1;
non_norm_matching_STC_rasters(session_vars,tunedLogical,registered,options,crossdir)

%% Extract performance fractions across sessions (respective laps)
%check if agree with manual analysis
%turn into table with future code upgrade

%which sessions to use
options.sessionSelect = [1 2 3 4 5 6];

[ses_perf,ses_lap_ct] = session_performance(session_vars,options);

%export session performance data
save(fullfile(crossdir,'ses_perf.mat'),'ses_perf','ses_lap_ct');


%% Plot smoothed event rate across track (function)

%% Raster spiral - prepare multi session

[plot_raster_vars] = prepare_inputs_raster_spiral_multi_ses(session_vars);


%% Plot spiral raster using inputs

ROI_categories.task_selective_ROIs = task_selective_ROIs;

%works with inputs
plot_raster_spiral_multi_ses(plot_raster_vars,session_vars,registered,ROI_zooms, ROI_outlines,ROI_categories,options)


%% Define tuned logical vectors

%flag to all A or B trial or only correct A or B trials
options.allCorrect = 0;
%select which session to use
options.sessionSelect = [1 2 3 4 5 6];
%returns struct of structs
[tunedLogical] = defineTunedLogicals(session_vars,options);

%% Calculate the transient rates in each of the place fields (integrate later) and recalculate centroids based on highest transient rate field

%funtion to calculate transient rate in field
%take raw event rate and divide by occupancy (s) transients/s
% TODO: see ifrecalculate to see if dividing by normalized occupancy (fractional 0-1)
%yields different result
%[field_event_rates,pf_vector] = transient_rate_in_field(session_vars);

options.selectSes = [4 5];
%continue to modify 
[field_event_rates,pf_vector,field_total_events, select_fields] = transient_rate_in_field_multi_ses(session_vars,options);


%% Tuning specificity differences pre-learning vs. post-learning (matching neurons)

TS_score_diff(session_vars,tunedLogical,registered)

%% Plot spatial tuning curves according to transient rate 

options.tuning_criterion = 'ts'; %si or ts
plot_STC_transient_rate(session_vars,tunedLogical,registered,field_event_rates, pf_vector,options)



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


%% Learning pre and post event spiral
%display event spiral
%dF/F across laps between unlearned and learned day
%ROI outlines that show matching between days - need day from match ROI
%script


%%% CODE FROM GOL TASK BELOW %%%%

%% Calculate the distance from centroid to reward (function)

%compare to which reward (obsolete)
%if 1 --> GOL Block 1 location; if 2 --> GOL Block 2
options.rewardBlock = 1;

[centroids] = centroid_diff_reward(session_vars, options);


%Danielson
%(for cells with multiple fields
%the field with the highest in-field transient rate was used to calculate the centroid)
% The place f e fraction of complete forward passes through the
% place field associated with a significant Ca2+ transient) was significantly higher in deep than in superficial (n = 14
% mice, p = 0.001, paired T-Test), but the (ii) place cell specificity (defined as the fraction of running-related transients
% occurring within the place field)

%% Plot mean centroid diff to reward for each session against performance (D vs. P) 

%fraction of licks in reward zone for each session
for ss=1:size(session_vars,2)
    frac_rew(ss) = session_vars{ss}.Behavior.performance.frac_rew;
end

%mean of centroid diff to reward of each session for TS sig neurons
for ss=1:size(session_vars,2)
    mean_Distance_TS_sig(ss) = nanmean(centroids.TS_sig.angleDiffReward{ss});
end

%make scatter plot
figure;
hold on
%ylim([0 0.8]);
%xlim([0.5 2.25]);
ylabel('Fraction of licks in reward zone')
xticks([0.7854, 1.1781,1.5708, 1.9635])
xticklabels({'\pi/4','3\pi/8','\pi/2','5\pi/8'})
xlabel('Mean distance to reward');
scatter(mean_Distance_TS_sig,frac_rew,'b');

%% Bar graph of distance to reward on 1st vs last day
%only tuned cells according to TS tuning criteria

figure;
hold on
ylabel('Distance to reward');
xticks([1 2 3 4 5 6 7]);
xticklabels({'RF Day 0','GOL 1 Day 1','GOL 1 Day 2','GOL 1 Day 3',...
    'GOL 2 Day 4','GOL 2 Day 5','GOL 2 Day 6'})
xtickangle(45)
yticks([0.7854,1.5708])
ylim([0.7854,2])
yticklabels({'\pi/4','\pi/2'})
bar(mean_Distance_TS_sig,'b')
refline(0,1.5708)


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

