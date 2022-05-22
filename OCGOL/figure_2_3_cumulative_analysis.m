%% Define experiment core directories
%screened and good datasets for analysis
%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_120618_1'}; %
%field rate error (after PSEM silencing
path_dir{1} = {'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1'};
path_dir{2} = {'G:\Figure_2_3_selective_remap\I42R_AB_d1_032118_1'};

path_dir{3} = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1'};
path_dir{4} = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_2'};

path_dir{5} = {'G:\Figure_2_3_selective_remap\I53LT_AB_sal_113018_1'}; 
path_dir{6} = {'G:\Figure_2_3_selective_remap\I56_RTLS_AB_prePost_sal_042419_1'};

path_dir{7} = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_113018_1'};
path_dir{8} = {'G:\Figure_2_3_selective_remap\I57_RTLS_AB_prePost_792_042519_1'};

path_dir{9} = {'G:\Figure_2_3_selective_remap\I45_RT_AB_d1_062018_1'};
path_dir{10} = {'G:\Figure_2_3_selective_remap\I46_AB_d1_062018_1'};
path_dir{11} = {'G:\Figure_2_3_selective_remap\I57_LT_ABrand_no_punish_042119_1'};

%path_dir = {'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I42R_AB_d1_032118_1'};

%path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_2'};
 
%path_dir = {'G:\Figure_2_3_selective_remap\I53LT_AB_sal_113018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I56_RTLS_AB_prePost_sal_042419_1'};

%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_113018_1'};  
%path_dir = {'G:\Figure_2_3_selective_remap\I57_RTLS_AB_prePost_792_042519_1'};

%path_dir = {'G:\Figure_2_3_selective_remap\I45_RT_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I46_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I57_LT_ABrand_no_punish_042119_1'};

%N =10; FOV = 11 (2 FOV from single animal)

%% Load the reward endzones calculated from script below (remapping_centroids_updated) into variable

%load in struct with reward zone start and end for each animal
load(fullfile('G:\Figure_2_3_selective_remap\cumulative_data_output','reward_zones_all_animals.mat'));



%% Reviewer AUC re-analysis and nb events for each class of of task-selective place_cells

%get matfile path and filename
for ii=1:numel(path_dir)
    cd(fullfile(string(path_dir{ii}),'output'))
    source_path{ii} = dir ('*.mat');
end

%all data loaded in for AUC/min and event/min analysis
parfor ii=1:numel(path_dir)
cum_data{ii} = load(fullfile(source_path{ii}.folder, source_path{ii}.name),'Place_cell','Events_split',...
    'Behavior_split','Imaging_split')
disp(ii)
end

%nb events for each cell on A trials
cum_data{1, 1}.Events_split{1, 1}.Run.properties.nb_events

%% AUC comparison for reviewers

% Isolate AUC of all neurons
% Isolate AUC of SI tuned neurons
% Isolate AUC of TS tuned neurons
% Get time in min for run and non-run laps for all correct A and B trials
% (1,2)

%idx of SI place cells

%G:\Figure_2_3_selective_remap\I57_LT_ABrand_no_punish_042119_1\cumul_analysis\select_ROI_criteria
%.mat

%Read in the idx of all SI/TS tuned task-selective neurons
for ee=1:size(path_dir,2)
    load_data_select_ROIs{ee} = fullfile(path_dir{ee},'cumul_analysis','select_ROI_criteria.mat');
    ROI_idx{ee} = load(string(load_data_select_ROIs{ee}));
end

%% ROI index extraction of place cells with shared place fields on A and B trials with single PF
%only place cells with single fields tuned by SI or TS criteria
%single PF common place cells (non-task-selective) by SI or TS criteria -
%extract the indices of these cells for each animaa;
for ee=1:size(path_dir,2)
    %SI common (non-selective place cells with single PF)
    ROI_sel{ee}.si.common = intersect(find(ROI_idx{ee}.tunedLogical.si.AandB_tuned ==1), ROI_idx{ee}.task_remapping_ROIs.common);
    %TS common (non-selective place cells with single PF)
    ROI_sel{ee}.ts.common = intersect(find(ROI_idx{ee}.tunedLogical.ts.AandB_tuned ==1), ROI_idx{ee}.task_remapping_ROIs.common);

end

% selective vs. non-selective place cells 
%selective

%task-selective A and B as defined by paper (Including A and B selective
%neurons)
ROI_idx{ee}.task_selective_ROIs.A.idx
%place cells with single fields by SI or TS criteria - control for reviewer
ROI_sel{ee}.si.common
ROI_sel{ee}.ts.common

%A&B tuned as controls to task selective figure for transient rates for
%Figure 2 addressing reviewers comments
ROI_idx{1, 1}.tunedLogical.tuned_AB_si_ts

%% Extract AUC/min for each place cell for common single 
%logical list of SI/TS only A or only B tuned neurons
for ee=1:size(path_dir,2)
    %A_sel
    A_sel_idx.si{ee} = ROI_idx{1, ee}.tunedLogical.si.onlyA_tuned;
    A_sel_idx.ts{ee} = ROI_idx{1, ee}.tunedLogical.ts.onlyA_tuned;
    %B_sel
    B_sel_idx.si{ee} = ROI_idx{1, ee}.tunedLogical.si.onlyB_tuned;
    B_sel_idx.ts{ee} = ROI_idx{1, ee}.tunedLogical.ts.onlyB_tuned;

    %A&B tuned controls for comparison
    AB_idx{ee} = ROI_idx{1, ee}.tunedLogical.tuned_AB_si_ts;
end

%load AUC values run and no run events
for ee=1:size(path_dir,2)
    %run
    AUC_A{ee}.run = cum_data{1,ee}.Events_split{1, 1}.Run.properties.AUC;
    AUC_B{ee}.run =  cum_data{1,ee}.Events_split{1, 2}.Run.properties.AUC;
    %no run
    AUC_A{ee}.norun = cum_data{1,ee}.Events_split{1, 1}.NoRun.properties.AUC;
    AUC_B{ee}.norun =  cum_data{1,ee}.Events_split{1, 2}.NoRun.properties.AUC;    
end

%constant imaging period
dt= 0.033427969000002;

%% Extract run time for run/no-run epochs in min
for ee=1:size(path_dir,2)
    %correct A trials
    time_length{ee}.total.A = size(cum_data{1, ee}.Behavior_split{1, 1}.run_ones,1).*dt./60;
    time_length{ee}.run.A = sum(cum_data{1, ee}.Behavior_split{1, 1}.run_ones).*dt./60;
    time_length{ee}.norun.A = sum(~cum_data{1, ee}.Behavior_split{1, 1}.run_ones).*dt./60;

    %correct B trials
    time_length{ee}.total.B = size(cum_data{1, ee}.Behavior_split{1, 2}.run_ones,1).*dt./60;
    time_length{ee}.run.B = sum(cum_data{1, ee}.Behavior_split{1, 2}.run_ones).*dt./60;
    time_length{ee}.norun.B = sum(~cum_data{1, ee}.Behavior_split{1, 2}.run_ones).*dt./60;
end

%% A and B selective (and task selective) neurons by separate SI and TS criteria  AUC extraction
%animal --> run/no_run --> index of ROI
%RUN A and B selective by SI criteria
for ee=1:size(path_dir,2)

    %%% A-selective neurons %%%
    %A-sel run SI
    AUC_A_sel{ee}.si.run.A = AUC_A{1, ee}.run(A_sel_idx.si{ee});
    AUC_A_sel{ee}.si.run.B = AUC_B{1, ee}.run(A_sel_idx.si{ee});

    %A-sel norun SI
    AUC_A_sel{ee}.si.norun.A = AUC_A{1, ee}.norun(A_sel_idx.si{ee});
    AUC_A_sel{ee}.si.norun.B = AUC_B{1, ee}.norun(A_sel_idx.si{ee});

    %A-sel run TS
    AUC_A_sel{ee}.ts.run.A = AUC_A{1, ee}.run(A_sel_idx.ts{ee});
    AUC_A_sel{ee}.ts.run.B = AUC_B{1, ee}.run(A_sel_idx.ts{ee});

    %A-sel norun TS
    AUC_A_sel{ee}.ts.norun.A = AUC_A{1, ee}.norun(A_sel_idx.ts{ee});
    AUC_A_sel{ee}.ts.norun.B = AUC_B{1, ee}.norun(A_sel_idx.ts{ee});
 
    %%% B-selective neurons %%%
    %B-sel run SI
    AUC_B_sel{ee}.si.run.A = AUC_A{1, ee}.run(B_sel_idx.si{ee});
    AUC_B_sel{ee}.si.run.B = AUC_B{1, ee}.run(B_sel_idx.si{ee});

    %B-sel norun SI
    AUC_B_sel{ee}.si.norun.A = AUC_A{1, ee}.norun(B_sel_idx.si{ee});
    AUC_B_sel{ee}.si.norun.B = AUC_B{1, ee}.norun(B_sel_idx.si{ee});

    %B-sel run TS
    AUC_B_sel{ee}.ts.run.A = AUC_A{1, ee}.run(B_sel_idx.ts{ee});
    AUC_B_sel{ee}.ts.run.B = AUC_B{1, ee}.run(B_sel_idx.ts{ee});

    %B-sel norun TS
    AUC_B_sel{ee}.ts.norun.A = AUC_A{1, ee}.norun(B_sel_idx.ts{ee});
    AUC_B_sel{ee}.ts.norun.B = AUC_B{1, ee}.norun(B_sel_idx.ts{ee});


    %A -task selective neurons per figure 2 (temporary holding vars)
    A_sel_temp = ROI_idx{ee}.task_selective_ROIs.A.idx;
    B_sel_temp = ROI_idx{ee}.task_selective_ROIs.B.idx;

    %A-task selective run
    AUC_A_task_sel{ee}.run.A = AUC_A{1, ee}.run(A_sel_temp);
    AUC_A_task_sel{ee}.run.B = AUC_B{1, ee}.run(A_sel_temp);

    %A-task selective norun 
    AUC_A_task_sel{ee}.norun.A = AUC_A{1, ee}.norun(A_sel_temp);
    AUC_A_task_sel{ee}.norun.B = AUC_B{1, ee}.norun(A_sel_temp);

    %B-task selective run
    AUC_B_task_sel{ee}.run.A = AUC_A{1, ee}.run(B_sel_temp);
    AUC_B_task_sel{ee}.run.B = AUC_B{1, ee}.run(B_sel_temp);

    %B-task selective norun 
    AUC_B_task_sel{ee}.norun.A = AUC_A{1, ee}.norun(B_sel_temp);
    AUC_B_task_sel{ee}.norun.B = AUC_B{1, ee}.norun(B_sel_temp);

    %A&B tuned neurons as controls - run
    AUC_AB{ee}.run.A = AUC_A{1, ee}.run(AB_idx{ee});
    AUC_AB{ee}.run.B = AUC_B{1, ee}.run(AB_idx{ee});

    %A&B tuned neurons as controls - norun
    AUC_AB{ee}.norun.A = AUC_A{1, ee}.norun(AB_idx{ee});
    AUC_AB{ee}.norun.B = AUC_B{1, ee}.norun(AB_idx{ee});
    
    %AB common control - run SI
    AUC_common{ee}.si.run.A = AUC_A{1, ee}.run(ROI_sel{ee}.si.common);
    AUC_common{ee}.si.run.B = AUC_B{1, ee}.run(ROI_sel{ee}.si.common);

    %AB common control - norun SI
    AUC_common{ee}.si.norun.A = AUC_A{1, ee}.norun(ROI_sel{ee}.si.common);
    AUC_common{ee}.si.norun.B = AUC_B{1, ee}.norun(ROI_sel{ee}.si.common);    

    %AB common control - run TS
    AUC_common{ee}.ts.run.A = AUC_A{1, ee}.run(ROI_sel{ee}.ts.common);
    AUC_common{ee}.ts.run.B = AUC_B{1, ee}.run(ROI_sel{ee}.ts.common);

    %AB common control - norun TS
    AUC_common{ee}.ts.norun.A = AUC_A{1, ee}.norun(ROI_sel{ee}.ts.common);
    AUC_common{ee}.ts.norun.B = AUC_B{1, ee}.norun(ROI_sel{ee}.ts.common);     
end


%% Calculate AUC/min for common place fields

for ee=1:size(path_dir,2)
    %SI run
    AUC_min_ABcom{ee}.si.run.A = cell2mat(cellfun(@(x) sum(x),AUC_common{1, ee}.si.run.A  ,'UniformOutput',false))./time_length{1, ee}.run.A;
    AUC_min_ABcom{ee}.si.run.B = cell2mat(cellfun(@(x) sum(x),AUC_common{1, ee}.si.run.B  ,'UniformOutput',false))./time_length{1, ee}.run.B;
    %SI norun
    AUC_min_ABcom{ee}.si.norun.A = cell2mat(cellfun(@(x) sum(x),AUC_common{1, ee}.si.norun.A  ,'UniformOutput',false))./time_length{1, ee}.norun.A;
    AUC_min_ABcom{ee}.si.norun.B = cell2mat(cellfun(@(x) sum(x),AUC_common{1, ee}.si.norun.B  ,'UniformOutput',false))./time_length{1, ee}.norun.B;

    %TS run
    AUC_min_ABcom{ee}.ts.run.A = cell2mat(cellfun(@(x) sum(x),AUC_common{1, ee}.ts.run.A  ,'UniformOutput',false))./time_length{1, ee}.run.A;
    AUC_min_ABcom{ee}.ts.run.B = cell2mat(cellfun(@(x) sum(x),AUC_common{1, ee}.ts.run.B  ,'UniformOutput',false))./time_length{1, ee}.run.B;
    %TS norun
    AUC_min_ABcom{ee}.ts.norun.A = cell2mat(cellfun(@(x) sum(x),AUC_common{1, ee}.ts.norun.A  ,'UniformOutput',false))./time_length{1, ee}.norun.A;
    AUC_min_ABcom{ee}.ts.norun.B = cell2mat(cellfun(@(x) sum(x),AUC_common{1, ee}.ts.norun.B  ,'UniformOutput',false))./time_length{1, ee}.norun.B;

end

%% A and B selective neurons transient event extraction (count number of AUC events)
%from AUC for each event - significant events by definition

%get number of sig events from each class
for ee=1:size(path_dir,2)
    %A-sel events SI run
    transients_A_sel{ee}.si.run.A = cell2mat(cellfun(@(x) size(x,2),AUC_A_sel{1, ee}.si.run.A  ,'UniformOutput',false));
    transients_A_sel{ee}.si.run.B = cell2mat(cellfun(@(x) size(x,2),AUC_A_sel{1, ee}.si.run.B,'UniformOutput',false));

    %A-sel events SI norun
    transients_A_sel{ee}.si.norun.A = cell2mat(cellfun(@(x) size(x,2),AUC_A_sel{1, ee}.si.norun.A  ,'UniformOutput',false));
    transients_A_sel{ee}.si.norun.B = cell2mat(cellfun(@(x) size(x,2),AUC_A_sel{1, ee}.si.norun.B,'UniformOutput',false));

    %A-sel events TS run
    transients_A_sel{ee}.ts.run.A = cell2mat(cellfun(@(x) size(x,2),AUC_A_sel{1, ee}.ts.run.A  ,'UniformOutput',false));
    transients_A_sel{ee}.ts.run.B = cell2mat(cellfun(@(x) size(x,2),AUC_A_sel{1, ee}.ts.run.B,'UniformOutput',false));

    %A-sel events TS norun
    transients_A_sel{ee}.ts.norun.A = cell2mat(cellfun(@(x) size(x,2),AUC_A_sel{1, ee}.ts.norun.A  ,'UniformOutput',false));
    transients_A_sel{ee}.ts.norun.B = cell2mat(cellfun(@(x) size(x,2),AUC_A_sel{1, ee}.ts.norun.B,'UniformOutput',false)); 


    %B-sel events SI run
    transients_B_sel{ee}.si.run.A = cell2mat(cellfun(@(x) size(x,2),AUC_B_sel{1, ee}.si.run.A  ,'UniformOutput',false));
    transients_B_sel{ee}.si.run.B = cell2mat(cellfun(@(x) size(x,2),AUC_B_sel{1, ee}.si.run.B,'UniformOutput',false));

    %B-sel events SI norun
    transients_B_sel{ee}.si.norun.A = cell2mat(cellfun(@(x) size(x,2),AUC_B_sel{1, ee}.si.norun.A  ,'UniformOutput',false));
    transients_B_sel{ee}.si.norun.B = cell2mat(cellfun(@(x) size(x,2),AUC_B_sel{1, ee}.si.norun.B,'UniformOutput',false));

    %B-sel events TS run
    transients_B_sel{ee}.ts.run.A = cell2mat(cellfun(@(x) size(x,2),AUC_B_sel{1, ee}.ts.run.A  ,'UniformOutput',false));
    transients_B_sel{ee}.ts.run.B = cell2mat(cellfun(@(x) size(x,2),AUC_B_sel{1, ee}.ts.run.B,'UniformOutput',false));

    %B-sel events TS norun
    transients_B_sel{ee}.ts.norun.A = cell2mat(cellfun(@(x) size(x,2),AUC_B_sel{1, ee}.ts.norun.A  ,'UniformOutput',false));
    transients_B_sel{ee}.ts.norun.B = cell2mat(cellfun(@(x) size(x,2),AUC_B_sel{1, ee}.ts.norun.B,'UniformOutput',false));

    %%% Task-selective neurons by paper criteria %%%     
    %A task-sel events run
    transients_A_task_sel{ee}.run.A = cell2mat(cellfun(@(x) size(x,2),AUC_A_task_sel{1, ee}.run.A  ,'UniformOutput',false));
    transients_A_task_sel{ee}.run.B = cell2mat(cellfun(@(x) size(x,2),AUC_A_task_sel{1, ee}.run.B,'UniformOutput',false));

    %A task-sel events norun
    transients_A_task_sel{ee}.norun.A = cell2mat(cellfun(@(x) size(x,2),AUC_A_task_sel{1, ee}.norun.A ,'UniformOutput',false));
    transients_A_task_sel{ee}.norun.B = cell2mat(cellfun(@(x) size(x,2),AUC_A_task_sel{1, ee}.norun.B,'UniformOutput',false));

    %B task-sel events run
    transients_B_task_sel{ee}.run.A = cell2mat(cellfun(@(x) size(x,2),AUC_B_task_sel{1, ee}.run.A  ,'UniformOutput',false));
    transients_B_task_sel{ee}.run.B = cell2mat(cellfun(@(x) size(x,2),AUC_B_task_sel{1, ee}.run.B,'UniformOutput',false));

    %B task-sel events norun
    transients_B_task_sel{ee}.norun.A = cell2mat(cellfun(@(x) size(x,2),AUC_B_task_sel{1, ee}.norun.A ,'UniformOutput',false));
    transients_B_task_sel{ee}.norun.B = cell2mat(cellfun(@(x) size(x,2),AUC_B_task_sel{1, ee}.norun.B,'UniformOutput',false));

    %AB events run
    transients_AB{ee}.run.A = cell2mat(cellfun(@(x) size(x,2),AUC_AB{1, ee}.run.A  ,'UniformOutput',false));
    transients_AB{ee}.run.B = cell2mat(cellfun(@(x) size(x,2),AUC_AB{1, ee}.run.B,'UniformOutput',false));

    %AB events norun
    transients_AB{ee}.norun.A = cell2mat(cellfun(@(x) size(x,2),AUC_AB{1, ee}.norun.A ,'UniformOutput',false));
    transients_AB{ee}.norun.B = cell2mat(cellfun(@(x) size(x,2),AUC_AB{1, ee}.norun.B,'UniformOutput',false));

    %AB common events run SI
    transients_ABcom{ee}.si.run.A = cell2mat(cellfun(@(x) size(x,2),AUC_common{1, ee}.si.run.A  ,'UniformOutput',false));
    transients_ABcom{ee}.si.run.B = cell2mat(cellfun(@(x) size(x,2),AUC_common{1, ee}.si.run.B,'UniformOutput',false));

    %AB common events norun SI
    transients_ABcom{ee}.si.norun.A = cell2mat(cellfun(@(x) size(x,2),AUC_common{1, ee}.si.norun.A ,'UniformOutput',false));
    transients_ABcom{ee}.si.norun.B = cell2mat(cellfun(@(x) size(x,2),AUC_common{1, ee}.si.norun.B,'UniformOutput',false));

    %AB common events run TS
    transients_ABcom{ee}.ts.run.A = cell2mat(cellfun(@(x) size(x,2),AUC_common{1, ee}.ts.run.A  ,'UniformOutput',false));
    transients_ABcom{ee}.ts.run.B = cell2mat(cellfun(@(x) size(x,2),AUC_common{1, ee}.ts.run.B,'UniformOutput',false));

    %AB common events norun TS
    transients_ABcom{ee}.ts.norun.A = cell2mat(cellfun(@(x) size(x,2),AUC_common{1, ee}.ts.norun.A ,'UniformOutput',false));
    transients_ABcom{ee}.ts.norun.B = cell2mat(cellfun(@(x) size(x,2),AUC_common{1, ee}.ts.norun.B,'UniformOutput',false));
    
end


%% Calculate transients per min for all classes (start with task-selective cells in paper first)

for ee=1:size(path_dir,2)
    %transient rate A task selective neurons per paper - run
    transient_rate.A_task_sel{ee}.run.A = transients_A_task_sel{1, ee}.run.A./time_length{1, ee}.run.A;
    transient_rate.A_task_sel{ee}.run.B = transients_A_task_sel{1, ee}.run.B./time_length{1, ee}.run.B;

    %transient rate A task selective neurons per paper - norun
    transient_rate.A_task_sel{ee}.norun.A = transients_A_task_sel{1, ee}.norun.A./time_length{1, ee}.norun.A;
    transient_rate.A_task_sel{ee}.norun.B = transients_A_task_sel{1, ee}.norun.B./time_length{1, ee}.norun.B;

    %transient rate B task selective neurons per paper -run
    transient_rate.B_task_sel{ee}.run.A = transients_B_task_sel{1, ee}.run.A./time_length{1, ee}.run.A;
    transient_rate.B_task_sel{ee}.run.B = transients_B_task_sel{1, ee}.run.B./time_length{1, ee}.run.B;

    %transient rate B task selective neurons per paper - norun
    transient_rate.B_task_sel{ee}.norun.A = transients_B_task_sel{1, ee}.norun.A./time_length{1, ee}.norun.A;
    transient_rate.B_task_sel{ee}.norun.B = transients_B_task_sel{1, ee}.norun.B./time_length{1, ee}.norun.B;

    %transient rate AB task selective neurons per paper -run
    transient_rate.AB{ee}.run.A = transients_AB{1, ee}.run.A./time_length{1, ee}.run.A;
    transient_rate.AB{ee}.run.B = transients_AB{1, ee}.run.B./time_length{1, ee}.run.B;

    %transient rate AB task selective neurons per paper - norun
    transient_rate.AB{ee}.norun.A = transients_AB{1, ee}.norun.A./time_length{1, ee}.norun.A;
    transient_rate.AB{ee}.norun.B = transients_AB{1, ee}.norun.B./time_length{1, ee}.norun.B;

    %transient rate AB common neurons per paper -run - SI
    transient_rate.ABcom{ee}.si.run.A = transients_ABcom{1, ee}.si.run.A./time_length{1, ee}.run.A;
    transient_rate.ABcom{ee}.si.run.B = transients_ABcom{1, ee}.si.run.B./time_length{1, ee}.run.B;

    %transient rate AB common neurons per paper - norun - SI
    transient_rate.ABcom{ee}.si.norun.A = transients_ABcom{1, ee}.si.norun.A./time_length{1, ee}.norun.A;
    transient_rate.ABcom{ee}.si.norun.B = transients_ABcom{1, ee}.si.norun.B./time_length{1, ee}.norun.B; 

    %transient rate AB common neurons per paper -run - TS
    transient_rate.ABcom{ee}.ts.run.A = transients_ABcom{1, ee}.ts.run.A./time_length{1, ee}.run.A;
    transient_rate.ABcom{ee}.ts.run.B = transients_ABcom{1, ee}.ts.run.B./time_length{1, ee}.run.B;

    %transient rate AB common neurons per paper - norun - TS
    transient_rate.ABcom{ee}.ts.norun.A = transients_ABcom{1, ee}.ts.norun.A./time_length{1, ee}.norun.A;
    transient_rate.ABcom{ee}.ts.norun.B = transients_ABcom{1, ee}.ts.norun.B./time_length{1, ee}.norun.B;
end

%% Calculate mean transient rate for task selective neurons - matrix output for plotting
for ee=1:size(path_dir,2)
    %A first then B in matrix - A task sel run - transients/min
    transient_rate_mean.A_task_sel.run(ee,1) = mean(transient_rate.A_task_sel{ee}.run.A);
    transient_rate_mean.A_task_sel.run(ee,2) = mean(transient_rate.A_task_sel{ee}.run.B);

    %A first then B in matrix - A task sel norun - transients/min
    transient_rate_mean.A_task_sel.norun(ee,1) = mean(transient_rate.A_task_sel{ee}.norun.A);
    transient_rate_mean.A_task_sel.norun(ee,2) = mean(transient_rate.A_task_sel{ee}.norun.B);
    
    %A first then B in matrix - B task sel run - transients/min
    transient_rate_mean.B_task_sel.run(ee,1) = mean(transient_rate.B_task_sel{ee}.run.A);
    transient_rate_mean.B_task_sel.run(ee,2) = mean(transient_rate.B_task_sel{ee}.run.B);

    %A first then B in matrix - B task sel norun - transients/min
    transient_rate_mean.B_task_sel.norun(ee,1) = mean(transient_rate.B_task_sel{ee}.norun.A);
    transient_rate_mean.B_task_sel.norun(ee,2) = mean(transient_rate.B_task_sel{ee}.norun.B);

    %A first then B in matrix - AB task sel run - transients/min
    transient_rate_mean.AB.run(ee,1) = mean(transient_rate.AB{ee}.run.A);
    transient_rate_mean.AB.run(ee,2) = mean(transient_rate.AB{ee}.run.B);

    %A first then B in matrix - AB task sel norun - transients/min
    transient_rate_mean.AB.norun(ee,1) = mean(transient_rate.AB{ee}.norun.A);
    transient_rate_mean.AB.norun(ee,2) = mean(transient_rate.AB{ee}.norun.B);

    %A first then B in matrix - AB com run - SI -  transients/min
    transient_rate_mean.ABcom.si.run(ee,1) = mean(transient_rate.ABcom{ee}.si.run.A);
    transient_rate_mean.ABcom.si.run(ee,2) = mean(transient_rate.ABcom{ee}.si.run.B);

    %A first then B in matrix - AB com norun - SI - transients/min
    transient_rate_mean.ABcom.si.norun(ee,1) = mean(transient_rate.ABcom{ee}.si.norun.A);
    transient_rate_mean.ABcom.si.norun(ee,2) = mean(transient_rate.ABcom{ee}.si.norun.B);

    %A first then B in matrix - AB com run - TS -  transients/min
    transient_rate_mean.ABcom.ts.run(ee,1) = mean(transient_rate.ABcom{ee}.ts.run.A);
    transient_rate_mean.ABcom.ts.run(ee,2) = mean(transient_rate.ABcom{ee}.ts.run.B);

    %A first then B in matrix - AB com norun - TS - transients/min
    transient_rate_mean.ABcom.ts.norun(ee,1) = mean(transient_rate.ABcom{ee}.ts.norun.A);
    transient_rate_mean.ABcom.ts.norun(ee,2) = mean(transient_rate.ABcom{ee}.ts.norun.B);

end

%% Calculate mean AUC/min for common place cell controls
for ee=1:size(path_dir,2)
    %SI run
    AUC_rate_mean.ABcom.si.run(ee,1) = mean(AUC_min_ABcom{ee}.si.run.A);
    AUC_rate_mean.ABcom.si.run(ee,2) = mean(AUC_min_ABcom{ee}.si.run.B);

    %SI norun
    AUC_rate_mean.ABcom.si.norun(ee,1) = mean(AUC_min_ABcom{ee}.si.norun.A);
    AUC_rate_mean.ABcom.si.norun(ee,2) = mean(AUC_min_ABcom{ee}.si.norun.B);

    %TS run
    AUC_rate_mean.ABcom.ts.run(ee,1) = mean(AUC_min_ABcom{ee}.ts.run.A);
    AUC_rate_mean.ABcom.ts.run(ee,2) = mean(AUC_min_ABcom{ee}.ts.run.B);

    %TS norun
    AUC_rate_mean.ABcom.ts.norun(ee,1) = mean(AUC_min_ABcom{ee}.ts.norun.A);
    AUC_rate_mean.ABcom.ts.norun(ee,2) = mean(AUC_min_ABcom{ee}.ts.norun.B);
    
end

%%% TO DO: generate plots for AUC; do stats for all AUC/min and
%%% transient/min comparisons


%% Calculate the mean of means and sem for means for plotting below

%% Get mean of means and sem for bar plotting
%A ---- B
%A-sel
%B-sel
%AB

%anonymous function for calculating sem
sem = @(x) std(x,0,1)./sqrt(size(x,1));

%Mean of means - task selective neurons for figure 2
%run - transients/min with control A&B place field tuned
grouped_mean.task_sel.run = [mean(transient_rate_mean.A_task_sel.run,1);
    mean(transient_rate_mean.B_task_sel.run,1);
    mean(transient_rate_mean.AB.run,1)];

grouped_sem.task_sel.run = [sem(transient_rate_mean.A_task_sel.run);
    sem(transient_rate_mean.B_task_sel.run)
    sem(transient_rate_mean.AB.run)];
% no run
grouped_mean.task_sel.norun = [mean(transient_rate_mean.A_task_sel.norun,1);
    mean(transient_rate_mean.B_task_sel.norun,1);
    mean(transient_rate_mean.AB.norun,1)];

grouped_sem.task_sel.norun = [sem(transient_rate_mean.A_task_sel.norun);
    sem(transient_rate_mean.B_task_sel.norun)
    sem(transient_rate_mean.AB.norun)];

%common run - SI transients/min
grouped_mean.ABcom.si.run = [mean(transient_rate_mean.ABcom.si.run,1)];
grouped_sem.ABcom.si.run = [sem(transient_rate_mean.ABcom.si.run)];

% common no run SI transients/min
grouped_mean.ABcom.si.norun = [mean(transient_rate_mean.ABcom.si.norun,1)];
grouped_sem.ABcom.si.norun = [sem(transient_rate_mean.ABcom.si.norun)];

%common run - TS transients/min
grouped_mean.ABcom.ts.run = [mean(transient_rate_mean.ABcom.ts.run,1)];
grouped_sem.ABcom.ts.run = [sem(transient_rate_mean.ABcom.ts.run)];

% common no run TS transients/min
grouped_mean.ABcom.ts.norun = [mean(transient_rate_mean.ABcom.ts.norun,1)];
grouped_sem.ABcom.ts.norun = [sem(transient_rate_mean.ABcom.ts.norun)];


%% Plot transient/min for common SI/TS controls for reviewer

paper_cmap = return_paper_colormap;

fig0 = figure;
fig0.Units = 'centimeters';
fig0.Position(1) = 8;
fig0.Position(2) = 1;
fig0.Position(3) = 18;
fig0.Position(4) = 12;

%tiledLayout (replaces subplot fxn)
gridSize = [2,2];
t0 = tiledlayout(fig0,gridSize(1),gridSize(2),'TileSpacing','compact','Padding','compact','Units','centimeters');
%nest tiles within master tile

% Mean plot - dummy plot no scatter data for this
nexttile(t0,1,[1,1])
hold on;
%axis square
title('S.I. - Run');
%bar the mean for each group
b = bar(1:1,grouped_mean.ABcom.si.run,'FaceColor', 'flat');
pause(0.1)
ylabel({'Transients/min'}) 
xlim([0.5 1.5])
ylim([0 9])
%plot the sem for each mean for each group
for ib = 1:numel(b)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    if ib ==1
        xData(1,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==2
        xData(2,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==3
        xData(3,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==4
        xData(4,:) = b(ib).XData + b(ib).XOffset;
    end
    errorbar(xData(ib,:),grouped_mean.ABcom.si.run(:,ib)',grouped_sem.ABcom.si.run(:,ib),'k.','LineWidth',1)
end

%set A group bars to blue
b(1).CData(1:1,:) =  repmat(paper_cmap(1,:),1,1);
%set B group bars to red
b(2).CData(1:1,:) =  repmat(paper_cmap(2,:),1,1);
%set B group bars to red
%b(3).CData(1:3,:) =  repmat(color_mat(3,:),3,1);
%set B group bars to red
%b(4).CData(1:3,:) =  repmat(color_mat(4,:),3,1);

xticks([1 2 3]);
xticklabels({'Non-selective place cell'});

legend('A laps','B laps')

set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)

nexttile(t0,2,[1,1])
hold on;
%axis square
title('T.S. - Run');
%bar the mean for each group
b = bar(1:1,grouped_mean.ABcom.ts.run,'FaceColor', 'flat');
pause(0.1)
ylabel({'Transients/min'}) 
xlim([0.5 1.5])
ylim([0 9])
%plot the sem for each mean for each group
for ib = 1:numel(b)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    if ib ==1
        xData(1,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==2
        xData(2,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==3
        xData(3,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==4
        xData(4,:) = b(ib).XData + b(ib).XOffset;
    end
    errorbar(xData(ib,:),grouped_mean.ABcom.ts.run(:,ib)',grouped_sem.ABcom.ts.run(:,ib),'k.','LineWidth',1)
end

%set A group bars to blue
b(1).CData(1:1,:) =  repmat(paper_cmap(1,:),1,1);
%set B group bars to red
b(2).CData(1:1,:) =  repmat(paper_cmap(2,:),1,1);
%set B group bars to red
%b(3).CData(1:3,:) =  repmat(color_mat(3,:),3,1);
%set B group bars to red
%b(4).CData(1:3,:) =  repmat(color_mat(4,:),3,1);

xticks([1 2 3]);
xticklabels({'Non-selective place cell'});

legend('A laps','B laps')

set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)

nexttile(t0,3,[1,1])
hold on;
%axis square
title('S.I. - No Run');
%bar the mean for each group
b = bar(1:1,grouped_mean.ABcom.si.norun,'FaceColor', 'flat');
pause(0.1)
ylabel({'Transients/min'}) 
xlim([0.5 1.5])
ylim([0 9])
%plot the sem for each mean for each group
for ib = 1:numel(b)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    if ib ==1
        xData(1,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==2
        xData(2,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==3
        xData(3,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==4
        xData(4,:) = b(ib).XData + b(ib).XOffset;
    end
    errorbar(xData(ib,:),grouped_mean.ABcom.si.norun(:,ib)',grouped_sem.ABcom.si.norun(:,ib),'k.','LineWidth',1)
end

%set A group bars to blue
b(1).CData(1:1,:) =  repmat(paper_cmap(1,:),1,1);
%set B group bars to red
b(2).CData(1:1,:) =  repmat(paper_cmap(2,:),1,1);
%set B group bars to red
%b(3).CData(1:3,:) =  repmat(color_mat(3,:),3,1);
%set B group bars to red
%b(4).CData(1:3,:) =  repmat(color_mat(4,:),3,1);

xticks([1 2 3]);
xticklabels({'Non-selective place cell'});

legend('A laps','B laps')

set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)

nexttile(t0,4,[1,1])
hold on;
%axis square
title('T.S. - No Run');
%bar the mean for each group
b = bar(1:1,grouped_mean.ABcom.ts.norun,'FaceColor', 'flat');
pause(0.1)
ylabel({'Transients/min'}) 
xlim([0.5 1.5])
ylim([0 9])
%plot the sem for each mean for each group
for ib = 1:numel(b)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    if ib ==1
        xData(1,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==2
        xData(2,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==3
        xData(3,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==4
        xData(4,:) = b(ib).XData + b(ib).XOffset;
    end
    errorbar(xData(ib,:),grouped_mean.ABcom.ts.norun(:,ib)',grouped_sem.ABcom.ts.norun(:,ib),'k.','LineWidth',1)
end

%set A group bars to blue
b(1).CData(1:1,:) =  repmat(paper_cmap(1,:),1,1);
%set B group bars to red
b(2).CData(1:1,:) =  repmat(paper_cmap(2,:),1,1);
%set B group bars to red
%b(3).CData(1:3,:) =  repmat(color_mat(3,:),3,1);
%set B group bars to red
%b(4).CData(1:3,:) =  repmat(color_mat(4,:),3,1);

xticks([1 2 3]);
xticklabels({'Non-selective place cell'});

legend('A laps','B laps')

set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',14, ...
'FontWeight','normal', 'LineWidth', 1.5,'layer','top')


%% Plot and export that data for reviewer - temp 
%temporary plot of data for development code

%return the colormaps used in paper and apply to bars
paper_cmap = return_paper_colormap;

fig0 = figure;
fig0.Units = 'centimeters';
fig0.Position(1) = 8;
fig0.Position(2) = 1;
fig0.Position(3) = 24;
fig0.Position(4) = 12;

%tiledLayout (replaces subplot fxn)
gridSize = [2,4];
t0 = tiledlayout(fig0,gridSize(1),gridSize(2),'TileSpacing','compact','Padding','compact','Units','centimeters');
%nest tiles within master tile

% Mean plot - dummy plot no scatter data for this
nexttile(t0,1,[2,2])
hold on
axis square
xlim([0 10])
ylim([0 10])
xticks(0:2:10)
yticks(0:2:10)
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
xlabel('AUC/min - A trials')
ylabel('AUC/min - B trials')
title('A, B selective RUN')

nexttile(t0,3,[1,2])
hold on;
%axis square
title('Run');
%bar the mean for each group
b = bar(1:3,grouped_mean.task_sel.run,'FaceColor', 'flat');
pause(0.1)
ylabel('Transients/min') 
xlim([0.5 3.5])
ylim([0 9])
%plot the sem for each mean for each group
for ib = 1:numel(b)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    if ib ==1
        xData(1,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==2
        xData(2,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==3
        xData(3,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==4
        xData(4,:) = b(ib).XData + b(ib).XOffset;
    end
    errorbar(xData(ib,:),grouped_mean.task_sel.run(:,ib)',grouped_sem.task_sel.run(:,ib),'k.','LineWidth',1)
end

%set A group bars to blue
b(1).CData(1:3,:) =  repmat(paper_cmap(1,:),3,1);
%set B group bars to red
b(2).CData(1:3,:) =  repmat(paper_cmap(2,:),3,1);
%set B group bars to red
%b(3).CData(1:3,:) =  repmat(color_mat(3,:),3,1);
%set B group bars to red
%b(4).CData(1:3,:) =  repmat(color_mat(4,:),3,1);

xticks([1 2 3]);
xticklabels({'A sel.','B sel.','A&B'});

legend('A laps','B laps')

set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)

nexttile(t0,7,[1,2])
hold on;
%axis square
title('No run');
%bar the mean for each group
b2 = bar(1:3,grouped_mean.task_sel.norun,'FaceColor', 'flat');
pause(0.1)
ylabel('Transients/min') 
xlim([0.5 3.5])
ylim([0 9])
%plot the sem for each mean for each group
for ib = 1:numel(b2)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    if ib ==1
        xData(1,:) = b2(ib).XData + b2(ib).XOffset;
    elseif ib ==2
        xData(2,:) = b2(ib).XData + b2(ib).XOffset;
    elseif ib ==3
        xData(3,:) = b2(ib).XData + b2(ib).XOffset;
    elseif ib ==4
        xData(4,:) = b2(ib).XData + b2(ib).XOffset;
    end
    errorbar(xData(ib,:),grouped_mean.task_sel.norun(:,ib)',grouped_sem.task_sel.norun(:,ib),'k.','LineWidth',1)
end

%set A group bars to blue
b2(1).CData(1:3,:) =  repmat(paper_cmap(1,:),3,1);
%set B group bars to red
b2(2).CData(1:3,:) =  repmat(paper_cmap(2,:),3,1);
%set B group bars to red
%b2(3).CData(1:3,:) =  repmat(color_mat(3,:),3,1);
%set B group bars to red
%b(4).CData(1:3,:) =  repmat(color_mat(4,:),3,1);

xticks([1 2 3]);
xticklabels({'A sel.','B sel.','A&B'});
ylabel('Transients/min');
legend('A laps','B laps')



%% Source data export struct - Fig 2 data

source_data_task_sel_remap = struct();

%% Fraction of place cell (A sel, B sel, A&B remapping,neither) - bar chart and pie chart
%Figure 2C
[fraction_place_plot_data,source_data_task_sel_remap] = fraction_place_cells(path_dir,source_data_task_sel_remap);

%% A/B/A&B/N cell SI/TS score distributions plots (Figure 2C/D)
%export these values and feed into master plotter

si_ts_score_dist_data = si_ts_score_distributions(path_dir);

%% Cumulative task-selective STCs

task_sel_STC_data = cumulative_task_sel_STC(path_dir);

%% Centroid distribution for A or B selective neurons - 
%how many bins to display
%1 - 20 (every 5) relative to 100 bins
%2 - 25 (every 4) relative to 100 bins
%3 - 50 (every 2) relative to 100 bins
options.bin_choose = 2;

[centroid_dist_data,source_data_task_sel_remap] = ...
    centroid_dist_selective(path_dir,reward_zones_all_animal,options,source_data_task_sel_remap);

%% AUC/min scatterplots of A vs B neurons for each animal - DATA FOR AUC analysis here FOR REVIEWER

[AUC_data,source_data_task_sel_remap] = auc_scatterplots(path_dir,source_data_task_sel_remap);

%% PV/TC correlation

[tc_corr_sel_data,source_data_task_sel_remap] = tc_pv_correlation_task_sel(path_dir,source_data_task_sel_remap);

%% Master plotter for Figure 2
% individual place maps derive from this code for reviewer
fig2_master_plotter(fraction_place_plot_data,si_ts_score_dist_data,...
        task_sel_STC_data,centroid_dist_data,AUC_data, tc_corr_sel_data);

%% Export source data Figure 2

%navigate to matlab summary stats directory
cd('E:\matlab_data_export_paper')

save('source_data_fig2.mat', 'source_data_task_sel_remap','-v7.3');


%% Source data export struct - Fig 3

source_data_AB_remap = struct();

%% Source data export struct - supplement data/figs for Fig 2 and Fig 3

source_data_sup_2_3 = struct();

%% FIGURE 3
%% Place field A&B distance separation (global remapper distance distribution - supplement)

[combined_distances,combined_dist_metric,global_pf_dist] = place_AB_distance_separation(path_dir);

%export the lap speed data and event speed data for cumulative analysis
save(fullfile('G:\Figure_2_3_selective_remap\cumulative_data_output','pf_distances_A&B_neurons.mat'),'combined_distances', 'combined_dist_metric');

%distances between global remapping neurons for sup figure data
source_data_sup_2_3.global_pf_dist= global_pf_dist;

%% Remapping ROIs based on correlation map criteria and return quantile for common field separation (Fig 3c)
%STCs of the different categories of neurons are here (Figure 3C)
%use this as definition for the separation of fields (common and partial
%field) for partially remapping neurons
[common_cutoffs_95,remap_rate_maps] = remapping_corr(path_dir);

save(fullfile('G:\Figure_2_3_selective_remap\cumulative_data_output','common_cutoff.mat'),'common_cutoffs_95');

%% Generate mean/SEM traces for rate remapping neurons - Figure 3 rate remapping (supplement activity remap and Fig 3c data)
%also calculate A-B/A+B index and compare against common neurons
%extract supplement data
[activity_remap,remap_sup_plot_data] = rate_remap_traces(path_dir);

%sup fig 7 - cdf plot for activity index of common vs activity remappers
%and AUC/min comparison bars in fig 7c
source_data_sup_2_3.remap_sup_plot_data = remap_sup_plot_data;


%% Activity remap supplement figure

fig_sup_activity_remap_master(remap_sup_plot_data)

%% Fractional remapping using updated criteria - fraction of each class

[frac_remapping,source_data_AB_remap] =frac_remapping_neurons_corr_criteria(path_dir,source_data_AB_remap);

%% Generate scatterplot with correlation against p-value for supplement (supplement data)
%return global correlation scores for each animal/FOV

[r_global] = corr_score_scatter(path_dir);

%export this value
save(fullfile('G:\Figure_2_3_selective_remap\cumulative_data_output','r_global.mat'),'r_global');

%% Modification to original code (updated using established criteria - majority of Figure 3 data is here)
%extract scatterplot points for for global remappers from here

%return bin position of the respective reward zones
[reward_zones_all_animal,partial_idx_by_animal_zone,remap_prop_figs,global_dist_scatter,source_data_AB_remap] = ...
remapping_centroids_updated(path_dir,source_data_AB_remap);

save(fullfile('G:\Figure_2_3_selective_remap\cumulative_data_output','partial_idx_by_animal.mat'),'partial_idx_by_animal_zone');

%save the reward zone endpoints (100 bins) into a common file for use by functions
%above
save(fullfile('G:\Figure_2_3_selective_remap\cumulative_data_output','reward_zones_all_animals.mat'),'reward_zones_all_animal');


%scatterplot data for global remappers (bin position (100)
source_data_sup_2_3.global_dist_scatter = global_dist_scatter;


%% Export source data Figure 3

%navigate to matlab summary stats directory
cd('E:\matlab_data_export_paper')

save('source_data_fig3.mat', 'source_data_AB_remap','-v7.3');

%% Master plotter for Figure 3
%load frac remapping data
fig3_master_plotter(remap_rate_maps,activity_remap,frac_remapping,remap_prop_figs)

%% Event vs speed scatterplot for task selective neurons for supplement (Sup 3 - distributions) (Ex fig 4)

event_speed_plot = event_speed_scatter(path_dir);

%scatterplot data for event speed difference for sup fig 4
source_data_sup_2_3.event_speed_plot = event_speed_plot;

%% Place field analysis (width and number for selective neurons) - supplement place field data (Ex fig 4)

[pf_prop_data, pf_width_pool, pf_count] = place_field_analysis(path_dir);

%scatterplot data for place field width and # fields stats ex fig. 4
source_data_sup_2_3.pf_width_pool = pf_width_pool;
source_data_sup_2_3.pf_count = pf_count;

%% Supplement Figure 4 speed of animal in place cells and place field props

fig_sup_speed_place_field_master(event_speed_plot,pf_prop_data)

%% Supplement plotter global remappers

fig_sup_global_remap_master(global_pf_dist,global_dist_scatter,reward_zones_all_animal)

%% Speed analysis for each animal used in Figure 2/3 (mean bin speed on A/B laps) (Sup 2) 
%added speed analysis for reviewers to this section of code

%0.29-0.3 - Reward B zone start bin
%0.69-0.7 - Reward A zone start bin

pre_post_rew_zone_speed_stats = lap_speed_by_animal(path_dir);
rew_sp_exp = pre_post_rew_zone_speed_stats;

speed_zone_stats.A = [rew_sp_exp.Atrial_Azone.p, rew_sp_exp.Atrial_Azone.stats.tstat,rew_sp_exp.Atrial_Azone.stats.df];
speed_zone_stats.B = [rew_sp_exp.Btrial_Bzone.p, rew_sp_exp.Btrial_Bzone.stats.tstat,rew_sp_exp.Btrial_Bzone.stats.df];

%table for export
t_speed = array2table([speed_zone_stats.A; speed_zone_stats.B],'VariableNames',["p-val","t-stat","df"],'RowNames',["A","B"]);
writetable(t_speed,'speed_pre_post_zone_fig2_3_sup.xlsx','WriteRowNames',true) 

%save(fullfile(path_dir{1},'cumul_analysis','lap_and_event_speed.mat'),'mean_bin_speed', 'lap_bin_split','mean_event_speed');

%% Export data for supplementary/extended figure analysis

%navigate to matlab summary stats directory
cd('E:\matlab_data_export_paper')

save('source_data_sup_2_3.mat', 'source_data_sup_2_3','-v7.3');


%% Centroid difference for A&B tuned neurons and centroid diff as fxn of max bin (OLD)
%scatterplot of centroid difference as a function of center between
%centroid of max place field - not used
%centroid_difference(path_dir)

%% Remapping centroids (OLD)
%Figure 3E/F ( all Figure 3 code is here - organize this)
%options.lowPVcorr = [6 7 8];
%organize this 
%remapping_centroids(path_dir,options)

%% Fractional distribution of remapping neuron subtypes - Figure 3C (OLD)

%frac_remapping_neurons(path_dir)
