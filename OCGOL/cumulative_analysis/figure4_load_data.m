function [CNMF_learn,reg_learn,reg_recall,reg_recall_long,session_vars_learn,...
        session_vars_recall,session_vars_recall_long,...
        short_term_learn,short_term_recall,long_term_recall] = figure4_load_data(options)

      %% Skip loading of SCE data for now
      %assign dummy outputs for now
      SCE_recall = [];
      SCE_learning = [];
      
      
%% Load in pre-defined experiments directories for all animals -learning and recall

%return struct with paths of individual sessions as well as crossSession
%directories for each animal 
[path_dir_learn,path_dir_recall,crossdir_learn,crossdir_recall,...
    path_dir_recall_long, crossdir_recall_long] = figure4_input_data();

%% Get folder directories for all sessions for each animal
%get mat directories in each output folder
%for each learning animal
for ii=1:size(path_dir_learn,2)
    %get matfile names for each session
    %matfiles{ii} = dir([path_dir{ii},'\output','\*.mat']);
    %input directories (.mat CNMF, .XML, .CSV)
    %for each session for that animal
    for ss =1:size(path_dir_learn{ii},2)
        %get input dir names learn
        inputfiles_learn{ii}{ss} = dir([path_dir_learn{ii}{ss},'\input','\*.mat']);
        %get output dir names learn
        outputfiles_learn{ii}{ss} = dir([path_dir_learn{ii}{ss},'\output','\*.mat']);
        %removed ROI directory
        removeROIfiles_learn{ii}{ss} = dir([path_dir_learn{ii}{ss},'\removedROI','\*.mat']);
    end
end


if options.skipSesVars == 0
    
    %learn session variables
    %for each animal
    for ii=1:size(path_dir_learn,2)
        disp(ii)
        %load in session varaibles
        for ss =1:size(path_dir_learn{ii},2)
            disp(ss)
            %decide which variables here do not need to be loaded
            session_vars_learn{ii}{ss} = load(fullfile(outputfiles_learn{ii}{ss}.folder, outputfiles_learn{ii}{ss}.name),...
                'Behavior_split', 'Imaging_split');
        end
    end
else %set blank variable for session vars
    session_vars_learn = [];
end

%for each recall animal
for ii=1:size(path_dir_recall,2)
    %get matfile names for each session
    %matfiles{ii} = dir([path_dir{ii},'\output','\*.mat']);
    %input directories (.mat CNMF, .XML, .CSV)
    %for each session for that animal
    for ss =1:size(path_dir_recall{ii},2)
        %input files for recall 
        inputfiles_recall{ii}{ss} = dir([path_dir_recall{ii}{ss},'\input','\*.mat']);
        %output files for recall
        outputfiles_recall{ii}{ss} = dir([path_dir_recall{ii}{ss},'\output','\*.mat']);
        %removed ROI directory
        removeROIfiles_recall{ii}{ss} = dir([path_dir_recall{ii}{ss},'\removedROI','\*.mat']);
    end
end

if options.skipSesVars == 0

%recall session variables
%for each animal
for ii=1:size(path_dir_recall,2)
    disp(ii)
    %load in session varaibles
    for ss =1:size(path_dir_recall{ii},2)
        disp(ss)
        %decide which variables here do not need to be loaded
        session_vars_recall{ii}{ss} = load(fullfile(outputfiles_recall{ii}{ss}.folder, outputfiles_recall{ii}{ss}.name),...
            'Behavior_split', 'Imaging_split');
    end
end

else %set blank variable for session vars
    session_vars_recall = [];

end

%% Long term session load variables
disp('Load in long term animals');
%for each recall animal
for ii=1:size(path_dir_recall_long,2)
    %get matfile names for each session
    %matfiles{ii} = dir([path_dir{ii},'\output','\*.mat']);
    %input directories (.mat CNMF, .XML, .CSV)
    %for each session for that animal
    for ss =1:size(path_dir_recall_long{ii},2)
        %input files for recall 
        inputfiles_recall_long{ii}{ss} = dir([path_dir_recall_long{ii}{ss},'\input','\*.mat']);
        %output files for recall
        outputfiles_recall_long{ii}{ss} = dir([path_dir_recall_long{ii}{ss},'\output','\*.mat']);
        %removed ROI directory
        removeROIfiles_recall_long{ii}{ss} = dir([path_dir_recall_long{ii}{ss},'\removedROI','\*.mat']);
    end
end

if options.skipSesVars == 0
    
    %recall session variables
    %for each animal
    for ii=1:size(path_dir_recall_long,2)
        disp(ii)
        %load in session varaibles
        for ss =1:size(path_dir_recall_long{ii},2)
            disp(ss)
            %decide which variables here do not need to be loaded
            session_vars_recall_long{ii}{ss} = load(fullfile(outputfiles_recall_long{ii}{ss}.folder, outputfiles_recall_long{ii}{ss}.name),...
                'Behavior_split', 'Imaging_split');
        end
    end
    
else %set blank variable for session vars
    session_vars_recall_long = [];
    
end

%% CNMF component outlines and templates - just learning for Figure 4 display (match components FOV)
for ii = 1:size(path_dir_learn,2)
    %load in relevant CNMF variables for visualization purposes
    %CNMF_vars{ii} = load(fullfile(inputfiles{ii}.folder,inputfiles{ii}.name),'A_keep','C_keep','Cn','Coor_kp','expDffMedZeroed','dims');
    for ss =1:size(path_dir_learn{ii},2)
        CNMF_vars_learn{ii}{ss} = load(fullfile(inputfiles_learn{ii}{ss}.folder,inputfiles_learn{ii}{ss}.name),'Coor_kp');
        %load in template (stack average)
        templates_learn{ii}{ss} = load(fullfile(path_dir_learn{ii}{ss},'template_nr.mat'),'template');
        %load in logical vector with selected and rejected ROIs
        removeROI_learn{ii}{ss} = load(fullfile(removeROIfiles_learn{ii}{ss}.folder,removeROIfiles_learn{ii}{ss}.name),'compSelect');
    end
end

%place variables into struct
CNMF_learn.CNMF_vars_learn = CNMF_vars_learn;
CNMF_learn.templates_learn = templates_learn;
CNMF_learn.removeROI_learn = removeROI_learn;

%% Load registered struct for each animal
%CNMF matched and manually filtered
%learning
disp('Load in registration data for SHORT term learning animals')
for ii = 1:size(path_dir_learn,2)
    reg_learn{ii} = load(fullfile(crossdir_learn{ii},'registered.mat'),'registered');
    %get dir path with wildcard match to .mat files
    filtered_ROI_dir_path_learn{ii} = subdir(fullfile(crossdir_learn{ii},'filtered_match_ROI','*.mat'));
    %load in temp var
    match_var_learn{ii} = load(filtered_ROI_dir_path_learn{ii}.name);
    %load in registered struct
    reg_learn{ii}.registered.multi.assigned_filtered = match_var_learn{ii}.ROI_assign_multi_filtered;
    disp('Learn Loaded.');
end

%short term recall
disp('Load in registration data for SHORT term recall animals')
%recall (short term 9 days)
for ii = 1:size(path_dir_recall,2)
    reg_recall{ii} = load(fullfile(crossdir_recall{ii},'registered.mat'),'registered');
        %get dir path with wildcard match to .mat files
    filtered_ROI_dir_path_recall{ii} = subdir(fullfile(crossdir_recall{ii},'filtered_match_ROI','*.mat'));
    %load in temp var
    match_var_recall{ii} = load(filtered_ROI_dir_path_recall{ii}.name);
    %load in registered struct
    reg_recall{ii}.registered.multi.assigned_filtered = match_var_recall{ii}.ROI_assign_multi_filtered;
    disp('Recall Loaded.');
end


%short term recall
disp('Load in registration data for LONG term recall animals')
%recall (short term 9 days)
for ii = 1:size(path_dir_recall_long,2)
    reg_recall{ii} = load(fullfile(crossdir_recall_long{ii},'registered.mat'),'registered');
        %get dir path with wildcard match to .mat files
    filtered_ROI_dir_path_recall_long{ii} = subdir(fullfile(crossdir_recall_long{ii},'filtered_match_ROI','*.mat'));
    %load in temp var
    match_var_recall_long{ii} = load(filtered_ROI_dir_path_recall_long{ii}.name);
    %load in registered struct
    reg_recall_long{ii}.registered.multi.assigned_filtered = match_var_recall_long{ii}.ROI_assign_multi_filtered;
    disp('Recall Loaded.');
end


%% read in short term recall data
for ss=1:size(crossdir_recall,2)
    %TC/PV correlations
    PV_TC_corr_recall(ss) = load(fullfile(crossdir_recall{ss},'PV_TC_corr.mat'));
    %performance data
    perf_recall{ss} = load(fullfile(crossdir_recall{ss},'ses_perf.mat'));
    %load SCE related data
    %SCE_recall{ss} = load(fullfile(crossdir_recall{ss},'SCE.mat'));
    %load Behavior and Behavior_split and Imaging struct
    
    %TC correlation for TS-event tuned/matched ROIs
    TC_corr_match_recall{ss} = load(fullfile(crossdir_recall{ss},'tc_corr_match.mat'));
    
    %fractions of neurons tuned in each cateogory across learning
    tuned_frac_recall{ss} = load(fullfile(crossdir_recall{ss},'tuned_fractions.mat'));

    %tuned logicals - which neurons are tuned by each category on given
    %session
    tuned_log_recall{ss} = load(fullfile(crossdir_recall{ss},'tuned_logicals.mat'));
    %pf vector max - vector of maximum firing place field for each neuron
    pf_vector_max_recall{ss} = load(fullfile(crossdir_recall{ss},'pf_vector_max.mat'));
    
    %recurrence data
    recur_recall{ss} = load(fullfile(crossdir_recall{ss},'recurrence.mat'));
    
end

%% read in long term recall data
for ss=1:size(crossdir_recall_long,2)
    %TC/PV correlations
    PV_TC_corr_recall_long(ss) = load(fullfile(crossdir_recall_long{ss},'PV_TC_corr.mat'));
    %performance data
    perf_recall_long{ss} = load(fullfile(crossdir_recall_long{ss},'ses_perf.mat'));
    %load SCE related data
    %SCE_recall{ss} = load(fullfile(crossdir_recall{ss},'SCE.mat'));
    %load Behavior and Behavior_split and Imaging struct
    
    %TC correlation for TS-event tuned/matched ROIs
    TC_corr_match_recall_long{ss} = load(fullfile(crossdir_recall_long{ss},'tc_corr_match.mat'));
    
    %fractions of neurons tuned in each cateogory across learning
    tuned_frac_recall_long{ss} = load(fullfile(crossdir_recall_long{ss},'tuned_fractions.mat'));

    %tuned logicals - which neurons are tuned by each category on given
    %session
    tuned_log_recall_long{ss} = load(fullfile(crossdir_recall_long{ss},'tuned_logicals.mat'));
    %pf vector max - vector of maximum firing place field for each neuron
    pf_vector_max_recall_long{ss} = load(fullfile(crossdir_recall_long{ss},'pf_vector_max.mat'));
    
    %recurrence data
    recur_recall_long{ss} = load(fullfile(crossdir_recall_long{ss},'recurrence.mat'));
    
end


%% read in short term learn data
for ss=1:size(crossdir_learn,2)
    %TC/PV correlations
    PV_TC_corr_learning(ss) = load(fullfile(crossdir_learn{ss},'PV_TC_corr.mat'));
    %performance data
    perf_learning{ss} = load(fullfile(crossdir_learn{ss},'ses_perf.mat'));
    %load SCE related data
    %SCE_learning{ss} = load(fullfile(crossdir_learn{ss},'SCE.mat'));
    
    %TC correlation for TS-event tuned/matched ROIs
    TC_corr_match_learning{ss} = load(fullfile(crossdir_learn{ss},'tc_corr_match.mat'));
    
    %fractions of neurons tuned in each cateogory across learning
    tuned_frac_learning{ss} = load(fullfile(crossdir_learn{ss},'tuned_fractions.mat'));
    
    %tuned logicals - which neurons are tuned by each category on given
    %session
    tuned_log_learning{ss} = load(fullfile(crossdir_learn{ss},'tuned_logicals.mat'));
    %pf vector max - vector of maximum firing place field for each neuron
    pf_vector_max_learning{ss} = load(fullfile(crossdir_learn{ss},'pf_vector_max.mat'));
    
    %recurrence data
    recur_learn{ss} = load(fullfile(crossdir_learn{ss},'recurrence.mat'));
    
    
end

%% Bundle learning and recall data into structs for export and analysis
%short term learning struct
short_term_learn.PV_TC_corr = PV_TC_corr_learning;
short_term_learn.perf = perf_learning;
short_term_learn.TC_corr_match = TC_corr_match_learning;
short_term_learn.tuned_frac = tuned_frac_learning;
%logical of tuned neurons (A/B/AB for TS/SI with min 5 events)
short_term_learn.tuned_log = tuned_log_learning;
%max tuning vectors
short_term_learn.pf_vector_max = pf_vector_max_learning;
%recurrence
short_term_learn.recurr = recur_learn;


%short term recall struct
short_term_recall.PV_TC_corr = PV_TC_corr_recall;
short_term_recall.perf = perf_recall;
short_term_recall.TC_corr_match = TC_corr_match_recall;
short_term_recall.tuned_frac = tuned_frac_recall;
%logical of tuned neurons (A/B/AB for TS/SI with min 5 events)
short_term_recall.tuned_log = tuned_log_recall;
%max tuning vectors
short_term_recall.pf_vector_max = pf_vector_max_recall;
%recurrence
short_term_recall.recurr = recur_recall;

%long term recall struct
long_term_recall.PV_TC_corr = PV_TC_corr_recall_long;
long_term_recall.perf = perf_recall_long;
long_term_recall.TC_corr_match = TC_corr_match_recall_long;
long_term_recall.tuned_frac = tuned_frac_recall_long;
%logical of tuned neurons (A/B/AB for TS/SI with min 5 events)
long_term_recall.tuned_log = tuned_log_recall_long;
%max tuning vectors
long_term_recall.pf_vector_max = pf_vector_max_recall_long;
%recurrence
long_term_recall.recurr = recur_recall_long;



% learn_comb_data.tuned_log_learning = tuned_log_learning;
% recall_comb_data.tuned_log_recall = tuned_log_recall;
% 
% %max tuning vectors
% learn_comb_data.pf_vector_max_learning = pf_vector_max_learning;
% recall_comb_data.pf_vector_max_recall = pf_vector_max_recall;


end

