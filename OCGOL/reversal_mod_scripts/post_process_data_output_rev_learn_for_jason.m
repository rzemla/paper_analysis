%Makes script bundling the output data for rev for Jason 


%% Input animal path directories
%MR1
 path_dir{1} = {'D:\OCGOL_reversal\MR1\MR1_Random_2022_02_28-001_1',...
     'D:\OCGOL_reversal\MR1\MR1_Random_2022_03_01-001_2',...
     'D:\OCGOL_reversal\MR1\MR1_Random_2022_03_02-002_3',...
     'D:\OCGOL_reversal\MR1\MR1_RevAB_2022_03_03-001_4',...
     'D:\OCGOL_reversal\MR1\MR1_RevRandom_2022_03_04-001_5',...
     'D:\OCGOL_reversal\MR1\MR1_RevRandom_2022_03_05-001_6',...
     'D:\OCGOL_reversal\MR1\MR1_RevRandom_2022_03_08-001_7',...
     'D:\OCGOL_reversal\MR1\MR1_RevRandom_2022_03_09-001_8',...
     'D:\OCGOL_reversal\MR1\MR1_RevRandom_2022_03_10-001_9'};
%cross session directory
crossdir{1} = 'D:\OCGOL_reversal\MR1\crossSession_update';


%MR2 reversal
path_dir{2} = {'D:\OCGOL_reversal\MR2\MR2_Random_2022_03_02-001_1',...
     'D:\OCGOL_reversal\MR2\MR2_Random_2022_03_03-001_2',...
     'D:\OCGOL_reversal\MR2\MR2_Random_2022_03_04-001_3',...
     'D:\OCGOL_reversal\MR2\MR2_RevAB_2022_03_05-001_4',...
     'D:\OCGOL_reversal\MR2\MR2_RevRandom_2022_03_06-001_5',...
     'D:\OCGOL_reversal\MR2\MR2_RevRandom_2022_03_07-001_6',...
     'D:\OCGOL_reversal\MR2\MR2_RevRandom_2022_03_12-001_7',...
     'D:\OCGOL_reversal\MR2\MR2_RevRandom_2022_03_13-001_8',...
     'D:\OCGOL_reversal\MR2\MR2_RevRandom_2022_03_14-001_9'};

%cross session directory
crossdir{2} = 'D:\OCGOL_reversal\MR2\crossSession_update';

%excluded bc no run epochs, ran several laps without doing task
%'D:\OCGOL_reversal\MR4\MR4_RevRandom_2022_03_11-001_7',...
%MR4
 path_dir{3} = {'D:\OCGOL_reversal\MR4\MR4_Random_2022_03_04-001_1',...
     'D:\OCGOL_reversal\MR4\MR4_Random_2022_03_05-001_2',...
     'D:\OCGOL_reversal\MR4\MR4_Random_2022_03_06-001_3',...
     'D:\OCGOL_reversal\MR4\MR4_RevAB_2022_03_07-002_4',...
     'D:\OCGOL_reversal\MR4\MR4_RevRandom_2022_03_08-001_5',...
     'D:\OCGOL_reversal\MR4\MR4_RevRandom_2022_03_09-001_6',...
     'D:\OCGOL_reversal\MR4\MR4_RevRandom_2022_03_12-001_8',...
     'D:\OCGOL_reversal\MR4\MR4_RevRandom_2022_03_13-001_9'};
%cross session directory
crossdir{3} = 'D:\OCGOL_reversal\MR4\crossSession_update';

%% Read in registered (reg_rev_learn)

%for each animal
for ii=1:size(path_dir,2)
    reg_rev_learn{ii} = load(fullfile(crossdir{ii},'registered.mat'));
end

%% Add the filtered ROI matches to the combined registered struct

%for each animal
for ii=1:size(path_dir,2)

    %get dir path with wildcard match to .mat files
    filtered_ROI_dir_path{ii} = subdir(fullfile(crossdir{ii},'filtered_match_ROI','*.mat'));
    %load in temp var
    match_var = load(filtered_ROI_dir_path{ii}.name);
    %load in registered struct
    reg_rev_learn{ii}.registered.multi.assigned_filtered = match_var.ROI_assign_multi_filtered;

end


%% Modify registered struct to remove MR 4, session 7 matches
       
    %assign registered temp variable to MR4
    registered = reg_rev_learn{3}.registered;
    %multi-matching session parsing
    registered.multi.assigned(:,7) = [];
    registered.multi.assigned_all(:,7) = [];
    registered.multi.assigned_filtered(:,7) = [];
    registered.multi.ROI_outlines(:,7) = [];
    registered.multi.ROI_zooms(:,7) = [];

    temp = registered.multi.assign_cell;
    temp(7,:) = [];
    temp(:,7) = [];
    registered.multi.assign_cell = temp;

    temp = registered.session.assign_cell;
    temp(7,:) = [];
    temp(:,7) = [];
    registered.session.assign_cell = temp;

%reassign registered struct back to MR4
reg_rev_learn{3}.registered = registered;

%% Read in vars for short_term_rev data struct

%for each animal
for ii=1:size(path_dir,2)
    %import PV TC corr data for each animal
    temp = load(fullfile(crossdir{ii},'PV_TC_corr.mat'));
    PV_TC_corr(ii) = temp;
    
    %performance data
    perf{ii} = load(fullfile(crossdir{ii},'ses_perf.mat'));
    
    %tuning curve matching data
    TC_corr_match{ii} = load(fullfile(crossdir{ii},'tc_corr_match.mat'));
    
    %fraction tuned 
    tuned_frac{ii} = load(fullfile(crossdir{ii},'tuned_fractions.mat'));
    
    %tuned logicals SI/TS
    tuned_log{ii} = load(fullfile(crossdir{ii},'tuned_logicals.mat'));
    
    %recurrence
    recurr{ii} = load(fullfile(crossdir{ii},'recurrence.mat'));
    
    %place field vectors
    pf_vector_max{ii} = load(fullfile(crossdir{ii},'pf_vector_max.mat'));
end

%% Assign into one struct and save output mat file
%merge into one struct
short_term_rev.perf = perf;
short_term_rev.recurr = recurr;
short_term_rev.pf_vector_max = pf_vector_max;
short_term_rev.tuned_frac = tuned_frac;
short_term_rev.tuned_log = tuned_log;
short_term_rev.PV_TC_corr = PV_TC_corr;
short_term_rev.TC_corr_match = TC_corr_match;

%save stuct as mat file
save('D:\OCGOL_reversal\tuned_match_ROI_rev_data.mat', "short_term_rev");
