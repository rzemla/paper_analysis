%% Define the directories for processing
%ANIMAL #1
%I56_RTLS
path_dir{1} = {'G:\OCGOL_learning_short_term\I56_RTLS\I56_RLTS_5AB_041019_1',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_5AB_041119_2',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_3A3B_041219_3',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_3A3B_041319_4',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_no_punish_041519_5',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_no_punish_041619_6',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_punish_041719_7'};

crossdir{1} = 'G:\OCGOL_learning_short_term\I56_RTLS\crossSession_update';
 
%ANIMAL #2
%I57_RTLS
path_dir{2} = {'G:\OCGOL_learning_short_term\I57_RTLS\I57_RLTS_5AB_041019_1',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_5AB_041119_2',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_3A3B_041219_3',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_1A1B_041319_4',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_no_punish_041519_5',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_no_punish_041619_6',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_punish_041719_7'};
 
%cross session directory
crossdir{2} = 'G:\OCGOL_learning_short_term\I57_RTLS\crossSession_update';

%ANIMAL #3
%I57_LT
path_dir{3} = {'G:\OCGOL_learning_short_term\I57_LT\I57_LT_5A5B_041619_1',...
     'G:\OCGOL_learning_short_term\I57_LT\I57_LT_5A5B_041719_2',...
     'G:\OCGOL_learning_short_term\I57_LT\I57_LT_3A3B_041819_3',...
     'G:\OCGOL_learning_short_term\I57_LT\I57_LT_3A3B_041919_4',...
     'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_no_punish_042019_5',...
     'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_no_punish_042119_6',...
     'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_punish_042219_7',...
     'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_punish_042319_8'};
%  
crossdir{3} = 'G:\OCGOL_learning_short_term\I57_LT\crossSession_update';


%ANIMAL #4
%I58 RT
%input directories to matching function
path_dir{4} = {'E:\OCGOL_learning_short_term\I58_RT\I58_RT_5A5B_073019_1',...
     'E:\OCGOL_learning_short_term\I58_RT\I58_RT_5A5B_073119_2',...
     'E:\OCGOL_learning_short_term\I58_RT\I58_RT_3A3B_080119_3',...
     'E:\OCGOL_learning_short_term\I58_RT\I58_RT_3A3B_080219_4',...
     'E:\OCGOL_learning_short_term\I58_RT\I58_RT_randAB_no_punish_080319_5',...
     'E:\OCGOL_learning_short_term\I58_RT\I58_RT_randAB_no_punish_080419_6',...
     'E:\OCGOL_learning_short_term\I58_RT\I58_RT_ABrand_punish_080519_7',...
     'E:\OCGOL_learning_short_term\I58_RT\I58_RT_randAB_punish_080619_8',...
     'E:\OCGOL_learning_short_term\I58_RT\I58_RT_randAB_punish_080719_9'};
% %cross session directory
crossdir{4} = 'E:\OCGOL_learning_short_term\I58_RT\crossSession_update';

%ANIMAL #5
%I58 LT
path_dir{5} = {'E:\OCGOL_learning_short_term\I58_LT\I58_LT_5A5B_080419_1',...
     'E:\OCGOL_learning_short_term\I58_LT\I58_LT_5A5B_080519_2',...
     'E:\OCGOL_learning_short_term\I58_LT\I58_LT_3A3B_080619_3',...
     'E:\OCGOL_learning_short_term\I58_LT\I58_LT_3A3B_080719_4',...
     'E:\OCGOL_learning_short_term\I58_LT\I58_LT_randAB_no_punish_080819_5',...
     'E:\OCGOL_learning_short_term\I58_LT\I58_LT_randAB_no_punish_080919_6',...
     'E:\OCGOL_learning_short_term\I58_LT\I58_LT_randAB_punish_081119_7'};
%cross session directory
crossdir{5} = 'E:\OCGOL_learning_short_term\I58_LT\crossSession_update';

%ANIMAL #6
%I58 RTLP
path_dir{6} = {'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_5A5B_080419_1',...
     'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_5A5B_080519_2',...
     'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_3A3B_080619_3',...
     'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_3A3B_080719_4',...
     'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_randAB_no_punish_080819_5',...
     'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_randAB_no_punish_080919_6'};

%cross session directory
crossdir{6} = 'E:\OCGOL_learning_short_term\I58_RTLP\crossSession_update';

%% Import variables and define options (set options to use to run global script) - LEARNING PARAMETERS

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

%change the the parameters below for whether it's a learning or recall
%experiment

%all A and B trials used for learning (sessions determined below)
%use [1 2] for recall
options.selectTrial = [4 5];

%for use in global workspace
%selectTrial = options.selectTrial;

%flag to all A or B trial or only correct A or B trials
%all correct = 0 ==> uses trials 4,5 (set for learning data)
%all correct = 1 ==> uses trials 1,2 (set for recall data)
options.allCorrect = 0;

%this is just for controlling the display of labels across sessions/days
%doesn't control anything else
options.learning_data = 1;

%% Run processing of data for generation of Figure 4

for ii=[5 6]
    disp(['Running: ', num2str(ii)])
    OCGOL_global_remapping_V0_update(path_dir{ii},crossdir{ii},options)
end







