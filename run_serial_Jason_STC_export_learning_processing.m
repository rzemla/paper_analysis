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

 crossdir{1} = 'G:\OCGOL_learning_short_term\I56_RTLS\crossSession';
 
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
crossdir{2} = 'G:\OCGOL_learning_short_term\I57_RTLS\crossSession';

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
crossdir{3} = 'G:\OCGOL_learning_short_term\I57_LT\crossSession';


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
crossdir{4} = 'E:\OCGOL_learning_short_term\I58_RT\crossSession';

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
crossdir{5} = 'E:\OCGOL_learning_short_term\I58_LT\crossSession';

%ANIMAL #6
%I58 RTLP
path_dir{6} = {'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_5A5B_080419_1',...
     'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_5A5B_080519_2',...
     'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_3A3B_080619_3',...
     'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_3A3B_080719_4',...
     'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_randAB_no_punish_080819_5',...
     'E:\OCGOL_learning_short_term\I58_RTLP\I58_RTLP_randAB_no_punish_080919_6'};

%cross session directory
crossdir{6} = 'E:\OCGOL_learning_short_term\I58_RTLP\crossSession';
 

%% Run the analysis for STC export for splitter cell analysis

for ii=[1 3 4 5 6]
    disp(['Running: ', num2str(ii)])
    OCGOL_learning_remapping_Jason_export_shortened_V0(path_dir{ii},crossdir{ii})
end

%% Run the export of correct vs incorrect trial analysis (does not include random selection of laps)

for ii=[1 2 3 4 5 6]
    disp(['Running: ', num2str(ii)])
    extract_STCs_for_split_session_analysis(path_dir{ii},crossdir{ii})
end




