function [path_dir_learn,path_dir_recall,crossdir_learn,crossdir_recall] = figure4_input_data()

%% List session directories for each animal (cells) - learning
 path_dir_learn{1} = {'G:\OCGOL_learning_short_term\I56_RTLS\I56_RLTS_5AB_041019_1',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_5AB_041119_2',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_3A3B_041219_3',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_3A3B_041319_4',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_no_punish_041519_5',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_no_punish_041619_6'};
% %cross session directory
 crossdir_learn{1} = 'G:\OCGOL_learning_short_term\I56_RTLS\crossSession';

% %I57_RTLS
 path_dir_learn{2} = {'G:\OCGOL_learning_short_term\I57_RTLS\I57_RLTS_5AB_041019_1',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_5AB_041119_2',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_3A3B_041219_3',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_1A1B_041319_4',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_no_punish_041519_5',...
     'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_no_punish_041619_6'};
     %'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_punish_041719_7'};
% %cross session directory
 crossdir_learn{2} = 'G:\OCGOL_learning_short_term\I57_RTLS\crossSession';

%I57_LT
 path_dir_learn{3} = {'G:\OCGOL_learning_short_term\I57_LT\I57_LT_5A5B_041619_1',...
     'G:\OCGOL_learning_short_term\I57_LT\I57_LT_5A5B_041719_2',...
     'G:\OCGOL_learning_short_term\I57_LT\I57_LT_3A3B_041819_3',...
     'G:\OCGOL_learning_short_term\I57_LT\I57_LT_3A3B_041919_4',...
     'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_no_punish_042019_5',...
     'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_no_punish_042119_6'};

 crossdir_learn{3} = 'G:\OCGOL_learning_short_term\I57_LT\crossSession';
 
%% List session directories for each animal (cells) - recall

path_dir_recall{1} = {'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d1_032118_1',...
     'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d2_032218_2',...
     'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d3_032318_3',...
     'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d6_032618_4_2',...
     'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d7_032718_5',...
     'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d8_032818_6',...
     'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d9_032918_7'};
%cross session directory
crossdir_recall{1} = 'G:\OCGOL_stability_recall\I42R_1\crossSession';
% 
path_dir_recall{2} = {'G:\OCGOL_stability_recall\I46\I46_AB_d1_062018_1',...
     'G:\OCGOL_stability_recall\I46\I46_AB_d2_062118_2',...
     'G:\OCGOL_stability_recall\I46\I46_AB_d3_062218_3',...
     'G:\OCGOL_stability_recall\I46\I46_AB_d6_062518_4',...
     'G:\OCGOL_stability_recall\I46\I46_AB_d7_062618_5',...
     'G:\OCGOL_stability_recall\I46\I46_AB_d8_062718_6',...
     'G:\OCGOL_stability_recall\I46\I46_AB_d9_062818_7'};
% %cross session directory
crossdir_recall{2} = 'G:\OCGOL_stability_recall\I46\crossSession';

path_dir_recall{3} = {'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d1_062018_1',...
     'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d2_062118_2',...
     'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d3_062218_3',...
     'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d6_062518_4',...
     'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d7_062618_5',...
     'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d8_062718_6',...
     'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d9_062818_7'};
% %cross session directory
crossdir_recall{3} = 'G:\OCGOL_stability_recall\I45_RT\crossSession';

 path_dir_recall{4} = {'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d1_032118_1',...
     'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d2_032218_2',...
     'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d3_032318_3',...
     'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d6_032618_4',...
     'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d7_032718_5',...
     'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d8_032818_6',...
     'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d9_032918_7'};
 
% %cross session directory
crossdir_recall{4} = 'G:\OCGOL_stability_recall\I42L_1\crossSession';
 
end
