%% Define the directories for processing

path_dir{1} = {'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1'};
path_dir{2} = {'G:\Figure_2_3_selective_remap\I42R_AB_d1_032118_1'};
path_dir{3} = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1'};
path_dir{4} = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_2'};
path_dir{5} = {'G:\Figure_2_3_selective_remap\I53LT_AB_sal_113018_1'}; %place field finder problem - adjust
path_dir{6} = {'G:\Figure_2_3_selective_remap\I56_RTLS_AB_prePost_sal_042419_1'}; %place field finder problem - adjust
path_dir{7} = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_113018_1'};
path_dir{8} = {'G:\Figure_2_3_selective_remap\I57_RTLS_AB_prePost_792_042519_1'};
path_dir{9} = {'G:\Figure_2_3_selective_remap\I45_RT_AB_d1_062018_1'};
path_dir{10} = {'G:\Figure_2_3_selective_remap\I46_AB_d1_062018_1'};
path_dir{11} = {'G:\Figure_2_3_selective_remap\I57_LT_ABrand_no_punish_042119_1'};

%% Run the analysis
for ii=1:size(path_dir,2)
    disp(['Running: ', num2str(ii)])
    figure_2_3_analysis(path_dir{ii})
    
end