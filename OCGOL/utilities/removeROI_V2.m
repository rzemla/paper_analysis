function [C_df_filtered,F_vars_filtered] = removeROI_V2(C_df, F_vars, directory_name)
%removes ROIs from F_vars struct and C_df variables

%% Load in the logical vector with selected ROIs

%get the name of the .mat file
removedROI_dir = dir([directory_name,'\','removedROI\','*selectedROIs.mat']);

%load the ROI selected logical variable
load(fullfile(removedROI_dir.folder, removedROI_dir.name),'compSelect');

%% Remove the ROIs based on logical vector in 

%C_df input variable to scripts
C_df_filtered = C_df(:,compSelect);

%raw F
F_vars_filtered.F = F_vars.F(compSelect,:);
%dF/F
F_vars_filtered.F_dff = F_vars.F_dff(compSelect,:);
%denoised dF/F
F_vars_filtered.F_dff_exp = F_vars.F_dff_exp(compSelect,:);
%F- F0_baseline
F_vars_filtered.Fd = F_vars.Fd(compSelect,:);

%background/baseline traces
F_vars_filtered.F0 = F_vars.F0(compSelect,:);
F_vars_filtered.F0_background = F_vars.F0_background(compSelect,:);
F_vars_filtered.F0_baseline = F_vars.F0_baseline(compSelect,:);

end

