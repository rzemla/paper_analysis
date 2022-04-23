%% Extended mutli session match GUI - 9 sessions vs previous 7 session select
clear
%% List of datasets

%input directories to matching function
%I56 RTLS
%  path_dir = {'G:\OCGOL_learning_short_term\I56_RTLS\I56_RLTS_5AB_041019_1',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_5AB_041119_2',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_3A3B_041219_3',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_3A3B_041319_4',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_no_punish_041519_5',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_no_punish_041619_6',...
%      'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_punish_041719_7'};
% %cross session directory
% crossdir = 'G:\OCGOL_learning_short_term\I56_RTLS\crossSession';

%  path_dir = {'G:\OCGOL_learning_short_term\I57_RTLS\I57_RLTS_5AB_041019_1',...
%      'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_5AB_041119_2',...
%      'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_3A3B_041219_3',...
%      'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_1A1B_041319_4',...
%      'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_no_punish_041519_5',...
%      'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_no_punish_041619_6',...
%      'G:\OCGOL_learning_short_term\I57_RTLS\I57_RTLS_ABrand_punish_041719_7'};
% %cross session directory
% crossdir = 'G:\OCGOL_learning_short_term\I57_RTLS\crossSession';

%I57_LT
%  path_dir = {'G:\OCGOL_learning_short_term\I57_LT\I57_LT_5A5B_041619_1',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_5A5B_041719_2',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_3A3B_041819_3',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_3A3B_041919_4',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_no_punish_042019_5',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_no_punish_042119_6',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_punish_042219_7',...
%      'G:\OCGOL_learning_short_term\I57_LT\I57_LT_ABrand_punish_042319_8'};
% %cross session directory
% crossdir = 'G:\OCGOL_learning_short_term\I57_LT\crossSession';

%  path_dir = {'G:\OCGOL_stability_recall\I47_LP\I47_LP_AB_d1_062018_1',...
%      'G:\OCGOL_stability_recall\I47_LP\I47_LP_AB_d2_062118_2',...
%      'G:\OCGOL_stability_recall\I47_LP\I47_LP_AB_d3_062218_3',...
%      'G:\OCGOL_stability_recall\I47_LP\I47_LP_AB_d6_062518_4',...
%      'G:\OCGOL_stability_recall\I47_LP\I47_LP_AB_d7_062618_5',...
%      'G:\OCGOL_stability_recall\I47_LP\I47_LP_AB_d8_062718_6',...
%      'G:\OCGOL_stability_recall\I47_LP\I47_LP_AB_d9_062818_7'};
% %cross session directory
% crossdir = 'G:\OCGOL_stability_recall\I47_LP\crossSession';

%  path_dir = {'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d1_032118_1',...
%      'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d2_032218_2',...
%      'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d3_032318_3',...
%      'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d6_032618_4',...
%      'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d7_032718_5',...
%      'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d8_032818_6',...
%      'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d9_032918_7'};
% %cross session directory
% crossdir = 'G:\OCGOL_stability_recall\I42L_1\crossSession';

%  path_dir = {'G:\OCGOL_stability_recall\I46\I46_AB_d1_062018_1',...
%      'G:\OCGOL_stability_recall\I46\I46_AB_d2_062118_2',...
%      'G:\OCGOL_stability_recall\I46\I46_AB_d3_062218_3',...
%      'G:\OCGOL_stability_recall\I46\I46_AB_d6_062518_4',...
%      'G:\OCGOL_stability_recall\I46\I46_AB_d7_062618_5',...
%      'G:\OCGOL_stability_recall\I46\I46_AB_d8_062718_6',...
%      'G:\OCGOL_stability_recall\I46\I46_AB_d9_062818_7'};
% %cross session directory
% crossdir = 'G:\OCGOL_stability_recall\I46\crossSession';
% 
%  path_dir = {'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d1_032118_1',...
%      'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d2_032218_2',...
%      'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d3_032318_3',...
%      'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d6_032618_4_2',...
%      'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d7_032718_5',...
%      'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d8_032818_6',...
%      'G:\OCGOL_stability_recall\I42R_1\I42R_AB_d9_032918_7'};
% %cross session directory
% crossdir = 'G:\OCGOL_stability_recall\I42R_1\crossSession';

%  path_dir = {'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d1_062018_1',...
%      'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d2_062118_2',...
%      'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d3_062218_3_2',...
%      'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d6_062518_4',...
%      'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d7_062618_5',...
%      'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d8_062718_6',...
%      'G:\OCGOL_stability_recall\I45_RT\I45_RT_AB_d9_062818_7'};
% %cross session directory
% crossdir = 'G:\OCGOL_stability_recall\I45_RT\crossSession';

%I58_RT
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

%I47 LP Long term
%  path_dir = {'D:\OCGOL_learning_long_term\I47_LP\I47_LP_AB_d1_062018_1',...
%      'D:\OCGOL_learning_long_term\I47_LP\I47_LP_AB_d6_062518_2',...
%      'D:\OCGOL_learning_long_term\I47_LP\I47_LP_AB_d16_070518_3',...
%      'D:\OCGOL_learning_long_term\I47_LP\I47_LP_AB_d20_070918_4',...
%      'D:\OCGOL_learning_long_term\I47_LP\I47_LP_AB_d25_071318_5',...
%      'D:\OCGOL_learning_long_term\I47_LP\I47_LP_AB_d30_071818_6'};
% %cross session directory
% crossdir = 'D:\OCGOL_learning_long_term\I47_LP\crossSession';

%I45 RT Long term
%  path_dir = {'D:\OCGOL_learning_long_term\I45_RT\I45_RT_AB_d1_062018_1',...
%      'D:\OCGOL_learning_long_term\I45_RT\I45_RT_AB_d6_062518_2',...
%      'D:\OCGOL_learning_long_term\I45_RT\I45_RT_AB_d16_070518_3',...
%      'D:\OCGOL_learning_long_term\I45_RT\I45_RT_AB_d20_070918_4',...
%      'D:\OCGOL_learning_long_term\I45_RT\I45_RT_AB_d25_071318_5',...
%      'D:\OCGOL_learning_long_term\I45_RT\I45_RT_AB_d30_071818_6'};
% %cross session directory
% crossdir = 'D:\OCGOL_learning_long_term\I45_RT\crossSession';

%I46 Long term
%  path_dir = {'D:\OCGOL_learning_long_term\I46\I46_AB_d1_062018_1',...
%      'D:\OCGOL_learning_long_term\I46\I46_AB_d6_062518_2',...
%      'D:\OCGOL_learning_long_term\I46\I46_AB_d16_070518_3',...
%      'D:\OCGOL_learning_long_term\I46\I46_AB_d20_070918_4',...
%      'D:\OCGOL_learning_long_term\I46\I46_AB_d25_071318_5',...
%      'D:\OCGOL_learning_long_term\I46\I46_AB_d30_071818_6'};
% %cross session directory
% crossdir = 'D:\OCGOL_learning_long_term\I46\crossSession';

% %MR2
%  path_dir = {'D:\OCGOL_reversal\MR2\MR2_Random_2022_03_02-001_1',...
%      'D:\OCGOL_reversal\MR2\MR2_Random_2022_03_03-001_2',...
%      'D:\OCGOL_reversal\MR2\MR2_Random_2022_03_04-001_3',...
%      'D:\OCGOL_reversal\MR2\MR2_RevAB_2022_03_05-001_4',...
%      'D:\OCGOL_reversal\MR2\MR2_RevRandom_2022_03_06-001_5',...
%      'D:\OCGOL_reversal\MR2\MR2_RevRandom_2022_03_07-001_6',...
%      'D:\OCGOL_reversal\MR2\MR2_RevRandom_2022_03_12-001_7',...
%      'D:\OCGOL_reversal\MR2\MR2_RevRandom_2022_03_13-001_8',...
%      'D:\OCGOL_reversal\MR2\MR2_RevRandom_2022_03_14-001_9'};
% %cross session directory
% crossdir = 'D:\OCGOL_reversal\MR2\crossSession_update';

% %MR1
%  path_dir = {'D:\OCGOL_reversal\MR1\MR1_Random_2022_02_28-001_1',...
%      'D:\OCGOL_reversal\MR1\MR1_Random_2022_03_01-001_2',...
%      'D:\OCGOL_reversal\MR1\MR1_Random_2022_03_02-002_3',...
%      'D:\OCGOL_reversal\MR1\MR1_RevAB_2022_03_03-001_4',...
%      'D:\OCGOL_reversal\MR1\MR1_RevRandom_2022_03_04-001_5',...
%      'D:\OCGOL_reversal\MR1\MR1_RevRandom_2022_03_05-001_6',...
%      'D:\OCGOL_reversal\MR1\MR1_RevRandom_2022_03_08-001_7',...
%      'D:\OCGOL_reversal\MR1\MR1_RevRandom_2022_03_09-001_8',...
%      'D:\OCGOL_reversal\MR1\MR1_RevRandom_2022_03_10-001_9'};
% %cross session directory
% crossdir = 'D:\OCGOL_reversal\MR1\crossSession_update';

%MR4
 path_dir = {'D:\OCGOL_reversal\MR4\MR4_Random_2022_03_04-001_1',...
     'D:\OCGOL_reversal\MR4\MR4_Random_2022_03_05-001_2',...
     'D:\OCGOL_reversal\MR4\MR4_Random_2022_03_06-001_3',...
     'D:\OCGOL_reversal\MR4\MR4_RevAB_2022_03_07-002_4',...
     'D:\OCGOL_reversal\MR4\MR4_RevRandom_2022_03_08-001_5',...
     'D:\OCGOL_reversal\MR4\MR4_RevRandom_2022_03_09-001_6',...
     'D:\OCGOL_reversal\MR4\MR4_RevRandom_2022_03_11-001_7',...
     'D:\OCGOL_reversal\MR4\MR4_RevRandom_2022_03_12-001_8',...
     'D:\OCGOL_reversal\MR4\MR4_RevRandom_2022_03_13-001_9'};
%cross session directory
crossdir = 'D:\OCGOL_reversal\MR4\crossSession_update';

%% Load cross session registered ROIs
%load auto component registration struct
load(fullfile(crossdir,'registered.mat'))

%all ROI assignent across all sessions (including nans)
ROI_assignments = registered.multi.assigned;

%% load in dataset variables

%load in struct table
for ss=1:size(path_dir,2) %all sessions
    %get directory names of mat files
    files(ss) = subdir(fullfile(path_dir{ss},'input','*.mat'));
    %load directly into memory (fast than memmap with smaller variables)
    m(ss) = load(files(ss).name,'A_keep','C_keep','Cn','Coor_kp','expDffMedZeroed','dims');
end

%find mat files with calcium data
for ss=1:size(path_dir,2) %all sessions

    %load removed ROIs
    removed_ROI_files(ss) = subdir(fullfile(path_dir{ss},'removedROI','*.mat'));
    %load the logical vector
    removed_ROI(ss) = load(removed_ROI_files(ss).name);
        
    %load template (from each session)
    template(ss) = load(fullfile(path_dir{ss},'template_nr.mat'));
    
    %assign template to m struct
    m(ss).template = template(ss).template;     
end

%number of sessions
nb_ses = size(path_dir,2);

%% Move the variables into a structure (all sessions)

%use motion correction template as the stack average for visualization
%purposes
for ss=1:size(path_dir,2) %all sessions
    selector_var(ss).template = m(ss).template;
    selector_var(ss).A_keep = m(ss).A_keep;
    selector_var(ss).C_keep = m(ss).C_keep;
    selector_var(ss).Cn = m(ss).Cn;
    selector_var(ss).Coor_kp = m(ss).Coor_kp;
    selector_var(ss).expDffMedZeroed = m(ss).expDffMedZeroed;
    selector_var(ss).centers = com(m(ss).A_keep,m(ss).dims(1,1),m(ss).dims(1,2));
end

%% Adjust all these variables for removed components 

for ss=1:size(path_dir,2) %all sessions
    
    selector_var(ss).A_keep_filt = selector_var(ss).A_keep(:,removed_ROI(ss).compSelect);
    selector_var(ss).C_keep_filt = selector_var(ss).C_keep(removed_ROI(ss).compSelect,:);
    selector_var(ss).Coor_kp_filt = selector_var(ss).Coor_kp(removed_ROI(ss).compSelect);
    selector_var(ss).expDffMedZeroed_filt = selector_var(ss).expDffMedZeroed(removed_ROI(ss).compSelect,:);
    selector_var(ss).centers_filt = selector_var(ss).centers(removed_ROI(ss).compSelect,:);
end

%% Prep variables for GUI input and run GUI

%remove 'single' assignment from ROI_assignmnents
nan_log_ROI = isnan(ROI_assignments);
remove_singles = find(sum(nan_log_ROI,2) == (nb_ses -1));

%without singles
ROI_assign_multi = ROI_assignments;
ROI_assign_multi(remove_singles,:) = [];

%generate 2D logical
assign_sel_log = ~isnan(ROI_assign_multi);

%start interactive GUI
[output_logical, last_idx] = multi_ses_match_selector_ex(selector_var,ROI_assign_multi,assign_sel_log);

%copy in case global overwrite from GUI
copy_output_logical = output_logical;

%% Update the match matrix

%set to nan filtered out matches
ROI_assign_multi_filtered = ROI_assign_multi;
ROI_assign_multi_filtered(~output_logical) = nan;

%find nans (single match or all nans)
filtered_nan_sum = sum(isnan(ROI_assign_multi_filtered),2);
filtered_idxs = find(filtered_nan_sum == nb_ses | filtered_nan_sum == (nb_ses -1));

%remove from filtered match ROI idx
ROI_assign_multi_filtered(filtered_idxs,:) = [];


%% Check the updated matches with the GUI
%start interactive GUI
%generate 2D logical
assign_sel_log_filtered = ~isnan(ROI_assign_multi_filtered);

multi_ses_match_selector_ex(selector_var,ROI_assign_multi_filtered,assign_sel_log_filtered);

%% Save and export

cd(crossdir) 

%make output directory if it does not exist
try
    mkdir('filtered_match_ROI');
    cd(fullfile(crossdir, ['\','filtered_match_ROI']))
catch
    disp('Directory already exists');
end

%get date
currentDate = date;
%replace dashes with underscores in date
dashIdx = regexp(currentDate,'\W');
currentDate(dashIdx) = '_';

%save all
save([currentDate,'_filtered_match_ROIs.mat'],'ROI_assign_multi_filtered','-v7.3');

%% REMOVE ALL BELOW

%{

%% start the selector GUI
%removedROI in the workspace will be cleared unless the removedROI variable
%is passed as an argument to the fxn
%removedROI - those set to 1 are selected for removal
%somaROI - those set to 1 are selected as soma
%dendriteROI - those set to 1 are selected as dendrites/non-soma

%split into 150 ROIs at a time to prevent slowdown

% % number of assignments from CNMF script
% nbROI = size(ROI_assignments,1);
% 
% %every how many ROIs to do split
% split_ROI = 20;
% 
% %how many splits
% splitNb = ceil(nbROI/split_ROI);
% 
% %ROI index range
% startIdx = 1:split_ROI:nbROI;
% endIdx = [startIdx(2:splitNb)-1,nbROI];
% 
% %split selector_var struct in cell of structs and iterate through
% %pre-define cell
% selector_var_split = cell(1,splitNb);

%split the selector_var struct
% for ii=1:splitNb
%     %copy the original struct
%     selector_var_split{ii} = selector_var;
%     %replace the variables with the split range
%     selector_var_split{ii}.A_keep = selector_var_split{ii}.A_keep(:,startIdx(ii):endIdx(ii));
%     selector_var_split{ii}.C_keep = selector_var_split{ii}.C_keep(startIdx(ii):endIdx(ii),:);
%     selector_var_split{ii}.Coor_kp = selector_var_split{ii}.Coor_kp(startIdx(ii):endIdx(ii));
%     selector_var_split{ii}.expDffMedZeroed = selector_var_split{ii}.expDffMedZeroed(startIdx(ii):endIdx(ii),:);
%     selector_var_split{ii}.centers = selector_var_split{ii}.centers(startIdx(ii):endIdx(ii),:);
% end

%preallocate logicals
% removedROI = false(nbROI,1);
% somaROI = false(nbROI,1);
% dendriteROI = false(nbROI,1);
% 
% %run through each batch of ROIs
% %for ii=1:splitNb
%     [removedROI_temp,somaROI_temp,dendriteROI_temp] = multi_ses_ROI_selector_GUI_V1(selector_var,ROI_assignments);
%     
%     removedROI(startIdx(ii):endIdx(ii)) = removedROI_temp;
%     somaROI(startIdx(ii):endIdx(ii)) = somaROI_temp;
%     dendriteROI(startIdx(ii):endIdx(ii)) = dendriteROI_temp;
%end

%%  Filter out removed soma

%remove ROIs that were accidentally selected in soma
filteredSoma = setdiff(find(somaROI == 1),find(removedROI ==1));

somaClean = zeros(size(somaROI,1),1);
somaClean(filteredSoma) = 1;

%% Revise component selection

%[removedROI,somaROI,dendriteROI] = ROI_selector_GUI_V1(selector_var,removedROI,somaROI,dendriteROI);

%% Update the components

%save component selections into struct
compKeep.removedROI = removedROI;
compKeep.somaROI = logical(somaClean);
compKeep.dendriteROI = dendriteROI;

%vector of ones to indicate the original 1's
originalROI = logical(ones(size(selector_var.A_keep,2),1));
% 
% %1's are those to keep
% compKeep.logicIdx.keepROI = xor(originalROI,compKeep.removedROI);
% 
% %keep those which are somata
% compKeep.logicIdx.somaROIkeep = and(keepROI,compKeep.somaROI);
% 
% %keep those which are non-somata
% compKeep.logicIdx.dendriteROIkeep = and(keepROI,compKeep.dendriteROI);

%% Preview the updated components

%select which components to view
compSelect = compKeep.somaROI;

%make a copy of the original structure
preview_var = selector_var;

%modify the variable to account for selected components
preview_var.A_keep = selector_var.A_keep(:,compSelect);
preview_var.C_keep = selector_var.C_keep(compSelect,:);
preview_var.Coor_kp = selector_var.Coor_kp(compSelect);
preview_var.expDffMedZeroed = selector_var.expDffMedZeroed(compSelect,:);
preview_var.centers = selector_var.centers(compSelect,:);

%preview the components in the GUI
[~,~,~] = ROI_selector_GUI_V1(preview_var);


%% Make separate mat files with updated components

%save the selected componenents into a mat file
% split_dir = split(dir_name,filesep);
% exp_dir = fullfile(split_dir{1:end-1});
% cd(exp_dir)

exp_dir = dir_name;
cd(exp_dir) 

%make output directory if it does not exist
try
    mkdir('removedROI');
    cd(fullfile(exp_dir, ['\','removedROI']))
catch
    disp('Directory already exists');
end

%get date
currentDate = date;
%replace dashes with underscores in date
dashIdx = regexp(currentDate,'\W');
currentDate(dashIdx) = '_';

%to exclude certain variables
%save([currentDate,'.mat'],'-regexp','^(?!(data)$).','-v7.3');

%save all
save([currentDate,'_selectedROIs.mat'],'compSelect','-v7.3');

%}
