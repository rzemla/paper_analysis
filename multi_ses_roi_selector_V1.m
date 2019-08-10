%%
clear
%% List of datasets

%input directories to matching function
 path_dir = {'G:\OCGOL_learning_short_term\I56_RTLS\I56_RLTS_5AB_041019_1',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_5AB_041119_2',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_3A3B_041219_3',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_3A3B_041319_4',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_no_punish_041519_5',...
     'G:\OCGOL_learning_short_term\I56_RTLS\I56_RTLS_ABrand_no_punish_041619_6'};
%cross session directory
crossdir = 'G:\OCGOL_learning_short_term\I56_RTLS\crossSession';

%% Load cross session registered ROIs
%load auto component registration struct
load(fullfile(crossdir,'registered.mat'))

%all ROI assignent across all sessions (including nans)
ROI_assignments = registered.multi.assigned;

%% load in dataset variables

%find mat files with calcium data
for ss=1:size(path_dir,2) %all sessions
    %get directory names of mat files
    files(ss) = subdir(fullfile(path_dir{ss},'input','*.mat'));
    
    %load directly into memory (fast than memmap with smaller variables)
    %m(ss) = load(files(ss).name,'A_keep','C_keep','Cn','Coor_kp','expDffMedZeroed','dims');
    
    %load removed ROIs
    removed_ROI_files(ss) = subdir(fullfile(path_dir{ss},'removedROI','*.mat'));
    %load the logical vector
    removed_ROI(ss) = load(removed_ROI_files(ss).name);
        
    %load template (from each session)
    template(ss) = load(fullfile(path_dir{ss},'template_nr.mat'));
    
    %assign template to m struct
    m(ss).template = template(ss).template;
end

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


%% start the selector GUI
%removedROI in the workspace will be cleared unless the removedROI variable
%is passed as an argument to the fxn
%removedROI - those set to 1 are selected for removal
%somaROI - those set to 1 are selected as soma
%dendriteROI - those set to 1 are selected as dendrites/non-soma

%split into 150 ROIs at a time to prevent slowdown

% number of assignments from CNMF script
nbROI = size(ROI_assignments,1);

%every how many ROIs to do split
split_ROI = 20;

%how many splits
splitNb = ceil(nbROI/split_ROI);

%ROI index range
startIdx = 1:split_ROI:nbROI;
endIdx = [startIdx(2:splitNb)-1,nbROI];

%split selector_var struct in cell of structs and iterate through
%pre-define cell
selector_var_split = cell(1,splitNb);

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
removedROI = false(nbROI,1);
somaROI = false(nbROI,1);
dendriteROI = false(nbROI,1);

%run through each batch of ROIs
%for ii=1:splitNb
    [removedROI_temp,somaROI_temp,dendriteROI_temp] = multi_ses_ROI_selector_GUI_V1(selector_var,ROI_assignments);
    
    removedROI(startIdx(ii):endIdx(ii)) = removedROI_temp;
    somaROI(startIdx(ii):endIdx(ii)) = somaROI_temp;
    dendriteROI(startIdx(ii):endIdx(ii)) = dendriteROI_temp;
%end
%% Working version


%remove 'single' assignment from ROI_assignmnents
nan_log_ROI = isnan(ROI_assignments);
remove_singles = find(sum(nan_log_ROI,2) == 5);

%without singles
ROI_assign_multi = ROI_assignments;
ROI_assign_multi(remove_singles,:) = [];

%generate 2D logical
assign_sel_log = ~isnan(ROI_assign_multi);

multi_ses_match_selector(selector_var,ROI_assign_multi,assign_sel_log)

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

