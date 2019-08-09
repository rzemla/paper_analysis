
%% load in dataset variables
clear

%select directory (location of experiment directory)
dir_name = uigetdir;

%find mat files with calcium data
files = subdir(fullfile(dir_name,'input','*.mat'));

%load directly into memory (fast than memmap with smaller variables)
m = load(files(1).name,'A_keep','C_keep','Cn','Coor_kp','expDffMedZeroed','dims');

%load template
load(fullfile(dir_name,'template_nr.mat'))

%check if template was saved, if not use stack average as template
m.template = template;

%% Move the variables into a structure

%use motion correction template as the stack average for visualization
%purposes
selector_var.template = m.template;
selector_var.A_keep = m.A_keep;
selector_var.C_keep = m.C_keep;
selector_var.Cn = m.Cn;
selector_var.Coor_kp = m.Coor_kp;
selector_var.expDffMedZeroed = m.expDffMedZeroed;
selector_var.centers = com(m.A_keep,m.dims(1,1),m.dims(1,2));


%% start the selector GUI
%removedROI in the workspace will be cleared unless the removedROI variable
%is passed as an argument to the fxn
%removedROI - those set to 1 are selected for removal
%somaROI - those set to 1 are selected as soma
%dendriteROI - those set to 1 are selected as dendrites/non-soma

%split into 150 ROIs at a time to prevent slowdown

%
nbROI = size(selector_var.A_keep,2);

%every how many ROIs to do split
split_ROI = 150;

%how many splits
splitNb = ceil(size(selector_var.A_keep,2)/split_ROI);

%ROI index range
startIdx = 1:split_ROI:nbROI;
endIdx = [startIdx(2:splitNb)-1,nbROI];

%split selector_var struct in cell of structs and iterate through
%pre-define cell
selector_var_split = cell(1,splitNb);

%split the selector_var struct
for ii=1:splitNb
    %copy the original struct
    selector_var_split{ii} = selector_var;
    %replace the variables with the split range
    selector_var_split{ii}.A_keep = selector_var_split{ii}.A_keep(:,startIdx(ii):endIdx(ii));
    selector_var_split{ii}.C_keep = selector_var_split{ii}.C_keep(startIdx(ii):endIdx(ii),:);
    selector_var_split{ii}.Coor_kp = selector_var_split{ii}.Coor_kp(startIdx(ii):endIdx(ii));
    selector_var_split{ii}.expDffMedZeroed = selector_var_split{ii}.expDffMedZeroed(startIdx(ii):endIdx(ii),:);
    selector_var_split{ii}.centers = selector_var_split{ii}.centers(startIdx(ii):endIdx(ii),:);
end

%preallocate logicals
removedROI = false(nbROI,1);
somaROI = false(nbROI,1);
dendriteROI = false(nbROI,1);

%run through each batch of ROIs
for ii=1:splitNb
    [removedROI_temp,somaROI_temp,dendriteROI_temp] = ROI_selector_GUI_V1(selector_var_split{ii});
    removedROI(startIdx(ii):endIdx(ii)) = removedROI_temp;
    somaROI(startIdx(ii):endIdx(ii)) = somaROI_temp;
    dendriteROI(startIdx(ii):endIdx(ii)) = dendriteROI_temp;
end

% removedROI = or(removedROI,removedROI_2);
% somaROI = or(somaROI,somaROI_2);

%view selected ROI

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

