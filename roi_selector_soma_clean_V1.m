
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

%load removed components
%remove file directory
files_rem = subdir(fullfile(dir_name,'removedROI','*.mat'));
%load removed logical vector
load(files_rem(1).name,'compSelect');


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

% %
% nbROI = size(selector_var.A_keep,2);
% 
% %every how many ROIs to do split
% split_ROI = 150;
% 
% %how many splits
% splitNb = ceil(size(selector_var.A_keep,2)/split_ROI);
% 
% %ROI index range
% startIdx = 1:split_ROI:nbROI;
% endIdx = [startIdx(2:splitNb)-1,nbROI];
% 
% %split selector_var struct in cell of structs and iterate through
% %pre-define cell
% selector_var_split = cell(1,splitNb);
% 
% %split the selector_var struct
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
% 
% %preallocate logicals
% removedROI = false(nbROI,1);
% somaROI = false(nbROI,1);
% dendriteROI = false(nbROI,1);
% 
% %run through each batch of ROIs
% for ii=1:splitNb
%     [removedROI_temp,somaROI_temp,dendriteROI_temp] = ROI_selector_GUI_V1(selector_var_split{ii});
%     removedROI(startIdx(ii):endIdx(ii)) = removedROI_temp;
%     somaROI(startIdx(ii):endIdx(ii)) = somaROI_temp;
%     dendriteROI(startIdx(ii):endIdx(ii)) = dendriteROI_temp;
% end
% 
% 
% %remove ROIs that were accidentally selected in soma
% filteredSoma = setdiff(find(somaROI == 1),find(removedROI ==1));
% 
% somaClean = zeros(size(somaROI,1),1);
% somaClean(filteredSoma) = 1;

%% Update the components

%save component selections into struct
% compKeep.removedROI = removedROI;
% compKeep.somaROI = logical(somaClean);
% compKeep.dendriteROI = dendriteROI;
% 
% %vector of ones to indicate the original 1's
% originalROI = logical(ones(size(selector_var.A_keep,2),1));

%% Preview and refine the updated components

%select which components to view
%compSelect = compKeep.somaROI;

%make a copy of the original structure
preview_var = selector_var;

%modify the variable to account for selected components
preview_var.A_keep = selector_var.A_keep(:,compSelect);
preview_var.C_keep = selector_var.C_keep(compSelect,:);
preview_var.Coor_kp = selector_var.Coor_kp(compSelect);
preview_var.expDffMedZeroed = selector_var.expDffMedZeroed(compSelect,:);
preview_var.centers = selector_var.centers(compSelect,:);

%%  split into 150 sessions and run

nbROI_clean = size(preview_var.A_keep,2);

%every how many ROIs to do split
split_ROI = 150;

%how many splits
splitNb_clean = ceil(size(preview_var.A_keep,2)/split_ROI);

%ROI index range
startIdx_clean = 1:split_ROI:nbROI_clean;
endIdx_clean = [startIdx_clean(2:splitNb_clean)-1,nbROI_clean];

%split selector_var struct in cell of structs and iterate through
%pre-define cell
selector_var_split_clean = cell(1,splitNb_clean);

%split the selector_var struct
for ii=1:splitNb_clean
    %copy the original struct
    selector_var_split_clean{ii} = preview_var;
    %replace the variables with the split range
    selector_var_split_clean{ii}.A_keep = selector_var_split_clean{ii}.A_keep(:,startIdx_clean(ii):endIdx_clean(ii));
    selector_var_split_clean{ii}.C_keep = selector_var_split_clean{ii}.C_keep(startIdx_clean(ii):endIdx_clean(ii),:);
    selector_var_split_clean{ii}.Coor_kp = selector_var_split_clean{ii}.Coor_kp(startIdx_clean(ii):endIdx_clean(ii));
    selector_var_split_clean{ii}.expDffMedZeroed = selector_var_split_clean{ii}.expDffMedZeroed(startIdx_clean(ii):endIdx_clean(ii),:);
    selector_var_split_clean{ii}.centers = selector_var_split_clean{ii}.centers(startIdx_clean(ii):endIdx_clean(ii),:);
end

%preallocate logicals
removedROI_clean = false(nbROI_clean,1);
somaROI_clean = false(nbROI_clean,1);
dendriteROI_clean = false(nbROI_clean,1);

%run through each batch of ROIs
for ii=1:splitNb_clean
    [removedROI_temp_clean,somaROI_temp_clean,dendriteROI_temp_clean] = ROI_selector_GUI_V1(selector_var_split_clean{ii});
    removedROI_clean(startIdx_clean(ii):endIdx_clean(ii)) = removedROI_temp_clean;
    somaROI_clean(startIdx_clean(ii):endIdx_clean(ii)) = somaROI_temp_clean;
    dendriteROI_clean(startIdx_clean(ii):endIdx_clean(ii)) = dendriteROI_temp_clean;
end


%% Check the soma cleaned componenets
%make a copy of the original structure
final_var = selector_var;

%modify the variable to account for selected components
final_var.A_keep = preview_var.A_keep(:,~removedROI_clean);
final_var.C_keep = preview_var.C_keep(~removedROI_clean,:);
final_var.Coor_kp = preview_var.Coor_kp(~removedROI_clean);
final_var.expDffMedZeroed = preview_var.expDffMedZeroed(~removedROI_clean,:);
final_var.centers = preview_var.centers(~removedROI_clean,:);

%preview the components in the GUI
[~,~,~] = ROI_selector_GUI_V1(final_var);


%% Make separate mat files with soma cleaned componenets for D1

%save the selected componenents into a mat file
% split_dir = split(dir_name,filesep);
% exp_dir = fullfile(split_dir{1:end-1});
% cd(exp_dir)

exp_dir = dir_name;
cd(exp_dir) 

%make output directory if it does not exist
try
    mkdir('removedROI_clean');
    cd(fullfile(exp_dir, ['\','removedROI_clean']))
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

save([currentDate,'_removedROI_clean.mat'],'removedROI_clean','-v7.3');

