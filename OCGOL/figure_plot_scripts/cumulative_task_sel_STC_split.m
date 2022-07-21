function [task_sel_STC_data] = cumulative_task_sel_STC_split(path_dir)


%% Load the STC data

%read relevant data
for aa=1:size(path_dir,2)
    load_data_path{aa} = fullfile(path_dir{aa},'cumul_analysis','task_sel_STC.mat');
    cumulative_STCs{aa} = load(string(load_data_path{aa}));
end

%% Generate a cumulative raster

%combine maps into single cell (animal x rasters)
%preallocate cell
STC_cell = cell(size(path_dir,2), 2);

for aa=1:size(path_dir,2)
    STC_cell(aa,:) = cumulative_STCs{aa}.task_sel_STC.maps{1};
end

%Split into A/B STC matrices
STC_A_sel = cell2mat(STC_cell(:,1));
STC_B_sel = cell2mat(STC_cell(:,2));

%% Sort each animal STC by max spatial tuning by A-sel and B-sel place cells

%for A selective place cells
for aa=1:size(path_dir,2)
    %output bin location where maximum norm value
    [~,maxBin_each_temp.A] = max(STC_cell{aa,1}(:,1:100)', [], 1);
    %sort by max spatial bin
    [~,sortOrder_each_temp{aa}.A] = sort(maxBin_each_temp.A,'ascend');
    %sort each A task-sel raster
    STC_task_sel_sorted_each{aa,1} =STC_cell{aa,1}(sortOrder_each_temp{aa}.A,:);
end

%for B selective place cells
for aa=1:size(path_dir,2)
    %output bin location where maximum norm value
    [~,maxBin_each_temp.B] = max(STC_cell{aa,2}(:,101:200)', [], 1);
    %sort by max spatial bin
    [~,sortOrder_each_temp{aa}.B] = sort(maxBin_each_temp.B,'ascend');
    %sort each A task-sel raster
    STC_task_sel_sorted_each{aa,2} =STC_cell{aa,2}(sortOrder_each_temp{aa}.B,:);
end

%% Global sort of A/B selective STCs

[~,maxBin_all.A] = max(STC_A_sel(:,1:100)', [], 1);
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,sortOrder_all.A] = sort(maxBin_all.A,'ascend');

[~,maxBin_all.B] = max(STC_B_sel(:,101:200)', [], 1);
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,sortOrder_all.B] = sort(maxBin_all.B,'ascend');

% figure
% imagesc(STC_A_sel(sortOrder_all.A,:))
% hold on
% 
% figure
% imagesc(STC_B_sel(sortOrder_all.B,:))
% hold on
%% Export data for external plotter

%reward A and B black marker lines
task_sel_STC_data.rewA_bin_pos = 30;
task_sel_STC_data.rewB_bin_pos = 70;
task_sel_STC_data.odor_bin_pos = 10;

%A sel on A laps
task_sel_STC_data.AselA = STC_A_sel(sortOrder_all.A,1:100);
%A sel on B laps
task_sel_STC_data.AselB = STC_A_sel(sortOrder_all.A,101:200);
%B sel on A laps
task_sel_STC_data.BselA = STC_B_sel(sortOrder_all.B,1:100);
%B sel on B laps
task_sel_STC_data.BselB = STC_B_sel(sortOrder_all.B,101:200);

%export sorted A and B selective STCs for each animal
task_sel_STC_data.split.STCs = STC_task_sel_sorted_each;
task_sel_STC_data.split.sort_order = sortOrder_each_temp;

end
