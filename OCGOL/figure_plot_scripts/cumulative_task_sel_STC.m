function [task_sel_STC_data] = cumulative_task_sel_STC(path_dir)


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



end

