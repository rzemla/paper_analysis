function [outputArg1,outputArg2] = cumulative_task_sel_STC(path_dir)


%% Load the stc data

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





%% Plot color-coded STCs

end

