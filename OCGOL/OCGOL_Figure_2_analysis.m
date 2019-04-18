%% Load in data from each animal

%session directories
%home
%session_dirs = {'F:\I46_AB_d1_062018'};
%workstation
session_dirs = {'G:\lec_paper_data\OCGOL_nonsilenced\I46_AB_d1_062018'};

%mat filenames containing data
%make sure only 1 mat file is present in output directory
for ii=1:size(session_dirs,1)
    file_names{ii} = dir(fullfile(session_dirs{ii},'output','*.mat'));
end

%load in processed data for each session (all variables)
for ii=1:size(session_dirs,1)
    animal_data{ii} = load(fullfile(file_names{ii}.folder, file_names{ii}.name));
end


%% Generate STCs neurons tuned in either session and plot side by side
%tuned in both sessions by SI score
%sorted by A trials




%% Pie chart of fraction of neurons tuned in each subset (A,B, both, neither)

%% 