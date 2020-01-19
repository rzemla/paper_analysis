function supplement_multi_day_recall_data_with_all_A_allB_Place_data(path_dir)
%Calculates the spatial properties (maps, occupancy etc.) on all A and B
%laps (correct and incorrect) and 
%Input: cell containing the main directory of the animal
%run per animal

%% Paths associated for each animal

%I42L 1
%path_dir = {'G:\OCGOL_stability_recall\I42L_1\I42L_AB_d1_032118_1'};


%% Load relevant variables for generating 4 and 5 trials for recall animals

%find the .mat file associated with the analysis in the output folder
mat_file_nam = dir([path_dir{1},'\output\*.mat']);

%check that there is only 1 mat file in the directory, if more than 1
%y throw error
if size(mat_file_nam,1) ~= 1 %works
    msg = 'More than 2 mat files detected. Stopping script!';
    %throw error
    error(msg)
end

%load in relevant variables
load(fullfile(mat_file_nam.folder, mat_file_nam.name),'Place_cell','Behavior_split','Events_split','Imaging_split','options');


%% !!!! Make back-up copy of the input Place_cell struct other 1:3 trials data will be lost !!!

%technically should be fine actually to pass in and only run for 4 and 5
%trials 
%make copy
Place_cell_original = Place_cell;

%delete original place cell variable
clear Place_cell

%% Identification of spatially-tuned cells - works RZ with new events - split this
% Set parameters
options.sigma_filter=3; % Sigma (in bin) of gaussian filter (Danielson et al.2016 = 3 bins)
options.smooth_span=3; % span for moving average filter on dF/F (Dombeck 2010 = 3) (Sheffield 2015 = 3)
options.minevents=3; % Min nb of events during session
options.Nbin=[2;4;5;8;10;20;25;100]; % Number of bins to test ([2;4;5;8;10;20;25;100] Danielson et al. 2016)
options.bin_spatial_tuning=100; % Number of bins to compute spatial tuning curve (rate map) -value must be in options.Nbin
options.Nshuffle=1000; % Nb of shuffle to perform
options.pvalue=0.05; % Min p value to be considered as significant
options.dispfig=1; % Display figure 

%use binned position rather than raw for tuning specificity calculation
options.binPosition = 1;

tic;
%for all type of trials --> current: 1 -A trials; 2 -B trials; 3 - all
%laps/trials
for ii=[4, 5]
    %work on this part
    %spatial binning, rate maps, and spatial tuning score
    disp('Calculate bin space, events, rate maps, generate STCs, and SI score')
    [Place_cell{ii}] = spatial_properties(Behavior_split{ii}, Events_split{ii}, Imaging_split{ii},options);
    %tuning specific calculation function here
    disp('Calculate turining specificity score')
    [Place_cell{ii}] = tuning_specificity_RZ_V2(Place_cell{ii},Behavior_split{ii},Events_split{ii},options);
    
end
toc;

%% Merge with original Place_cell struct and overwrite only this variable

%load in original 1:3 data
Place_cell_merge(1:3) = Place_cell_original(1:3);

%load in data generate for 4:5 data
Place_cell_merge(4:5) = Place_cell(4:5);

%% Overwrite the Place_cell struct with the merge one now

Place_cell = Place_cell_merge;

%% Overwrite the place cell struct in the output folder


save(fullfile(mat_file_nam.folder, mat_file_nam.name),'Place_cell','-append');

disp('Successfully replace Place_cell struct variable!')
end
