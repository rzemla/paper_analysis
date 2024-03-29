%% Load in data from each animal

%session directories
%home
session_dirs = {'F:\I46_AB_d1_062018'};
%workstation
%session_dirs = {'G:\lec_paper_data\OCGOL_nonsilenced\I46_AB_d1_062018'};

%mat filenames containing data
%make sure only 1 mat file is present in output directory
for ii=1:size(session_dirs,1)
    file_names{ii} = dir(fullfile(session_dirs{ii},'output','*.mat'));
end

%load in processed data for each session (all variables)
for ii=1:size(session_dirs,1)
    animal_data{ii} = load(fullfile(file_names{ii}.folder, file_names{ii}.name));
end

%% Defined tuned logical vectors

[tunedLogical] = defineTunedLogicals(animal_data);

%% Generate STCs neurons tuned in either session and plot side by side
%customize to add options
%tuned in both sessions by SI score
%sorted by A trials

plotSTC_OCGOL(animal_data,tunedLogical)


%% Spiral plots of individual ROIs

%where to make avi video of selected plotd
options.makeVideo = 0;
options.videoName = 'cellEvents_all';

%speed at which each event transitions to the next
options.plotSpeed = 0.001;
%plot all on top of one another vs each individually
options.hold = 0;

%add silence figure option to code
options.suppressFigure = 1;

%to manually move through the events
options.manualAdvance = 0;

% choose cells with sig tuning fields
%sigPF = Place_cell{1}.placeField.sig_ROI  & Place_cell{2}.placeField.sig_ROI;
%select single place fields
%sigPF_1a = find(Place_cell{1}.placeField.nb ==1);
%sigPF_1b = find(Place_cell{2}.placeField.nb ==1);

%single place field (regardless of signifncance)
%ROIrange = intersect(sigPF_1a,sigPF_1b);
ROIrange = find(AandB_tuned == 1);
ROIrange = find(onlyA_tuned == 1);

spiralEvents = event_spiral(animal_data{1}, ROIrange,options);

%% Pie chart of fraction of neurons tuned in each subset (A,B, both, neither)


%% Scatter plot of calcium properties of each neuron between A and B trials
%AUC/min, event rate, 

%temporary AUC/min calculation during run epochs
%move to this event properties script

%A
run_AUC{1} = animal_data{1}.Events_split{1}.Run.properties.AUC;
%B
run_AUC{2} = animal_data{1}.Events_split{2}.Run.properties.AUC;

%A and B run time (minutes)
runtime(1) = (size(animal_data{1}.Imaging_split{1}.time_restricted,1)*animal_data{1}.Imaging.dt)./60;
runtime(2) = (size(animal_data{1}.Imaging_split{2}.time_restricted,1)*animal_data{1}.Imaging.dt)./60;

%for each trial type
for tt=1:2
    %for each  ROI
    for rr=1:size(run_AUC{tt},2)
        AUC_min{tt}(rr) = sum(run_AUC{tt}{rr})./runtime(tt);
    end
end

%% Plot the scatterplots - AUC/min between trial types

event_scatterplots(AUC_min,tunedLogical);


%% Tuning specificity, fraction on cells in tuned in each trial by SI and TS scores

%% Remapping tuning specificity distance (like done in GOL)









