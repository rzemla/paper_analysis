%% Set parameters for the pipeline

options.defineDir = 1;

%setDir = 'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1';

%setDir = 'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1';

%setDir = 'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1';

%setDir = 'G:\Figure_2_3_selective_remap\I52RT_AB_sal_120618_1';
setDir = 'G:\Figure_2_3_selective_remap\I53LT_AB_sal_113018_1';


%whether to define experiment directory or use GUI to select
%1 = define in variable, 0 = GUI select

%performance/texture code works - backup this dataset!
%setDir = 'E:\I42L_AB_1\I42L_AB_d1_032118_1';

%performance/texture code works - backup this dataset!
%setDir = 'E:\I42R_AB_1\I42R_AB_d1_032118_1';


%setDir = 'G:\OCGOL_training\I56_RTLS_ABrand_no_punish_041519';
%setDir = 'G:\OCGOL_training\I57_RLTS_5AB_041019';

%setDir = 'G:\lec_paper_data\OCGOL_nonsilenced\I46_AB_d1_062018';
%setDir = 'F:\OCGOL_training\I56_RLTS_041019\5A5B';
%setDir = 'F:\OCGOL_training\I56_RLTS_041019\ABrand_no_punish_041619';



%performance/texture code works - backup this dataset!
%check missed reward collection feature
%setDir = 'G:\lec_paper_data\OCGOL_nonsilenced\I47_LP_AB_d1_062018';

%performance/texture code works - backup this dataset!
%setDir = 'G:\lec_paper_data\OCGOL_nonsilenced\I46_AB_d1_062018';

%performance/texture code works - backup this dataset!
%setDir = 'G:\lec_paper_data\OCGOL_nonsilenced\I45_RT_AB_d1_062018';

%performance/texture code works - backup this dataset!
%setDir ='G:\I52RT_AB_sal_PSEM_113018_120118\1';

%performance/texture code works
%setDir ='G:\I52RT_AB_sal_PSEM_113018_120118\2';
%performance/texture code works
%setDir ='G:\I52RT_AB_sal_PSEM_113018_120118\3';

%performance/texture code works
%setDir ='G:\I52RT_AB_sal_PSEM_113018_120118\4';

%performance/texture code works
%setDir ='G:\I53LT_AB_sal_PSEM_113018_120118\1';
%performance/texture code works
%setDir ='G:\I53LT_AB_sal_PSEM_113018_120118\2';
%performance/texture code works
%setDir ='G:\I53LT_AB_sal_PSEM_113018_120118\3';
%performance/texture code works
%setDir ='G:\I53LT_AB_sal_PSEM_113018_120118\4';

%setDir = 'E:\I42L_AB_1\I42L_AB_d1_032118_1';

%setDir = 'G:\GOL\I55_RTLS_RF_GOL_022219\1'; %(RF of GOL)

%setDir = 'G:\GOL\I55_RTLS_RF_GOL_022219\2'; % (Block 1 day 1)

%setDir = 'G:\GOL\I55_RTLS_RF_GOL_022219\6'; %(Block 2 day 5)

%home workstation directory
%setDir = 'F:\I55_RTLS_RF_GOL_022219\2'; % (Block 1 day 1)

%setDir = 'F:\I52RT_AB_sal_PSEM_113018_120118\1';

%setDir = 'G:\GOL\I55_RTLS_RF_GOL_022219\3';
%setDir = 'G:\GOL\I55_RTLS_RF_GOL_022219\4';
%setDir = 'G:\GOL\I55_RTLS_RF_GOL_022219\5';

 %setDir = 'G:\GOL\I55_RTLS_RF_GOL_022219\7';
 
 %setDir = 'G:\GOL\LEC_silenced\I56_RTLS_RF_GOL_032019\3';
 %setDir = 'G:\GOL\LEC_silenced\I56_RTLS_RF_GOL_032019\4';
 %setDir = 'G:\GOL\LEC_silenced\I57_RTLS_RF_GOL_032019\3';
  
%setDir = 'G:\GOL\I54_LTRS_RF_GOL_022219\1';
%setDir = 'G:\GOL\I54_LTRS_RF_GOL_022219\2';
%setDir = 'G:\GOL\I54_LTRS_RF_GOL_022219\2';
%setDir = 'G:\GOL\I54_LTRS_RF_GOL_022219\3';
%setDir = 'G:\GOL\I54_LTRS_RF_GOL_022219\7';

%setDir = 'G:\lec_paper_data\OCGOL_nonsilenced\I45_RT_AB_d1_062018';
%setDir = 'G:\lec_paper_data\OCGOL_nonsilenced\I46_AB_d1_062018';

%see an effect
%setDir = 'G:\I57_RTLS_RF_bas_sal_792_031419\3';
%setDir = 'G:\I57_RTLS_RF_bas_sal_792_031419\1';

%setDir = 'G:\I57_LT_RF_bas_sal_792_031419\1';


%setDir = 'G:\I56_RTLS_RF_bas_sal_792_031419\1';
%setDir = 'G:\I56_RTLS_RF_bas_sal_792_031419\3';
%setDir = 'G:\I56_RTLS_RF_bas_sal_792_031419\2';
%setDir ='G:\I56_RT_RF_bas_sal_792_031519\3';

%whether to load in existing XML and CSV behavioral data save in workspace
%1 - load from saved workspace
%0 - read and load from raw XML and CSV files
options.loadBehaviorData = 0;

%whether to load in previously read imaging data
options.loadImagingData = 0;

%choose the behavior that the animal ran
% RF, GOL-RF (GOL day 0), GOL, OCGOL
options.BehaviorType = 'OCGOL';

%type of calcium data
options.calcium_data_input = 'CNMF';

%% TODO

%insert selector here for type of behavior; RF, GOL, OCGOL

%% Define input directory for experiment (.mat, .csv, .xml data)

if options.defineDir == 0
    %select with GUI
    directory_name =  uigetdir;
elseif options.defineDir == 1
    %assign from setDir
    directory_name = setDir;
end
%save dir_name into options struct
options.dir_name = directory_name;

%% Import behavioral data and store as MATLAB workspace or load from workspace

tic;
if options.loadBehaviorData == 0
    [CSV,XML] = readBehaviorData_V1(directory_name);
elseif options.loadBehaviorData == 1
    [CSV,XML] = loadBehaviorData_V1(directory_name);
end
toc;

%% Import imaging data and compute initial dF/F - all ROIs
%output struct of dF/F-related variables
%takes around ~3-4 min

tic;
if options.loadImagingData == 0
    %read CNMF data
    [CNMF_output] = readImagingData_V1(directory_name,options);
    %compute dF/F 
    [F_vars] = compute_dFF_V1(CNMF_output,directory_name);
elseif options.loadImagingData == 1
    %read in the same data from file
    [F_vars,CNMF_output] = loadImagingData_V1(directory_name);
end
toc;

%select the calcium imaging data to be used in subsequent analysis
C_df = F_vars.F_dff_exp';

%% Include only pre-selected somatic ROI:

%whether to remove ROI
options.removeROI = 1;

if options.removeROI == 1
    %remove the respective ROIs from F_vars struct and C_df
    [C_df,F_vars] = removeROI_V2(C_df,F_vars,directory_name);
    disp('Non-somatic ROIs removed');
end

%% Extract position and laps

%minimum distance between sequential laps tag (cm)
options.mindist=10;

%behavior sampling rate (10 kHz) in Hz
options.acqHz=10000;

options.dispfig=1; % Display figure (lap, RFID)

[Behavior] = extractPositionAndLaps(CSV,options);


%% Extract texture/reward locations - there is an out of bounds issue with this

%whether to extract textures (and other signals)
options.textures = true;

if options.textures == true
    %[Behavior] = extractTextures(CSV, Behavior, options);
    %all signals, not just texture related signals
    switch options.BehaviorType
        
        case 'GOL-RF' %special case for day 1 of GOL
            [Behavior] = extractTextures_GOL_RF(CSV, Behavior, options);
        case 'GOL'
            [Behavior] = extractTextures_GOL(CSV, Behavior, options);
        case 'OCGOL'
            [Behavior] = extractTextures_OCGOL(CSV, Behavior, options);
    end
end

%% After here CSV raw data matrix should not be necessary 

%% Extract behavior performance

%add selector here depending on type of behavior run
switch options.BehaviorType
    case 'RF' %same as GOL for now
        [Behavior] = GOL_performance(Behavior, CSV);
    case 'GOL-RF'
        [Behavior] = GOL_RF_performance_new_inputs(Behavior);
    case 'GOL'
        %update lick struct in behavior to retain info about restricted
        %licks and median reward position
        [Behavior] = GOL_performance_new_inputs(Behavior);
    case 'OCGOL'
        %works well for (I52RT AB PSEM/sal 113018)
        %fault with missed reward (I53LT AB PSEM/sal 113018)
        %in licks plot shade the reward zones
        %for reward collected as well
        [Behavior] = OCGOL_performance_new_inputs(Behavior);
end

%% Save Behavior struct temporarily here - later do at end

% fprintf('Saving behavioral data...');
% save(fullfile(directory_name,'output','Behavior.mat'),'Behavior');
% fprintf('Done\n');

%% Restrict data to complete laps

% Restrict calcium data to selected lap -  restrict trace to only full laps
options.restrict=1; 

%Select from which starting lap
options.startlap= 'first'; 

%Select the end lap
options.endlap='last';

%display the calcium trace aligned to behavior figure
options.dispfig=1; % Display figure

%neurons to display
options.ROI =40; 

[Imaging, Behavior, F_vars] = restrict_data(C_df, Behavior, XML, options,F_vars);

%% Event detection -initial detection with 2 sigma onset - WORKS (for restricted as well)!

options.update_events = 0;
options.imaging_rate = 30; %in Hz
options.mindurevent = 1; %in s (250 ms Zaremba 2017 / 1 s Danielson 2016)
options.restrict = true; %only analyze data from complete laps 

%which ROI to plot
options.ROI = 40;
%how many iterations (after initial detection)
options.it = 2;
%set to zero before entering iteration stage in the later cell (initialization)
options.iterationNb = 0;

%go through this and the recalculation function 
%make sure that the recalculation function uses the same inputs as the
%initial function used to compute the dF/F
tic;
[Events] = detect_events(options, Imaging);
toc; 

%% Remove events, recalculate baseline and dF/F, baseline sigma and update events iteratively - WORKS (for restricted as well)!

%add a plot summary for this

%whether to use recalculate sigma and mean to detect events
options.update_events = 1;
%temporary
options.dff_type = 'Rolling median';

%check; imported from OCGOL
%add pass for dF/F calculation parameters to function
tic
for ii =1:options.it
    %which iteration the loop is on - used for duration filter
    %only filters on last iteration of the loop
    options.iterationNb = ii;
    %Update dff and sigma
    updated_dff = update_dff_revised(Events.onset_offset, F_vars, options);
    
    % Update events
    [Events] = detect_events(options, Imaging, updated_dff);
    
end
toc

%% Update imaging struct with recalculated dF/F with event masking - TODO

%% Resample behavioral data
%todo - add resampling of updated dF/F in imaging struct

[Behavior] = resample_behavioral_data(Behavior,Imaging,options);

%% Resample data and determine run epochs

%Danielson or Cossart method
options.method='peak'; %'peak' (Danielson) or 'speed'

% threshold running epochs based on :
% min peak speed(Danielson et al. 2016a, b, 2016) 
% OR average speed (Cossart) 
options.moving_window=10; % window width for moving mean filter and speed filter - Cossart too
%options.minspeed=2; %minimum speed (cm/s)  -- only if 'speed' - Cossart
%too
%merging is also done for Cossart in the script
%minimum duration criterion applies as well

%Danielson criteria
options.minpeak=5;  % minimum peak speed (cm/s) -- only if 'peak' %Danielson 2016b that the animal has to reach within epoch
options.mindur=1; %Minimum duration (s) for running epoch Danielson 2016b
options.merge=0.5; %Merge consecutive running epochs separated by less than (s) Danielson 2016b
options.minspeed = 0; %minimum speed for run epoch threshold (can raise in increments when animal does micromovements around stop point)
options.dispfig=1; % Display figure

%which ROI to display
options.c2plot = 40;

%TODO - merge updated dF/F from above 
[Events, Behavior] = determine_run_epochs(Events,Behavior,Imaging, updated_dff, options);

%% Split into behavior, imaging, event data into individual laps - work on this!

%modify to split laps irrespective of behavior - save into a single struct
%[Behavior_split,Imaging_split,Events_split,Behavior_split_lap,Events_split_lap] = split_laps(Behavior,Imaging,CSV,Events);

%% Split trials based on OCGOL performance

%TODO - split also the updated dF/F traces

switch options.BehaviorType
    case 'RF' %work on this
        %[Behavior_split,Imaging_split,Events_split,Behavior_split_lap,Events_split_lap] = split_laps(Behavior,Imaging,CSV,Events);
    case 'GOL-RF' %work on this
        
    case 'GOL' %work on this
        
    case 'OCGOL'
        [Behavior_split,Imaging_split,Events_split,Behavior_split_lap,Events_split_lap] = split_trials_OCGOL(Behavior,Imaging,Events,options);
end

%% Extract calcium event properties for split and all laps - check this

%for OCGOL data
%for each of 3 conditions = 1 - all, 2 - run, 3 = norun (inside function)
%for each trial type (trial type)


%TODO - update this using the updated dF/F imaging struct
for ii=1:5
    Events_split{ii} = event_properties(updated_dff, Events_split{ii},options);
end

%without split on all lap restricted data
Events = event_properties(updated_dff, Events,options);

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
for ii=1:3
    %work on this part
    %spatial binning, rate maps, and spatial tuning score
    disp('Calculate bin space, events, rate maps, generate STCs, and SI score')
    [Place_cell{ii}] = spatial_properties(Behavior_split{ii}, Events_split{ii}, Imaging_split{ii},options);
    %tuning specific calculation function here
    disp('Calculate turining specificity score')
    [Place_cell{ii}] = tuning_specificity_RZ_V2(Place_cell{ii},Behavior_split{ii},Events_split{ii},options);
    
end
toc; 

%% Save relevant variables for shuffle processing

%save(fullfile(directory_name,))

%% Run place cell shuffle for spatial information
%offload this to HPC for processing
%check what the difference is between 
for ii=1:3
    [Place_cell{ii}] = shuffle_place_cell(Place_cell{ii},Behavior_split{ii},Events_split{ii},options);
end
%alternative way of calculating shuffle for SI and TS - not finalized
%[Place_cell{ii}]=shuffle_place_cell_spatial_V2_RZ(Place_cell{ii},Behavior_split{ii},Events_split{ii},options);

%% Create tuned ROI binary mask for tuned cells and add to Place_cell struct
%add this to shuffle script and remove from there
%nb of ROIs
ROInb = size(Imaging.trace,2);

%for A, B, and all laps
for ii=1:size(Place_cell,2)
    tunedROImask = zeros(1,ROInb);
    tunedROImask(Place_cell{1,ii}.Tuned_ROI) = 1;
    Place_cell{1,ii}.Tuned_ROI_mask =  tunedROImask;
end


%% save the workspace

%make output folder containing the saved variables in the experiments
%directory
save_path = directory_name;

try
    cd(save_path)
    mkdir('output');
    cd(fullfile(save_path, ['\','output']))
catch
    disp('Directory already exists');
    cd(fullfile(save_path, ['\','output']))
end

%get date
currentDate = date;
%replace dashes with underscores in date
dashIdx = regexp(currentDate,'\W');
currentDate(dashIdx) = '_';

%to exclude certain variables
%save([currentDate,'.mat'],'-regexp','^(?!(data)$).','-v7.3');

%save relevant variable for further analysis
tic
save([currentDate,'_ca_analysis.mat'],...
    'directory_name','Place_cell','Events_split','Behavior_split','Imaging_split',...
    'Behavior','Imaging', 'Events', 'updated_dff',...
    'F_vars', 'options', 'Behavior_split_lap',...
    'Events_split_lap','-v7.3');
toc

disp('Done saving.');

%% For modifcation below
%{

%% Place field finder (Dombeck/Barthos 2018)
%Place field finder - Dombeck/Harvey 2010
tic;
%imaging period
options.dt = Imaging.dt;
%number of shuffles
options.nbShuffle = 1000;

%detect place fields and return - real
for ii=1
    [Place_cell{ii}] =  placeFieldDetection_barthosBeta(Place_cell{ii});

    
    if 0
    %shuffled the dF/F traces
    %1 min for 1000 shuffles of indices
    trace_shuffled_mat = shuffle_dFF(Imaging_split{ii},options);
    
    %recalculate the inputs parameters for each shuffle run
    %preallocate shuffle matrix with # of place fields per ROI
    placeField_nb_shuffle = zeros(options.nbShuffle,171);
    PF_inputs = cell(1,options.nbShuffle);
    
    tic;
    parfor n=1:options.nbShuffle
        %generate inputs based on shuffle
        PF_inputs{n} = activity_map_generator_PF(Behavior_split{ii}, Events_split{ii}, trace_shuffled_mat(:,:,n),options);
        %detect place fields and return - shuffled
        [placeField_nb_shuffle(n,:)] =  placeFieldDetection_barthosBeta_shuffle(PF_inputs{n});
    end
    toc;
    
    %get p_value for each ROI
    pVal = sum(logical(placeField_nb_shuffle))/options.nbShuffle;
    
    %assign which ROIs have sig fields
    Place_cell{ii}.placeField.sig_ROI = pVal < 0.05;
    end
    
end
toc;

%% Centroid difference/shift (Danielson 2016)
%use this
%calculate the angle between the the tuning specificity vectors (not
%difference)

%whether to plot figure
options.plotFigure = 0;

%centroid threshold - how much angular separation between tuning vectors
options.centroidThres = pi/10; %(pi/5 - 20 cm); pi/10 - 10 cm

centroids = centroidDiffV2_newEvents(Place_cell,options);


%% Event spiral and mean dF/F by lap plots - Zaremba 2017 based
%modify this for single trials

% plots each event for given ROI during running epochs along lap across laps for

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

ROIrange = 242;


spiralEvents = spiral_eventsV6_newSplit_PFtest(Behavior_split_lap, Events_split_lap, Behavior_split, Events_split, Imaging_split, Place_cell, Behavior, ROIrange,options);


%% Plot Spatial tuning curves (STC) with selected ROI, PV correlation, TC correlation, PF distribution
%need input from ROI selector function

%spatial tuning curves
%which trials to sort by
options.sortOrder = 1;
plotSTC(selectROI_PF_idx,Place_cell,options)

%place field distributions, place field counts, and PF distributions
plotPF(Place_cell)

%spatial correlation function goes here
[PVcorr,TCcorr] = spatialCorr(Place_cell);

%}

