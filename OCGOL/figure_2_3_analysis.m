%% Import variables and define options

%lab workstation
%input directories to matching function
%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_120618_1'}; % field rate error
%path_dir = {'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I42R_AB_d1_032118_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_2'};
%path_dir = {'G:\Figure_2_3_selective_remap\I53LT_AB_sal_113018_1'}; %place field finder problem - adjust
%path_dir = {'G:\Figure_2_3_selective_remap\I56_RTLS_AB_prePost_sal_042419_1'}; %place field finder problem - adjust
%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_113018_1'};
path_dir = {'G:\Figure_2_3_selective_remap\I57_RTLS_AB_prePost_792_042519_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I45_RT_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I46_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I57_LT_ABrand_no_punish_042119_1'};

%load place cell variables for each session
%get mat directories in each output folder
for ii=1:size(path_dir,2)
    %get matfile names for each session
    matfiles{ii} = dir([path_dir{ii},'\output','\*.mat']);
    %input directories (.mat CNMF, .XML, .CSV)
    inputfiles{ii} = dir([path_dir{ii},'\input','\*.mat']);
    %removed ROI directory
    removeROIfiles{ii} = dir([path_dir{ii},'\removedROI','\*.mat']);
end

%load in place cell variables (and others later)
for ii = 1:size(path_dir,2)
    %add event variables
    session_vars{ii} = load(fullfile(matfiles{ii}.folder,matfiles{ii}.name),'Place_cell', 'Behavior',...
        'Behavior_split_lap','Behavior_split','Events_split','Events_split_lap', 'Imaging_split');
    %load in relevant CNMF variables for visualization purposes
    CNMF_vars{ii} = load(fullfile(inputfiles{ii}.folder,inputfiles{ii}.name),'A_keep','C_keep','Cn','Coor_kp','expDffMedZeroed','dims');
    %load in template (stack average)
    templates{ii} = load(fullfile(path_dir{ii},'template_nr.mat'),'template');
    %load in logical vector with selected and rejected ROIs
    removeROI{ii} = load(fullfile(removeROIfiles{ii}.folder,removeROIfiles{ii}.name),'compSelect');
end

%% Define tuned logical vectors
%flag to all A or B trial or only correct A or B trials
options.allCorrect = 1; %1  = A correct; 2 = B correct
%returns struct of structs
[tunedLogical] = defineTunedLogicals(session_vars,options);

%% Add a filter for logical selection of A/B selective neurons here (Figure 2)

%% Add a filter for logical selection of A&B remapping neurons here (Figure 3)

%% Plot fraction of each neuron tuned 

%plot pie chart for each and return counts
[tuned_fractions] = fractionTuned(tunedLogical);

%export for cumulative analysis
%make cumul_analysis folder
mkdir(path_dir{1},'cumul_analysis')
%save the fractions output data
save(fullfile(path_dir{1},'cumul_analysis','frac_tuned.mat'),'tuned_fractions');

%% Find  place fields
%use rate map - number of event onsets/ occupancy across all laps
options.gSigma = 3;
%which place cell struct to do placefield extraction on
%iterate through place_cell cells of interest
%4 - all A regardless if correct
%5 - all B regardless if correct

%for each session
for ss=1:size(session_vars,2)
    %for ii =[4,5] %all A or B 
    for ii =[1 2 ] %only correct A or B
        %works for 1 (A) laps
        %373 A trial not merged
        %449 - third field seems not caught
        %clipping vectors/cells at end
        options.place_struct_nb = ii;
        [session_vars{ss}.Place_cell] = place_field_finder_gaussian(session_vars{ss}.Place_cell,options);
    end
end

%for overlapping area merge (calculate number of place fields)
%https://www.mathworks.com/matlabcentral/answers/361760-area-between-two-overlapping-plots

%% Comparison of AUC/min rate of exclusive A tuned or exclusive B tuned neurons

options.tuning_criterion = 'si'; %si or ts
[total_AUC_min] = AUC_scatter(tunedLogical,session_vars,options);
save(fullfile(path_dir{1},'cumul_analysis','auc.mat'),'total_AUC_min');


%% Number of place fields and widths for each sub-class of neurons

options.tuning_criterion = 'si'; %si or ts
[placeField_dist] = placeField_properties(session_vars, tunedLogical,options);
%save the place field distributions output data
save(fullfile(path_dir{1},'cumul_analysis','placeField_dist.mat'),'placeField_dist');

%% Calculate the transient rates in each of the place fields (integrate later) and recalculate centroids based on highest transient rate field

%funtion to calculate transient rate in field
%take raw event rate and divide by occupancy (s) transients/s
% TODO: see ifrecalculate to see if dividing by normalized occupancy (fractional 0-1)
%yields different result

%which sessions to use to calculate the in field transient rate
%[1 2] - only correct A B trials
%[4 5] - all A B trials
options.selectSes = [1 2];

[field_event_rates,pf_vector] = transient_rate_in_field_single_ses(session_vars,options);

%% Plot spatial tuning curves according to transient rate and return max bin position (split into 2 fxn later)
%normalized to each neuron across 

options.tuning_criterion = 'si'; %si or ts
%AandB, AorB, onlyA, onlyB, neither all
%all won't work b/c no defined fields to sort by
options.trialTuning = 'onlyB';
options.selectSes = [1 2];
%sort according to which trial 1 2 4 5
options.sortTrial = 2;
[max_bin_rate] = plot_STC_transient_rate_single_ses(session_vars,tunedLogical,field_event_rates, pf_vector,options);

%% PV and TC correlation matrices for each class of tuned neurons

options.tuning_criterion = 'ts';
[correlation] = PV_TC_correlation_single_ses(session_vars,tunedLogical,options);

%save the fractions output data
save(fullfile(path_dir{1},'cumul_analysis','corr.mat'),'correlation');

%% Centroid distribution across lap for A tuned and B tuned neurons
%use tuning spec criterion for this
options.tuning_criterion = 'ts'; %si or ts
[centroid_ct] = centroid_dist(tunedLogical, max_bin_rate,options);

%save the fractions output data
save(fullfile(path_dir{1},'cumul_analysis','centroid.mat'),'centroid_ct');


%% Generate STC maps of neurons tuned in either session and plot side by side
%customize to add options
%tuned in both sessions by SI score
%sorted by A trials

%set option as to how to select neurons for plots
options.tuning_criterion = 'si'; %si or ts
%normalized across both sessions

%PV correlation embedded inside of this function
plot_STC_OCGOL_singleSes(session_vars,tunedLogical,options)


%% Get percentage correct in each trial type and

trialOrder = session_vars{1}.Behavior.performance.trialOrder;
trialCorrect = session_vars{1}.Behavior.performance.trialCorrect;

fracA_corr = size(find(trialOrder == 2),1) / size(find(trialOrder == 2 | trialOrder == 20),1)

fracB_corr = size(find(trialOrder == 3),1) / size(find(trialOrder == 3 | trialOrder == 30),1)


%% Show outlines of selected and discarded neurons over FOV
%only 1 session
ii=1;
%get idx's of selected and rejected ROIs
selectedROI_idx = find(removeROI{ii}.compSelect == 1)';
rejectedROI_idx = find(removeROI{ii}.compSelect == 0)';

%plot BW outline of the component
figure
imagesc(templates{ii}.template);
hold on
axes(gca);
axis square
xticks(gca,[])
yticks(gca,[])
grayMap = brighten(gray,0.6);
colormap(gca,grayMap)

%plot all selected ROIs as green
for ROI = selectedROI_idx
    %plot componenet outline
    plot(CNMF_vars{ii}.Coor_kp{ROI}(1,:),CNMF_vars{ii}.Coor_kp{ROI}(2,:),'g', 'LineWidth',1);
    %pause(0.01)
end

%plot all selected ROIs as green
for ROI = rejectedROI_idx
    %plot componenet outline
    plot(CNMF_vars{ii}.Coor_kp{ROI}(1,:),CNMF_vars{ii}.Coor_kp{ROI}(2,:),'r', 'LineWidth',1);
    %pause(0.01)
end







