%% Import variables and define options

%lab workstation
%input directories to matching function
%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_120618_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I42R_AB_d1_032118_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_2'};
%path_dir = {'G:\Figure_2_3_selective_remap\I53LT_AB_sal_113018_1'}; 
%path_dir = {'G:\Figure_2_3_selective_remap\I56_RTLS_AB_prePost_sal_042419_1'}; %bug with remapper select ROI code
%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_113018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I57_RTLS_AB_prePost_792_042519_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I45_RT_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I46_AB_d1_062018_1'};
path_dir = {'G:\Figure_2_3_selective_remap\I57_LT_ABrand_no_punish_042119_1'};

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
        'Behavior_split_lap','Behavior_split','Events_split','Events_split_lap', 'Imaging_split','updated_dff');
    %load in relevant CNMF variables for visualization purposes
    CNMF_vars{ii} = load(fullfile(inputfiles{ii}.folder,inputfiles{ii}.name),'A_keep','C_keep','Cn','Coor_kp','expDffMedZeroed','dims');
    %load in template (stack average)
    templates{ii} = load(fullfile(path_dir{ii},'template_nr.mat'),'template');
    %load in logical vector with selected and rejected ROIs
    removeROI{ii} = load(fullfile(removeROIfiles{ii}.folder,removeROIfiles{ii}.name),'compSelect');
end

%% Define tuned logical vectors

%flag to all A or B trial or only correct A or B trials
options.allCorrect = 1;
%select which session to use (single session for these animals)
options.sessionSelect = [1];
%returns struct of structs
[tunedLogical] = defineTunedLogicals(session_vars,options);

%% Plot example neuron (for display/figure) - trace, position, event spiral, TS arrows, A vs. B
%I42L example plots ROIs
%59, 66 (nice), 111(nice),162,180, 208,214 234 (nice),252,280  - A specific
%22,34, 46 (one used previously),62,207,246 (nice),282 (nice) - B specific
%10,30, 49, 82, 89 (nice), 104 (nice),158,161 (nice), 167,171,182,197 (used before) -common
%237
%currently displaying in Figure 2
%[234,246,197] 
%break down in simpler code in the future;

%generate figure used in task selecitve figure
options.plotFigure2 = 1;
raster_spiral_single_ses(session_vars,CNMF_vars,removeROI,templates,options)

%% Plot fraction of each neuron tuned 
%TODO: add filter for classfing whether each field is significant (min 5 events)

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
    for ii =[1 2] %only correct A or B
        %works for 1 (A) laps
        %373 A trial not merged
        %449 - third field seems not caught
        %clipping vectors/cells at end
        options.place_struct_nb = ii;
        disp(['Running trial type: ', num2str(ii)]);
        [session_vars{ss}.Place_cell] = place_field_finder_gaussian(session_vars{ss}.Place_cell,options);
    end
end

%parse place fields significance by checking that a minimum of 5 events
%occured in the field

%for overlapping area merge (calculate number of place fields)
%https://www.mathworks.com/matlabcentral/answers/361760-area-between-two-overlapping-plots

%% Calculate the transient rates in each of the place fields (integrate later) and recalculate centroids based on highest transient rate field

%funtion to calculate transient rate in field
%take raw event rate and divide by occupancy (s) transients/s
% TODO: see ifrecalculate to see if dividing by normalized occupancy (fractional 0-1)
%yields different result

options.sessionSelect = [1];
%A correct/B correct or all
options.selectTrial = [1 2];
%continue to modify 
[field_event_rates,pf_vector,field_total_events, select_fields] = transient_rate_in_field_multi_ses(session_vars,options);

%% Plot spatial tuning curves according to transient rate and return max bin position (split into 2 fxn later)
%normalized to each neuron across 

options.tuning_criterion = 'si'; %si or ts
%AandB, AorB, onlyA, onlyB, neither all
%all won't work b/c no defined fields to sort by
options.trialTuning = 'onlyA';
options.selectSes = [1 2];
%sort according to which trial 1 2 4 5
options.sortTrial = 1;

%outputs:
%max_transient_peak - peak idx where event rate is highest for correct A/B
%trials
[max_bin_rate,max_transient_peak] = plot_STC_transient_rate_single_ses(session_vars,tunedLogical,field_event_rates, pf_vector,options);

%% Calculate centroid difference between A&B tuned neurons (max in field transient rate)

options.tuning_criterion = 'ts';
[cent_diff_AandB, pf_vector_max] = centroid_diff_single_ses(session_vars,tunedLogical, pf_vector,field_event_rates,options);

%save the fractions output data
save(fullfile(path_dir{1},'cumul_analysis','centroid_diff.mat'),'cent_diff_AandB');

%% Split neurons by A or B task selective category - A or B selective (exclusive)
%which criterion to use for task-selective ROIs
%ts or both - ts selects only selective neurons based on TS tuning
%criterion
%both - uses both SI and TS criterion to select selectiven neurons
options.tuning_criterion = 'ts';
%display events vs position for each task selective neuron in A or B
options.dispFigure = 0;
[task_selective_ROIs] = task_selective_categorize(cent_diff_AandB, tunedLogical, pf_vector_max, session_vars, max_transient_peak,options);


%% Number of place fields and widths for each sub-class of neurons
%add filter for classfing whether each field is significant (min 5 events)

options.tuning_criterion = 'si'; %si or ts
%A correct/B correct or all
options.selectTrial = [1 2];

[placeField_dist] = placeField_properties(session_vars,tunedLogical,select_fields,task_selective_ROIs,options);
%save the place field distributions output data
save(fullfile(path_dir{1},'cumul_analysis','placeField_dist.mat'),'placeField_dist');


%% Split A&B neurons by remapping category - common, partial, global, rate remapping
%which criterion to use for task-selective ROIs
options.tuning_criterion = 'ts';
%display events vs position for each task selective neuron in A or B
options.dispFigure = 0;
%make sure that this function does not overwrite the the previous
%task_selective_ROIs structure
[task_remapping_ROIs] = remapping_categorize(cent_diff_AandB, tunedLogical, pf_vector_max, session_vars, max_transient_peak,options);

%% PV and TC correlation matrices for each class of tuned neurons

options.tuning_criterion = 'ts';
[correlation] = PV_TC_correlation_single_ses(session_vars,tunedLogical,task_selective_ROIs,options);

%save the fractions output data
save(fullfile(path_dir{1},'cumul_analysis','corr.mat'),'correlation');


%% Centroid distribution across lap for A tuned and B tuned neurons
%use tuning spec criterion for this
options.tuning_criterion = 'selective_filtered'; %si or ts or selective_filtered
[centroid_ct] = centroid_dist_task_selective(tunedLogical,task_selective_ROIs, max_bin_rate,options);

%save the fractions output data
save(fullfile(path_dir{1},'cumul_analysis','centroid.mat'),'centroid_ct');

%% Generate STC maps of neurons tuned in either session and plot side by side
%customize to add options
%tuned in both sessions by SI score
%sorted by A trials

%set option as to how to select neurons for plots
options.tuning_criterion = 'selective_filtered'; %si or ts or selective_filtered
%normalized across both sessions

plot_STC_OCGOL_singleSes_task_selective(session_vars,tunedLogical,task_selective_ROIs,options);

%% Comparison of AUC/min rate of exclusive A tuned or exclusive B tuned neurons

options.tuning_criterion = 'ts'; %si or ts
[total_AUC_min] = AUC_scatter(tunedLogical,task_selective_ROIs,session_vars,options);
save(fullfile(path_dir{1},'cumul_analysis','auc.mat'),'total_AUC_min');


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

%% Synchronized calcium event analysis




