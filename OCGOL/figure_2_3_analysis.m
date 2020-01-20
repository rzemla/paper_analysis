function figure_2_3_analysis(path_dir)
%% Set to not display figures when running in serial export fashion

set(0,'DefaultFigureVisible','off'); 

%% Import variables and define options

%input directories to matching function
%this one was excluded b/c analyzed after PSEM silecing occurred

%include in processing (detailed notes below for each dataset)
%1
%path_dir = {'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1'};
%2
%path_dir = {'G:\Figure_2_3_selective_remap\I42R_AB_d1_032118_1'};
%3
%path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1'};
%4
%path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_2'};
%5 
%path_dir = {'G:\Figure_2_3_selective_remap\I53LT_AB_sal_113018_1'};
%6
%path_dir = {'G:\Figure_2_3_selective_remap\I56_RTLS_AB_prePost_sal_042419_1'};
%7
%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_113018_1'};  
%8
%path_dir = {'G:\Figure_2_3_selective_remap\I57_RTLS_AB_prePost_792_042519_1'};
%9
%path_dir = {'G:\Figure_2_3_selective_remap\I45_RT_AB_d1_062018_1'};
%10
%path_dir = {'G:\Figure_2_3_selective_remap\I46_AB_d1_062018_1'};
%11
%path_dir = {'G:\Figure_2_3_selective_remap\I57_LT_ABrand_no_punish_042119_1'};

%DON'T INCLUDE (only one)
%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_120618_1'};
 %near flat PV until end - check PSAM experiment order - late experiment - animal
%already exposed to PSEM - use I52RT_AB_sal_113018 (no PSEM exposure yet
%during task)- processed below

%I52_RT - I52RT_ABp_112218 - last imaged training day

% path_dir = {'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1'};
% path_dir = {'G:\Figure_2_3_selective_remap\I42R_AB_d1_032118_1'};
% path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1'};
% %path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_2'};

% %path_dir = {'G:\Figure_2_3_selective_remap\I53LT_AB_sal_113018_1'}; 
% %good -%no PSEM exposure

% path_dir = {'G:\Figure_2_3_selective_remap\I56_RTLS_AB_prePost_sal_042419_1'}; %near flat PV until end %bug with idx remapper code; run parser (global) global difference
%good first PSEM sliencing exp; 

%I56_RTLS_ABrand_punish_041719 - last exp - part of learning before silence

% path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_113018_1'};  
% %near flat PV until end %bug with idx remapper code; run parser (global) global difference
% 
% path_dir = {'G:\Figure_2_3_selective_remap\I57_RTLS_AB_prePost_792_042519_1'};
%near flat PV until end %bug with idx remapper code; run parser (global) global difference
%OK acquired the day of silencing but before; can also use session before
%and from learning - late learning session in Figure 4 learning datasets

% path_dir = {'G:\Figure_2_3_selective_remap\I45_RT_AB_d1_062018_1'};
% 
% path_dir = {'G:\Figure_2_3_selective_remap\I46_AB_d1_062018_1'};
% 
% path_dir = {'G:\Figure_2_3_selective_remap\I57_LT_ABrand_no_punish_042119_1'}; %OK - well trained - from one of learning days
%I57_LT_AB_prePost_sal_050619 - last session before next day silencing with PSAM 

%whether to load place field data processed below
options.loadPlaceField_data = 1;

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
        'Behavior_split_lap','Behavior_split','Events_split','Events_split_lap', 'Imaging_split','updated_dff','Imaging');
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

%skip plots
if 0
    %% Plot example neuron (for display/figure) - trace, position, event spiral, TS arrows, A vs. B
    %I42L example plots ROIs
    %59, 66 (nice), 111(nice),162,180, 208,214 234 (nice),252,280  - A specific
    %22,34, 46 (one used previously),62,207,246 (nice),282 (nice) - B specific
    %10,30, 49, 82, 89 (nice), 104 (nice),158,161 (nice), 167,171,182,197 (used before) -common
    %237
    %currently displaying in Figure 2
    %[234,246,197]
    
    
    %generate figure used in task selecitve figure
    options.plotFigure2 = 1;
    raster_spiral_single_ses(session_vars,CNMF_vars,removeROI,templates,options)
    
    %% Generate partial remapper display for supplement figure
    
    options.plotFigure2 = 1;
    options.animalIdx = 8;
    options.partial = 1;
    raster_spiral_partial_maps(session_vars,CNMF_vars,removeROI,templates,partial_idx_by_animal_zone,options)    
    
    %% Preprare inputs for raster spiral plotter
    %add option above to prevent from re-calculating values on every run
    tic;
    [plot_raster_vars] = prepare_inputs_raster(session_vars,CNMF_vars,removeROI);
    toc;
    
    %% Plot raster here only
    %take the processed inputs for spiral plots, ROI zooms above and plot for
    %figure presentation
    options.plotFigure2 = 1;
    plot_raster_spiral_only(plot_raster_vars,session_vars,templates,task_remapping_ROIs,path_dir,options)
    
end

%% Get mean norm position of each reward (A or B)
[reward_positions] = extract_reward_position(session_vars);

%save the fractions output data %
save(fullfile(path_dir{1},'cumul_analysis','reward_positions.mat'),'reward_positions');

%% Export rate maps from all ROIs (visualization supplement figures)

[rate_maps_display] = export_rate_maps_for_display(session_vars);

%save the fractions output data %
save(fullfile(path_dir{1},'cumul_analysis','rate_maps_display.mat'),'rate_maps_display');


%% Find place fields
if options.loadPlaceField_data == 0
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
    
    %save whole place cell struct and load in and replace for each session in
    %the future
    %make post-processing directory (postProcess)
    mkdir(path_dir{1},'postProcess')
    %for each Place_cell session extract placeField struct
    %use trial types here
    for tt=1:2
        session_pf(tt).placeField = session_vars{1}.Place_cell{tt}.placeField;
    end
    
    %save Place_cell struct in that directory
    save(fullfile(path_dir{1},'postProcess','placeField_upd_struct.mat'),'session_pf')
else
    tic;
    disp('Loading place field data')
    load(fullfile(path_dir{1},'postProcess','placeField_upd_struct.mat'));
    toc
    %replace the Place_cell struct in the session_vars cell
    for tt=1:2
        session_vars{1}.Place_cell{tt}.placeField = session_pf(tt).placeField;
    end
end


%for overlapping area merge (calculate number of place fields)
%https://www.mathworks.com/matlabcentral/answers/361760-area-between-two-overlapping-plots

%% Extract speed according to selected trial types
%(1) only A corr
%(2) only B corr
%(3) all laps
%(4) all A trials
%(5) all B trials

%extract speed on each trial and load into session_vars struct
session_vars = extract_speed_for_each_trial_set(session_vars);

%% Calculate the transient rates in each of the place fields (integrate later) and recalculate centroids based on highest transient rate field

%funtion to calculate transient rate in field
%take raw event rate and divide by occupancy (s) transients/s

options.sessionSelect = [1];
%A correct/B correct or all
options.selectTrial = [1 2];
%continue to modify 
[field_event_rates,pf_vector,field_total_events, select_fields] = transient_rate_in_field_multi_ses(session_vars,options);

%% Get max transient peak here (developed from multi-ses analysis)
%get field event rates of max peak

%select whether to use TS vector of adjusted vector for cells with single fields
%1 - adjusts based on max transient field and corresponding pf_vector (for
%pf with single fields)
%0 - uses the tuning vector from TS calculation (single fields only)
options.select_adj_vec = 1;

[max_bin_rate,max_transient_peak] = max_transient_rate_multi_ses(session_vars,field_event_rates,pf_vector,options);

%% Calculate centroid difference between A&B tuned neurons (max in field transient rate) - alt - take centroid of PF of single PF place cells
%returns same variables as above
options.tuning_criterion = 'ts';
[cent_diff,cent_diff_AandB, pf_vector_max] = centroid_diff_single_ses_single_PF(session_vars,tunedLogical, pf_vector,field_event_rates,select_fields,options);

%save the fractions output data %
save(fullfile(path_dir{1},'cumul_analysis','centroid_diff.mat'),'cent_diff_AandB');
%save the cent diff for all neurons (TS) tuned
save(fullfile(path_dir{1},'cumul_analysis','centroid_diff_all.mat'),'cent_diff');

%% Split neurons by A or B task selective category - A or B selective (exclusive)
%which criterion to use for task-selective ROIs
%ts or both - ts selects only selective neurons based on TS tuning
%criterion
%both - uses both SI and TS criterion to select selectiven neurons
options.tuning_criterion = 'both';
%display events vs position for each task selective neuron in A or B
options.dispFigure = 0;
%[task_selective_ROIs] = task_selective_categorize(cent_diff_AandB, tunedLogical, pf_vector_max, session_vars, max_transient_peak,options);
%yields same output - correct edge
%QC checked
[task_selective_ROIs_ses] = task_selective_categorize_multi_ses(tunedLogical,session_vars, max_transient_peak,options);

%convert task_selective ROIs output to remove cell session
task_selective_ROIs.A = task_selective_ROIs_ses{1}.A;
task_selective_ROIs.B = task_selective_ROIs_ses{1}.B;

%% Generate STC maps - task-selective neurons - STC Figure 2 generator
%customize to add options
%tuned in both sessions by SI score
%sorted by A trials

%export A/B trial-normalized STCs for task-selective neurons (using both option above - selective by S.I. or T.S criteria)

%set option as to how to select neurons for plots
options.tuning_criterion = 'selective_filtered'; %si or ts or selective_filtered

%normalized across both sessions
[task_sel_STC] = plot_STC_OCGOL_singleSes_task_selective(session_vars,tunedLogical,task_selective_ROIs,options);

%save task-selective STCs for cumulative Figure 2 plots
save(fullfile(path_dir{1},'cumul_analysis','task_sel_STC.mat'),'task_sel_STC');

%% Return significant place field logical
%on each trial type, returns logical where at least 1 sig place field was
%detected
options.selectTrial = [1 2];
[pf_count_filtered_log] = place_field_logical(select_fields,options);

%% Filtered tuning classes (same output as that from fraction tuned, but needed for place field script below)
%QC checked
[ROI_idx_tuning_class] = return_tuned_classes(tunedLogical,pf_count_filtered_log);

%% Return all A&B tuned neurons (by SI or TS category tuned)

Atuned_all_si = union(ROI_idx_tuning_class.si.Aonly, ROI_idx_tuning_class.si.AB);
Btuned_all_si = union(ROI_idx_tuning_class.si.Bonly, ROI_idx_tuning_class.si.AB);

Atuned_all_ts = union(ROI_idx_tuning_class.ts.Aonly, ROI_idx_tuning_class.ts.AB);
Btuned_all_ts = union(ROI_idx_tuning_class.ts.Bonly, ROI_idx_tuning_class.ts.AB);

Atuned_all_si_ts = union(Atuned_all_si,Atuned_all_ts);
Btuned_all_si_ts = union(Btuned_all_si,Btuned_all_ts);

%A&B tuned neurons by either SI or TS criteria - for determining
%unclassified class in remapping neuron category
ABtuned_all_si_ts = intersect(Atuned_all_si_ts, Btuned_all_si_ts);

%save A&B tuned neurons to work out fraction of category
save(fullfile(path_dir{1},'cumul_analysis','AB_tuned_si_ts.mat'),'ABtuned_all_si_ts');

%% Number of place fields and widths for each sub-class of neurons
%add filter for classfing whether each field is significant (min 5 events)

%no tuning criterion - return parameters for both tuning params
%options.tuning_criterion = 'si'; %si or ts
%A correct/B correct or all
options.selectTrial = [1 2];

[placeField_dist, pf_count_filtered_log, pf_count_filtered] = placeField_properties(session_vars,tunedLogical,select_fields,task_selective_ROIs,ROI_idx_tuning_class,options);
%save the place field distributions output data
save(fullfile(path_dir{1},'cumul_analysis','placeField_dist.mat'),'placeField_dist');

%% Plot fraction of each neuron tuned 

%plot pie chart for each and return counts
[tuned_fractions,~] = fractionTuned(tunedLogical,pf_count_filtered_log);

%export for cumulative analysis
%make cumul_analysis folder
mkdir(path_dir{1},'cumul_analysis')
%save the fractions output data
save(fullfile(path_dir{1},'cumul_analysis','frac_tuned.mat'),'tuned_fractions');

%% Export S.I/T.S scores for each neuron tuned by each criterion as well as score for A_sel/B_sel neurons

[tuning_scores] = si_ts_score_category(session_vars,ROI_idx_tuning_class,task_selective_ROIs,options);

%save task-selective STCs for cumulative Figure 2 plots
save(fullfile(path_dir{1},'cumul_analysis','tuning_scores.mat'),'tuning_scores');

%% PV and TC correlation matrices for each class of tuned neurons

%correlation values (1); spatial tuning curves (2)
[correlation,STC_export] = PV_TC_correlation_single_ses(session_vars,tunedLogical,task_selective_ROIs,ROI_idx_tuning_class,options);

%save the tuning correlation output
save(fullfile(path_dir{1},'cumul_analysis','corr.mat'),'correlation');

%save spatial tuning curve output for each ROI
save(fullfile(path_dir{1},'cumul_analysis','STC.mat'),'STC_export')

%% Centroid distribution across lap for A/B selective tuned neurons (only - modify inputs in future for rest of neurons if necessary)
%QC checked

%use tuning spec criterion for this
options.tuning_criterion = 'selective_filtered'; %si or ts or selective_filtered
[centroid_ct,centroid_bins] = centroid_dist_task_selective(tunedLogical,task_selective_ROIs, max_bin_rate,options);

%save the fractions output data
save(fullfile(path_dir{1},'cumul_analysis','centroid.mat'),'centroid_ct','centroid_bins');

%% Extract AUC/min rate of selective A/B tuned neurons (and other classes)
%RUN, NO RUN FOR NOW for A sel, Bsel, A&B (si/ts tuned)
[total_AUC_min] = AUC_scatter(tunedLogical,task_selective_ROIs,session_vars,ROI_idx_tuning_class,options);

save(fullfile(path_dir{1},'cumul_analysis','auc.mat'),'total_AUC_min');


%% Split A&B neurons by remapping category - common, partial, global, rate remapping
%which criterion to use for task-selective ROIs
options.tuning_criterion = 'ts';
%display events vs position for each task selective neuron in A or B
options.dispFigure = 0;
%number of degrees of centroid difference (input these into analysis below)
%45 deg ~25 cm; 
%36 deg ~20 cm;
%27 deg ~15 cm;
%18 deg ~10 cm;
%25 deg ~13 cm;
options.deg_thres = 18;
%ranges for splitting the global remappers
%0-10 cm; 10 - 30cm; 30+ cm
options.deg_ranges = [0 18 54];
%degree threshold for partial remappers
options.partial_deg_thres = [25 25];
%choice between KS test of unpaired Mann Whitney U (later)
%either 'ranksum' or ks
options.AUC_test = 'ranksum';
%significance level of test
options.p_sig = 0.05;
%make sure that this function does not overwrite the the previous
%task_selective_ROIs structure

[task_remapping_ROIs,partial_field_idx,event_5_min_and_occup_filtered_ROI] = remapping_categorize(cent_diff, tunedLogical, pf_vector_max ,pf_vector, session_vars,...
                        max_transient_peak,pf_count_filtered_log, pf_count_filtered,select_fields,options);

%ROIs that are single fields by TS; min 5 events; 80% occupnancy on opposing trials 
event_5_min_and_occup_filtered_ROI;
                                        
%% Split remapping categories by stat sig of rate map correlations (global vs non global)

[remapping_corr_idx,tun_curve_corr] = remapping_correlations(session_vars,tunedLogical,task_selective_ROIs,ROI_idx_tuning_class, task_remapping_ROIs,pf_count_filtered, options);

%tun_curve_corr - contains the r and p values for the rate map correlations
%of each neuron

%export the r values and p values associated with each rate map correlation
save(fullfile(path_dir{1},'cumul_analysis','tun_curve_corr.mat'),'tun_curve_corr'); 

%% Speed data for each lap (extract speed in each bin) - insert speed data

[mean_bin_speed, lap_bin_split] =task_sel_speed(tunedLogical,task_selective_ROIs,session_vars,ROI_idx_tuning_class,options);


%% Compare speed and AUC of events for rate remapping neurons and DEFINE  remapping each category
%determine the other subcategories of the remapping neurons
%run 2 way anova for each neuron to determine which one is rate remapper

[remapping_corr_idx,remap_idx_traces,com_idx_traces,AUC_remappers] = speed_AUC_comparison(task_remapping_ROIs, remapping_corr_idx, lap_bin_split, session_vars, max_transient_peak,...
    STC_export, event_5_min_and_occup_filtered_ROI,ABtuned_all_si_ts, options);

%get rate remapping neurons as well with 2-way ANOVA
%save the neurons parsed by correlation remapping criteria
save(fullfile(path_dir{1},'cumul_analysis','remap_corr_idx.mat'),'remapping_corr_idx'); 

%export remapping traces for supp data figure
save(fullfile(path_dir{1},'cumul_analysis','remap_traces.mat'),'remap_idx_traces','com_idx_traces'); 

%export the AUC values for analysis
save(fullfile(path_dir{1},'cumul_analysis','AUC_remapping.mat'),'AUC_remappers'); 


%% Event vs. speed analysis

[mean_event_speed] = event_vs_speed(session_vars, task_selective_ROIs,ROI_idx_tuning_class,...
                select_fields,max_transient_peak,mean_bin_speed,lap_bin_split,options);

%export the lap speed data and event speed data for cumulative analysis
save(fullfile(path_dir{1},'cumul_analysis','lap_and_event_speed.mat'),'mean_bin_speed', 'lap_bin_split','mean_event_speed');

%% Calculate centroid difference between neurons with single place fields in each A and B trials 
%max_transient_peak - index of max transient peak
%pf_count_filtered - stat. significant place fields for each neuron 
%cent_diff - between max fields
%ROI_idx_tuning_class - category of each ROI (filtered by events and place
%field)

[ts_bin_conv_diff,pf_distance_metric_ts,final_global_remap_bin_diff] = place_field_AandB_distribution(session_vars,ROI_idx_tuning_class,remapping_corr_idx,cent_diff,pf_count_filtered,max_transient_peak);

%export distances
save(fullfile(path_dir{1},'cumul_analysis','place_field_AB_distances.mat'),'ts_bin_conv_diff','pf_distance_metric_ts',...
    'final_global_remap_bin_diff');

%% Calculate centroid difference/distasnce between neurons with common place fields for setting
 
[common_bin_conv_diff,common_pf_distance_metric] = place_field_AandB_distribution_common_neurons(session_vars,ROI_idx_tuning_class,remapping_corr_idx,cent_diff,pf_count_filtered,max_transient_peak);

%export distances for common neurons defined as having similar
%correlation maps 
save(fullfile(path_dir{1},'cumul_analysis','place_field_common_distances.mat'),'common_bin_conv_diff','common_pf_distance_metric');

%% Generate STC maps of neurons tuned in either session and plot side by side
%customize to add options
%tuned in both sessions by SI score
%sorted by A trials

%add option here to switch between the display of different categories of
%cells

%set option as to how to select neurons for plots
options.tuning_criterion = 'remapping_filtered'; %si or ts or selective_filtered
%normalized across both sessions

plot_STC_OCGOL_singleSes_task_remapping(session_vars,tunedLogical,task_remapping_ROIs,path_dir,options);


%% Export dF/F values for rate remapping neurons

[rate_mean_dFF] =export_dFF_rasters(session_vars,remapping_corr_idx);

save(fullfile(path_dir{1},'cumul_analysis','rate_mean_dFF.mat'),'rate_mean_dFF');

%% Extract spatial bins for each class of remapping neurons
%add extraction for ROIs gathered with correlation based remapping

bin_center = extract_centroid_bins_remap(cent_diff,task_remapping_ROIs, remapping_corr_idx, select_fields,partial_field_idx);

save(fullfile(path_dir{1},'cumul_analysis','place_field_centers_remap.mat'),'bin_center');

%plot remapping as scatter
%green dots - common ; red - partial field 

%distribution of partial fields
figure
hold on
histogram(bin_center.partial_far(1,:),0:50:100)
histogram(bin_center.partial_far(2,:),0:50:100)

%% Save the idxs of each class of remapping neurons, placeFields that match min 5 event criteria and tuned logical

save(fullfile(path_dir{1},'cumul_analysis','select_ROI_criteria.mat'),'select_fields','tunedLogical',...
                'task_remapping_ROIs','task_selective_ROIs');

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

%% Show outlines of only select neurons for visulization purposes

if 0
ii=1;
%get idx's of selected and rejected ROIs

%filter coor kp for selected ROIs
Coor_kp_soma_rem = CNMF_vars{ii}.Coor_kp(removeROI{ii}.compSelect);

%show only select neurons for visualization purposes
selectVis =1;

%show only select neurons for localization purposes
selectedROI_idx =[627 416 269 524 458];

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
    plot(Coor_kp_soma_rem{ROI}(1,:),Coor_kp_soma_rem{ROI}(2,:),'g', 'LineWidth',1);
    pause()
    
end
end

%% Extract frame indices for figure 2A/3A

%skip
if 0
% MATLAB:
% 1) split motion corrected hdf stack in A and B trials
% 2) downsampled 5x timewise (interpolate /average over interval)
% 3) convert to uint16 type (from double)
% 4) save as tif stacks
% Imagej/Fiji:
% 1) Load downsample tif stacks
% 2) Get stack average and max projection for each
% 3) Subtract average from max projectionf or each
% 4) Apply blue LUT for A and red LUT for A
% 5) Window/level both until estimated background looks near black
% 6) Window/level stack average (in gray LUT) to get nice/clear contrast
% 7) Use merge colors feature with gray channel being Average stack, red is B, and blue is A.
% 8) Convert to 8bit and save as tif.

%get absolute frames of correct A and and correct B trials (no run
%intervals included)
[~,frames_A,~] = intersect(session_vars{1}.Imaging.time ,session_vars{1}.Imaging_split{1}.time_restricted);
[~,frames_B,~] = intersect(session_vars{1}.Imaging.time ,session_vars{1}.Imaging_split{2}.time_restricted);

%extract logical for selecting run frames only
logical_run_A = session_vars{1}.Behavior_split{1}.run_ones == 1;
logical_run_B = session_vars{1}.Behavior_split{2}.run_ones == 1;

%indices of run frames only
frames_A_run = frames_A(logical_run_A);
frames_B_run = frames_B(logical_run_B);
%indices of no run frames
frames_A_norun = frames_A(~logical_run_A);
frames_B_norun = frames_B(~logical_run_B);

%read in the motion correct stack
Y = read_file(fullfile(path_dir{1}, '1_nr.h5'));

%if 1, then only run frames
%if 2, then no run, 
%if 3 then all
runOnly = 2;

%select run only or all
if runOnly == 1
    %split Y in into A and B frames
    Ya = Y(:,:,frames_A_run);
    Yb = Y(:,:,frames_B_run);
elseif runOnly == 2
    Ya = Y(:,:,frames_A_norun);
    Yb = Y(:,:,frames_B_norun);  
elseif runOnly == 3
    Ya = Y(:,:,frames_A);
    Yb = Y(:,:,frames_B);
end

%clear Y to clear memory
clear Y

%downsample 5x 
Ya_ds = downsample_data(Ya,'time',5);
Yb_ds = downsample_data(Yb,'time',5);

%convert to uint16
Ya_ds_16 = uint16(Ya_ds);
Yb_ds_16 = uint16(Yb_ds);

%navigate to exp dir
cd(path_dir{1})

if runOnly == 1
    saveastiff(Ya_ds_16,'A_DS_nr_run.tif');
    saveastiff(Yb_ds_16,'B_DS_nr_run.tif');
elseif runOnly == 2
    saveastiff(Ya_ds_16,'A_DS_nr_norun.tif');
    saveastiff(Yb_ds_16,'B_DS_nr_norun.tif');
elseif runOnly == 3
    saveastiff(Ya_ds_16,'A_DS_nr.tif');
    saveastiff(Yb_ds_16,'B_DS_nr.tif');
end

end

%% Calculate centroid difference between A&B tuned neurons (max in field transient rate)
%TODO: put conditional here
%also generate the cent diff for either A&B tuned by either criteria

% options.tuning_criterion = 'ts';
% [cent_diff,cent_diff_AandB, pf_vector_max] = centroid_diff_single_ses(session_vars,tunedLogical, pf_vector,field_event_rates,select_fields,options);
% 
% %save the fractions output data %
% save(fullfile(path_dir{1},'cumul_analysis','centroid_diff.mat'),'cent_diff_AandB');
% %save the cent diff for all neurons (TS) tuned
% save(fullfile(path_dir{1},'cumul_analysis','centroid_diff_all.mat'),'cent_diff');
%%

    %save to directory as 16 bit tiff
    % imwrite(uint16(meanStk{ii}),'AVG_nr.tif','tif')
    
    %take max of the downsampled stack
    %     maxStk{ii} = max(Yds_16,[],3);
    %
    %     %save max projection
    %     imwrite(uint16(maxStk{ii}),'MAX_nr.tif','tif')
    %
    % meanStk{1} = mean(Y,3);

    
%%%%%% NOT USED %%%%%%%%%%%% OLD
%% Plot spatial tuning curves according to transient rate and return max bin position 
%normalized to each neuron across 

% options.tuning_criterion = 'si'; %si or ts
% %AandB, AorB, onlyA, onlyB, neither all
% %all won't work b/c no defined fields to sort by
% options.trialTuning = 'onlyA';
% options.selectSes = [1 2];
% %sort according to which trial 1 2 4 5
% options.sortTrial = 1;
% 
% %outputs:
% %placed into separate function below
% %max_transient_peak - peak idx where event rate is highest for correct A/B
% %trials
% [~,~] = plot_STC_transient_rate_single_ses(session_vars,tunedLogical,field_event_rates, pf_vector,options);
% 

end



