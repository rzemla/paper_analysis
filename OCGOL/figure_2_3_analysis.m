%% Import variables and define options

%lab workstation
%input directories to matching function
%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_120618_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I42R_AB_d1_032118_1'};
path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1'};

%load place cell variables for each session
%get mat directories in each output folder
for ii=1:size(path_dir,2)
    %get matfile names for each session
    matfiles{ii} = dir([path_dir{ii},'\output','\*.mat']);
end
%load in place cell variables (and others later)
for ii = 1:size(path_dir,2)
    %add event variables
    session_vars{ii} = load(fullfile(matfiles{ii}.folder,matfiles{ii}.name),'Place_cell', 'Behavior',...
        'Behavior_split_lap','Behavior_split','Events_split','Events_split_lap', 'Imaging_split');
end

%% Define tuned logical vectors
%flag to all A or B trial or only correct A or B trials
options.allCorrect = 1; %1  = A correct; 2 = B correct
%returns struct of structs
[tunedLogical] = defineTunedLogicals(session_vars,options);

%% Generate STC maps of neurons tuned in either session and plot side by side
%customize to add options
%tuned in both sessions by SI score
%sorted by A trials

%set option as to how to select neurons for plots
options.tuning_criterion = 'si'; %si or ts
%normalized across both sessions

%PV correlation embedded inside of this function
plot_STC_OCGOL_singleSes(session_vars,tunedLogical,options)

%% Plot spatial tuning curves according to transient rate 

options.tuning_criterion = 'ts'; %si or ts
%plot_STC_transient_rate(session_vars,tunedLogical,registered,field_event_rates, pf_vector,options)


%% PV correlation matrix Aand B tuned neurons

%% Plot fraction of each neuron tuned 

tuned_si(1) = size(find(tunedLogical.si.onlyA_tuned ==1),2);
tuned_si(2) = size(find(tunedLogical.si.onlyB_tuned ==1),2);
tuned_si(3) = size(find(tunedLogical.si.AandB_tuned ==1),2);
tuned_si(4) = size(find(tunedLogical.si.neither ==1),2); 

fracTuned_si = tuned_si/sum(tuned_si);

figure
p = pie(fracTuned_si,{['A ', num2str(round(100*fracTuned_si(1))), '%'],...
                        ['B ', num2str(round(100*fracTuned_si(2))), '%'],...
                        ['A&B ', num2str(round(100*fracTuned_si(3))), '%'],...
                        ['     Neither ', num2str(round(100*fracTuned_si(4))), '%']});
hold on
title('Percentage of active neurons tuned');
colormap([0 0 1; 1 0 0; 1 0 1; 0.5 0.5 0.5])

%% Get percentage correct in each trial type and

trialOrder = session_vars{1}.Behavior.performance.trialOrder;
trialCorrect = session_vars{1}.Behavior.performance.trialCorrect;

fracA_corr = size(find(trialOrder == 2),1) / size(find(trialOrder == 2 | trialOrder == 20),1)

fracB_corr = size(find(trialOrder == 3),1) / size(find(trialOrder == 3 | trialOrder == 30),1)



