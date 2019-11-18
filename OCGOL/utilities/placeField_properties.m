function [placeField_dist,pf_count_filtered_log, pf_count_filtered] = placeField_properties(session_vars,tunedLogical,...
                                            select_fields,task_selective_ROIs,ROI_idx_tuning_class,options)

%% Bin to distance conversion factor
binToCm = 1.96;

%% Select trial tuned classes of neurons (use logicals as input)
%prefiltered for min 5 events on distinct laps

%S.I.
%for each session
for ss =1:size(session_vars,2)
    %spatial information criterion - regardless if tuned in other session
    Atuned.si{ss} = ROI_idx_tuning_class.si.log.Aonly | ROI_idx_tuning_class.si.log.AB;
    Btuned.si{ss} = ROI_idx_tuning_class.si.log.Bonly | ROI_idx_tuning_class.si.log.AB;

    onlyA_tuned.si{ss} = ROI_idx_tuning_class.si.log.Aonly;
    onlyB_tuned.si{ss} = ROI_idx_tuning_class.si.log.Bonly;    
    AandB_tuned.si{ss} = ROI_idx_tuning_class.si.log.AB;
    neither_tuned.si{ss} = ROI_idx_tuning_class.si.log.N;
    
end

%T.S.
for ss =1:size(session_vars,2)
    %regardless if tuned in other session
    Atuned.ts{ss} = ROI_idx_tuning_class.ts.log.Aonly | ROI_idx_tuning_class.ts.log.AB;
    Btuned.ts{ss} = ROI_idx_tuning_class.ts.log.Bonly | ROI_idx_tuning_class.ts.log.AB;

    onlyA_tuned.ts{ss} = ROI_idx_tuning_class.ts.log.Aonly;
    onlyB_tuned.ts{ss} = ROI_idx_tuning_class.ts.log.Bonly;    
    AandB_tuned.ts{ss} = ROI_idx_tuning_class.ts.log.AB;
    neither_tuned.ts{ss} = ROI_idx_tuning_class.ts.log.N;
end


for ss =1:size(session_vars,2)
    %blank all neuron idx logical vector
    all_neurons{ss} = true(size(Atuned.si{ss}));
    %blank logical equal to size of number of neurons
    blank_log = false(1,size(all_neurons{ss},2));
end

%A/B selective logicals; all A, all B, A&B by either criterion
for ss=1:size(session_vars,2)
    
    %A selective and B selective neurons
    A_sel{ss} = blank_log;
    A_sel{ss}(task_selective_ROIs.A.idx) = 1;
    
    B_sel{ss} = blank_log;
    B_sel{ss}(task_selective_ROIs.B.idx) = 1;
    
    %either criterion - A all,B all, A&B all
    Atuned.si_ts{ss} = Atuned.si{ss} | Atuned.ts{ss};
    Btuned.si_ts{ss} = Btuned.si{ss} | Btuned.ts{ss};
    AandB_tuned.si_ts{ss} = Atuned.si_ts{ss} & Btuned.si_ts{ss}; 
end

%% Convert to indices from logicals for input

%selective
A_sel_idx = find(A_sel{1} ==1);
B_sel_idx = find(B_sel{1} ==1);

%all A or B by either criterion or both
%S.I.
Atuned_idx.si = find(Atuned.si{1} ==1);
Btuned_idx.si = find(Btuned.si{1} ==1);
%T.S.
Atuned_idx.ts = find(Atuned.ts{1} ==1);
Btuned_idx.ts = find(Btuned.ts{1} ==1);
%S.I. or T.S
Atuned_idx.si_ts = find(Atuned.si_ts{1} ==1);
Btuned_idx.si_ts = find(Btuned.si_ts{1} ==1);

%AandB - either criterion or both
%S.I.
AandB_tuned_idx.si = find(AandB_tuned.si{1} ==1);
%T.S.
AandB_tuned_idx.ts = find(AandB_tuned.ts{1} ==1);
%S.I. or T.S.
AandB_tuned_idx.si_ts = find(AandB_tuned.si_ts{1} ==1);

%only A or B by either criterion
%S.I.
onlyA_tuned_idx.si = find(onlyA_tuned.si{1} ==1);
onlyB_tuned_idx.si = find(onlyB_tuned.si{1} ==1);
%T.S.
onlyA_tuned_idx.ts = find(onlyA_tuned.ts{1} ==1);
onlyB_tuned_idx.ts = find(onlyB_tuned.ts{1} ==1);


%% Count number of detected PFs based on filter for 5 lap-distinct events in each field (for each ROI) - global (no selection)

%for each trial type
for tt =1:size(options.selectTrial,2)
    for rr =1:size(select_fields{1}{tt},2)
        %count the place fields for each roi
        %row is trial type
        %column is each ROI
        pf_count_filtered(tt,rr)= sum(select_fields{1}{tt}{rr});
    end
end

%create logical where there is at least 1 field with at least 5
%lap-distinct events
pf_count_filtered_log = ~(pf_count_filtered == 0);

%% Get place field width and number of fields for task-selective place cells
%input ROIs are already filtered so no need to do it here

%count the field bin width of task-selective neurons (includes both sig and
%in-sig field - must be filtered below)
task_sel.A.bin_width_not_filt = session_vars{1}.Place_cell{options.selectTrial(1)}.placeField.width(task_selective_ROIs.A.idx);
task_sel.B.bin_width_not_filt = session_vars{1}.Place_cell{options.selectTrial(2)}.placeField.width(task_selective_ROIs.B.idx);

%task selective neurons (only fields that are significant) - filtered
task_sel.A.field_count = pf_count_filtered(1,task_selective_ROIs.A.idx);
task_sel.B.field_count = pf_count_filtered(2,task_selective_ROIs.B.idx);

%select fields that match the place field crtieria (logicals corresponding
%to sign field 
task_sel.A.field_select = select_fields{1}{options.selectTrial(1)}(task_selective_ROIs.A.idx);
task_sel.B.field_select = select_fields{1}{options.selectTrial(2)}(task_selective_ROIs.B.idx);

%count number of PFs cumulatively (inputs are filtered)
%number of ROIs with 1,2,3+ place fields
for ii=1:3
    if ii < 3
        task_sel.A.field_count_total(ii) = size(find(task_sel.A.field_count == ii),2);
        task_sel.B.field_count_total(ii) = size(find(task_sel.B.field_count == ii),2);
        
    elseif ii ==3 %3 or more fields
        task_sel.A.field_count_total(ii) = size(find(task_sel.A.field_count >= ii),2);
        task_sel.B.field_count_total(ii) = size(find(task_sel.B.field_count >= ii),2);
    end
end

%extract the field widths for task selective neurons for A and B

%A selective
%for each ROI in class
for rr=1:size(task_sel.A.bin_width_not_filt,2)
    %extract the relevant fields
    task_sel.A.field_bin_width_filt{rr} = task_sel.A.bin_width_not_filt{rr}(task_sel.A.field_select{rr});
end

%B trials
%for each ROI in class
for rr=1:size(task_sel.B.bin_width_not_filt,2) 
    task_sel.B.field_bin_width_filt{rr} = task_sel.B.bin_width_not_filt{rr}(task_sel.B.field_select{rr});
end

%convert bin width to cm from selective ROIs
%combine widths into single vector and multiply by conversion factor
task_sel.A.width_cm = cell2mat(task_sel.A.field_bin_width_filt)*binToCm;
task_sel.B.width_cm = cell2mat(task_sel.B.field_bin_width_filt)*binToCm;

%QC check with histogram
figure
subplot(1,2,1)
hold on
ylim([0 0.7])
histogram(task_sel.A.width_cm, 10,'Normalization','probability')

subplot(1,2,2)
hold on
ylim([0 0.7])
histogram(task_sel.B.width_cm, 10,'Normalization','probability')

mean(task_sel.A.width_cm)
mean(task_sel.B.width_cm)

%% Get place field width and number of fields for other classes of place cells (code above made universal here)
%turn this info function based on above
%input - idx, trial type (1 2 4 5) , session_vars, select_fields, pf_count_filtered

%A only
%SI
[other_classes.si.Aonly] = extract_pf_width_count(onlyA_tuned_idx.si, 1, session_vars, select_fields, pf_count_filtered,binToCm);
%TS
[other_classes.ts.Aonly] = extract_pf_width_count(onlyA_tuned_idx.ts, 1, session_vars, select_fields, pf_count_filtered,binToCm);

% A all
%SI
[other_classes.si.Atuned] = extract_pf_width_count(Atuned_idx.si, 1, session_vars, select_fields, pf_count_filtered,binToCm);
%TS
[other_classes.ts.Atuned] = extract_pf_width_count(Atuned_idx.ts, 1, session_vars, select_fields, pf_count_filtered,binToCm);
%TS or SI
[other_classes.si_ts.Atuned] = extract_pf_width_count(Atuned_idx.si_ts, 1, session_vars, select_fields, pf_count_filtered,binToCm);


%B only
%SI
[other_classes.si.Bonly] = extract_pf_width_count(onlyB_tuned_idx.si, 2, session_vars, select_fields, pf_count_filtered,binToCm);
%TS
[other_classes.ts.Bonly] = extract_pf_width_count(onlyB_tuned_idx.ts, 2, session_vars, select_fields, pf_count_filtered,binToCm);

% B all
%SI
[other_classes.si.Btuned] = extract_pf_width_count(Btuned_idx.si, 2, session_vars, select_fields, pf_count_filtered,binToCm);
%TS
[other_classes.ts.Btuned] = extract_pf_width_count(Btuned_idx.ts, 2, session_vars, select_fields, pf_count_filtered,binToCm);
%TS or SI
[other_classes.si_ts.Btuned] = extract_pf_width_count(Btuned_idx.si_ts, 2, session_vars, select_fields, pf_count_filtered,binToCm);

%A&B
%SI
%A laps
[other_classes.si.AB.A] = extract_pf_width_count(AandB_tuned_idx.si, 1, session_vars, select_fields, pf_count_filtered,binToCm);
%B laps
[other_classes.si.AB.B] = extract_pf_width_count(AandB_tuned_idx.si, 2, session_vars, select_fields, pf_count_filtered,binToCm);
%TS
%A laps
[other_classes.ts.AB.A] = extract_pf_width_count(AandB_tuned_idx.ts, 1, session_vars, select_fields, pf_count_filtered,binToCm);
%B laps
[other_classes.ts.AB.B] = extract_pf_width_count(AandB_tuned_idx.ts, 2, session_vars, select_fields, pf_count_filtered,binToCm);
%SI or TS
%A laps
[other_classes.si_ts.AB.A] = extract_pf_width_count(AandB_tuned_idx.si_ts, 1, session_vars, select_fields, pf_count_filtered,binToCm);
%B laps
[other_classes.si_ts.AB.B] = extract_pf_width_count(AandB_tuned_idx.si_ts, 2, session_vars, select_fields, pf_count_filtered,binToCm);


%% Plot PF count for SI/TS A/B neurons and histogram of PF widths

%place field counts
figure;
hold on
title('PF counts - task-selective A and B');
xticks([1 2 3])
xticklabels( {'1','2','3+'});
bar([task_sel.A.field_count_total',task_sel.B.field_count_total'])

%plot histograms of pf width
figure;
subplot(1,2,1)
hold on
histogram(task_sel.A.width_cm ,0:5:70,'Normalization','probability')
title('Task selective A')
xlabel('Width [cm]')
ylabel('Normalized density')
xlim([15 70]);
ylim([0 0.7]);

subplot(1,2,2)
hold on
title('Task selective B')
histogram(task_sel.B.width_cm ,0:5:70,'Normalization','probability')
xlabel('Width [cm]')
ylabel('Normalized density')
xlim([15 70]);
ylim([0 0.7]);

%% Export structure

%place fields for all A and B tuned neurons by either TS or SI with place
%field and event filters in place
placeField_dist.other_classes = other_classes;
placeField_dist.task_sel = task_sel;


%% OLD CODE BELOW - preserve for for sake of compatibility
%RESUME WITH EXPORT OF VARIABLES TO BE USED IN REMAPPING FILTER
%{
%% Count neurons with 1,2,3+ number of place fields from filtered and tuning criteria selected groups

%count # of neurons with 1 2 or 3+ place fields
for ii=1:3
    if ii < 3
        field_count_A(ii) = size(find(nb_fields_A == ii),2);
        field_count_B(ii) = size(find(nb_fields_B == ii),2);
    elseif ii ==3
        field_count_A(ii) = size(find(nb_fields_A >= ii),2);
        field_count_B(ii) = size(find(nb_fields_B >= ii),2);
    end
end


%% Get counts of place field numbers for each sub-class

%place centers of each neuron
centers_Aonly = session_vars{1}.Place_cell{1}.placeField.center(onlyA_tuned{1});
centers_Bonly = session_vars{1}.Place_cell{2}.placeField.center(onlyB_tuned{1});

%% Count the fields

%get # of place fields for each neuron
nb_fields_A = cellfun('size', centers_Aonly,1);
nb_fields_B = cellfun('size', centers_Bonly,1);

%count # of neurons with 1 2 or 3+ place fields
for ii=1:3
    if ii < 3
        field_count_A(ii) = size(find(nb_fields_A == ii),2);
        field_count_B(ii) = size(find(nb_fields_B == ii),2);
    elseif ii ==3
        field_count_A(ii) = size(find(nb_fields_A >= ii),2);
        field_count_B(ii) = size(find(nb_fields_B >= ii),2);
    end
end
}%

%% Bar plot
figure;
hold on;
title('Place field counts');
xticks([1 2 3])
xticklabels( {'1','2','3+'});
bar([field_count_A',field_count_B'])

%% Get width distributions here

%100 bins -> each bin ~2 cm - get exact number in future
width_Aonly = session_vars{1}.Place_cell{1}.placeField.width(onlyA_tuned{1});
width_Bonly = session_vars{1}.Place_cell{2}.placeField.width(onlyB_tuned{1});

%get cm widths cumulative and plot histogram
width_cm_Aonly = cell2mat(width_Aonly)*binToCm;
width_cm_Bonly = cell2mat(width_Bonly)*binToCm;

%plot histograms
figure;
subplot(1,2,1)
hold on
histogram(width_cm_Aonly,0:5:70,'Normalization','probability')
xlabel('Width [cm]')
ylabel('Normalized density')
xlim([15 70]);
ylim([0 0.4]);

subplot(1,2,2)
hold on
histogram(width_cm_Bonly,0:5:70,'Normalization','probability')
xlabel('Width [cm]')
ylabel('Normalized density')
xlim([15 70]);
ylim([0 0.4]);

% placeField_dist.width_cm_Aonly = width_cm_Aonly;
% placeField_dist.width_cm_Bonly = width_cm_Bonly;
% placeField_dist.field_count_A = field_count_A;
% placeField_dist.field_count_B = field_count_B;

%}
end
