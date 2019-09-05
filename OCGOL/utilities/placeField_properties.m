function [placeField_dist,pf_count_filtered_log, pf_count_filtered] = placeField_properties(session_vars, tunedLogical,select_fields,task_selective_ROIs,options)

%% Bin to distance conversion factor
binToCm = 2;

%% Select trial tuned classes of neurons

switch options.tuning_criterion
    case 'si' %spatial information
        %for each session
        for ss =1:size(session_vars,2)
            %spatial information criterion
            Atuned{ss} = tunedLogical(ss).si.Atuned;
            Btuned{ss} = tunedLogical(ss).si.Btuned;
            
            AandB_tuned{ss} =  tunedLogical(ss).si.AandB_tuned;
            AorB_tuned{ss} = tunedLogical(ss).si.AorB_tuned;
            onlyA_tuned{ss} = tunedLogical(ss).si.onlyA_tuned;
            onlyB_tuned{ss} = tunedLogical(ss).si.onlyB_tuned;
            AxorB_tuned{ss} =  AorB_tuned{ss} & ~AandB_tuned{ss};
            %all tuned neurons - logical
            all_neurons{ss} = true(size(Atuned{ss}));
            %tuned to neither environment - logical
            neither_tuned{ss} = ~((onlyA_tuned{ss} | onlyB_tuned{ss}) | AandB_tuned{ss});
            
        end
    case 'ts' %spatial information
        for ss =1:size(session_vars,2)
            %spatial information criterion
            Atuned{ss} = tunedLogical(ss).ts.Atuned;
            Btuned{ss} = tunedLogical(ss).ts.Btuned;
            
            AandB_tuned{ss} =  tunedLogical(ss).ts.AandB_tuned;
            AorB_tuned{ss} = tunedLogical(ss).ts.AorB_tuned;
            onlyA_tuned{ss} = tunedLogical(ss).ts.onlyA_tuned;
            onlyB_tuned{ss} = tunedLogical(ss).ts.onlyB_tuned;
            AxorB_tuned{ss} =  AorB_tuned{ss} & ~AandB_tuned{ss};
            
            all_neurons{ss} = true(size(Atuned{ss}));
            %tuned to neither environment - logical
            neither_tuned{ss} = ~((onlyA_tuned{ss} | onlyB_tuned{ss}) | AandB_tuned{ss});
        end
end

%bypass and get all SI or TS tuned (first session - adopt later to
%multisession)
tuned{1}.A.ts = tunedLogical(ss).ts.Atuned;
tuned{1}.B.ts = tunedLogical(ss).ts.Btuned;

tuned{1}.A.si = tunedLogical(ss).si.Atuned;
tuned{1}.B.si = tunedLogical(ss).si.Btuned;

%% Count number of PFs based on filter for 5 lap-distinct events in each field (for each ROI) - global (no selection)
for tt =1:size(options.selectTrial,2)
    for rr =1:size(select_fields{1}{tt},2)
        %count the place fields for each roi
        pf_count_filtered(tt,rr)= sum(select_fields{1}{tt}{rr});
    end
end

%create logical where there is at least 1 field with at least 5
%lap-distinct events
pf_count_filtered_log = ~(pf_count_filtered == 0);

%% Parse place field properties for task-selective place cells
%input ROIs are already filtered so no need to do it here

%task selective neurons
task_sel.A.field_count = pf_count_filtered(1,task_selective_ROIs.A.idx);
task_sel.B.field_count = pf_count_filtered(2,task_selective_ROIs.B.idx);

%count number of PFs cumulatively
for ii=1:3
    if ii < 3
        task_sel.A.field_count_total(ii) = size(find(task_sel.A.field_count == ii),2);
        task_sel.B.field_count_total(ii) = size(find(task_sel.B.field_count == ii),2);
    elseif ii ==3
        task_sel.A.field_count_total(ii) = size(find(task_sel.A.field_count >= ii),2);
        task_sel.B.field_count_total(ii) = size(find(task_sel.B.field_count >= ii),2);
    end
end

%count the field width of task-selective neurons
task_sel.A.width = session_vars{1}.Place_cell{options.selectTrial(1)}.placeField.width(task_selective_ROIs.A.idx);
task_sel.B.width = session_vars{1}.Place_cell{options.selectTrial(2)}.placeField.width(task_selective_ROIs.B.idx);

task_sel.A.field_select = select_fields{1}{options.selectTrial(1)}(task_selective_ROIs.A.idx);
task_sel.B.field_select = select_fields{1}{options.selectTrial(2)}(task_selective_ROIs.B.idx);

%extract the field widths for task selective neurons for A and B
%(cumulative)
%A trials
for pp=1:size(task_sel.A.width,2) 
    task_sel.A.field_bin_width{pp} = task_sel.A.width{pp}(task_sel.A.field_select{pp});
end
%B trials
for pp=1:size(task_sel.B.width,2) 
    task_sel.B.field_bin_width{pp} = task_sel.B.width{pp}(task_sel.B.field_select{pp});
end
%convert to cm and vector from all ROIs
task_sel.A.width_cm = cell2mat(task_sel.A.field_bin_width)*binToCm;
task_sel.B.width_cm = cell2mat(task_sel.B.field_bin_width)*binToCm;

%% Parse place field properties for task-selective place cells

%remove neurons (logical idxs) with no id'd sig. place fields
tuned{1}.A.ts_PFadj = ~(pf_count_filtered_log(1,:) == 0) & tuned{1}.A.ts;
tuned{1}.B.ts_PFadj = ~(pf_count_filtered_log(2,:) == 0) & tuned{1}.B.ts;
tuned{1}.A.si_PFadj = ~(pf_count_filtered_log(1,:) == 0) & tuned{1}.A.si;
tuned{1}.B.si_PFadj = ~(pf_count_filtered_log(2,:) == 0) & tuned{1}.B.si;

%ts/si tuned select min event filtered field counts
all.A.ts.field_count = pf_count_filtered(1,tuned{1}.A.ts_PFadj);
all.B.ts.field_count = pf_count_filtered(2,tuned{1}.B.ts_PFadj);
all.A.si.field_count = pf_count_filtered(1,tuned{1}.A.si_PFadj);
all.B.si.field_count = pf_count_filtered(2,tuned{1}.B.si_PFadj);

%count number of PFs cumulatively
for ii=1:3
    if ii < 3
        all.A.ts.field_count_total(ii) = size(find(all.A.ts.field_count == ii),2);
        all.B.ts.field_count_total(ii) = size(find(all.B.ts.field_count == ii),2);
        all.A.si.field_count_total(ii) = size(find(all.A.si.field_count == ii),2);
        all.B.si.field_count_total(ii) = size(find(all.B.si.field_count == ii),2);        
    elseif ii ==3
        all.A.ts.field_count_total(ii) = size(find(all.A.ts.field_count >= ii),2);
        all.B.ts.field_count_total(ii) = size(find(all.B.ts.field_count >= ii),2);
        all.A.si.field_count_total(ii) = size(find(all.A.si.field_count == ii),2);
        all.B.si.field_count_total(ii) = size(find(all.B.si.field_count == ii),2);         
    end
end

%count the field width of task-selective neurons
%ts
all.A.ts.width = session_vars{1}.Place_cell{options.selectTrial(1)}.placeField.width(tuned{1}.A.ts_PFadj);
all.B.ts.width = session_vars{1}.Place_cell{options.selectTrial(2)}.placeField.width(tuned{1}.B.ts_PFadj);
%si
all.A.si.width = session_vars{1}.Place_cell{options.selectTrial(1)}.placeField.width(tuned{1}.A.si_PFadj);
all.B.si.width = session_vars{1}.Place_cell{options.selectTrial(2)}.placeField.width(tuned{1}.B.si_PFadj);
%ts
all.A.ts.field_select = select_fields{1}{options.selectTrial(1)}(tuned{1}.A.ts_PFadj);
all.B.ts.field_select = select_fields{1}{options.selectTrial(2)}(tuned{1}.B.ts_PFadj);
%si
all.A.si.field_select = select_fields{1}{options.selectTrial(1)}(tuned{1}.A.si_PFadj);
all.B.si.field_select = select_fields{1}{options.selectTrial(2)}(tuned{1}.B.si_PFadj);

%extract the field widths for task selective neurons for A and B
%(cumulative)
%A trials
%ts
for pp=1:size(all.A.ts.width,2) 
    all.A.ts.field_bin_width{pp} = all.A.ts.width{pp}(all.A.ts.field_select{pp});
end
%si
for pp=1:size(all.A.si.width,2) 
    all.A.si.field_bin_width{pp} = all.A.si.width{pp}(all.A.si.field_select{pp});
end

%B trials
%ts
for pp=1:size(all.B.ts.width,2) 
    all.B.ts.field_bin_width{pp} =all.B.ts.width{pp}(all.B.ts.field_select{pp});
end
%si
for pp=1:size(all.B.si.width,2) 
    all.B.si.field_bin_width{pp} =all.B.si.width{pp}(all.B.si.field_select{pp});
end

%convert to cm and vector from all ROIs
all.A.ts.width_cm = cell2mat(all.A.ts.field_bin_width)*binToCm;
all.B.ts.width_cm = cell2mat(all.B.ts.field_bin_width)*binToCm;
all.A.si.width_cm = cell2mat(all.A.si.field_bin_width)*binToCm;
all.B.si.width_cm = cell2mat(all.B.si.field_bin_width)*binToCm;

%% Plot PF count for SI/TS A/B neurons and histogram of PF widths

figure;
subplot(1,2,1)
hold on
title('Place field counts for task-selective neurons - SI');
xticks([1 2 3])
xticklabels( {'1','2','3+'});
bar([all.A.si.field_count_total',all.B.si.field_count_total'])
subplot(1,2,2)
hold on
title('Place field counts for task-selective neurons - TS');
xticks([1 2 3])
xticklabels( {'1','2','3+'});
bar([all.A.ts.field_count_total',all.B.ts.field_count_total'])

%plot histograms
figure;
subplot(2,2,1)
hold on
histogram(all.A.ts.width_cm ,0:5:70,'Normalization','probability')
title('A tuned- TS')
xlabel('Width [cm]')
ylabel('Normalized density')
xlim([15 70]);
ylim([0 0.4]);

subplot(2,2,2)
hold on
title('B tuned- TS')
histogram(all.B.ts.width_cm ,0:5:70,'Normalization','probability')
xlabel('Width [cm]')
ylabel('Normalized density')
xlim([15 70]);
ylim([0 0.4]);

subplot(2,2,3)
hold on
histogram(all.A.si.width_cm ,0:5:70,'Normalization','probability')
title('A tuned- SI')
xlabel('Width [cm]')
ylabel('Normalized density')
xlim([15 70]);
ylim([0 0.4]);

subplot(2,2,4)
hold on
title('B tuned- SI')
histogram(all.B.si.width_cm ,0:5:70,'Normalization','probability')
xlabel('Width [cm]')
ylabel('Normalized density')
xlim([15 70]);
ylim([0 0.4]);

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
%% Export stuct
% 

%place fields for all A and B tuned neurons by either TS or SI with place
%field and event filters in place
placeField_dist.all = all;
placeField_dist.task_sel = task_sel;

end
