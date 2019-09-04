function [placeField_dist] = placeField_properties(session_vars, tunedLogical,select_fields,task_selective_ROIs,options)


%bin to distance conversion factor
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


%% Get counts of place field numbers for each sub-class

%place centers of each neuron
centers_Aonly = session_vars{1}.Place_cell{1}.placeField.center(onlyA_tuned{1});
centers_Bonly = session_vars{1}.Place_cell{2}.placeField.center(onlyB_tuned{1});

%% Count number of PFs based on filter for 5 lap-distinct events in each field (for each ROI) - global (no selection)
for tt =1:size(options.selectTrial,2)
    for rr =1:size(select_fields{1}{tt},2)
        %count the place fields for each roi
        pf_count_filtered(tt,rr)= sum(select_fields{1}{tt}{rr});
    end
end

%% Select which categories of neurons to plot/compare

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

%% All A / B tuned neurons by either ST or TS score



%% Plot PF count for A/B selective and histogram of PF widths

figure;
hold on
title('Place field counts for task-selective neurons');
xticks([1 2 3])
xticklabels( {'1','2','3+'});
bar([task_sel.A.field_count_total',task_sel.B.field_count_total'])

%plot histograms
figure;
subplot(1,2,1)
hold on
histogram(task_sel.A.width_cm  ,0:5:70,'Normalization','probability')
title('A selective')
xlabel('Width [cm]')
ylabel('Normalized density')
xlim([15 70]);
ylim([0 0.4]);

subplot(1,2,2)
hold on
title('B selective')
histogram(task_sel.B.width_cm  ,0:5:70,'Normalization','probability')
xlabel('Width [cm]')
ylabel('Normalized density')
xlim([15 70]);
ylim([0 0.4]);

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

%% Export stuct

placeField_dist.width_cm_Aonly = width_cm_Aonly;
placeField_dist.width_cm_Bonly = width_cm_Bonly;
placeField_dist.field_count_A = field_count_A;
placeField_dist.field_count_B = field_count_B;

end
