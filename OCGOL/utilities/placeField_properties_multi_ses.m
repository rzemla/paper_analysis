function [placeField_dist,pf_count_filtered_log, pf_count_filtered] = placeField_properties_multi_ses(session_vars, tunedLogical,select_fields,task_selective_ROIs,options)

%% Define variables
sessionSelect = options.sessionSelect;
selectTrial = options.selectTrial;


%% Bin to distance conversion factor
binToCm = 1.96;

%% Select trial tuned classes of neurons

switch options.tuning_criterion
    case 'si' %spatial information
        %for each session
        for ss =sessionSelect
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
        for ss =sessionSelect
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

%bypass and get all SI or TS tuned 

for ss =sessionSelect
    tuned{ss}.A.ts = tunedLogical(ss).ts.Atuned;
    tuned{ss}.B.ts = tunedLogical(ss).ts.Btuned;
    
    tuned{ss}.A.si = tunedLogical(ss).si.Atuned;
    tuned{ss}.B.si = tunedLogical(ss).si.Btuned;
end


%% Count number of PFs based on filter for 5 lap-distinct events in each field (for each ROI) - global (no selection)

for ss =sessionSelect
    for tt =selectTrial
        for rr =1:size(select_fields{ss}{tt},2)
            %count the place fields for each roi
            pf_count_filtered{ss}(tt,rr)= sum(select_fields{ss}{tt}{rr});
        end
    end
end

%create logical where there is at least 1 field with at least 5
%lap-distinct events
for ss =sessionSelect
    pf_count_filtered_log{ss} = ~(pf_count_filtered{ss} == 0);
end

%% Parse place field properties for task-selective place cells
%input ROIs are already filtered so no need to do it here
for ss =sessionSelect
    %task selective neurons
    task_sel{ss}.A.field_count = pf_count_filtered{ss}(selectTrial(1),task_selective_ROIs{ss}.A.idx);
    task_sel{ss}.B.field_count = pf_count_filtered{ss}(selectTrial(2),task_selective_ROIs{ss}.B.idx);
    
    
    %count number of PFs cumulatively
    for ii=1:3
        if ii < 3
            task_sel{ss}.A.field_count_total(ii) = size(find(task_sel{ss}.A.field_count == ii),2);
            task_sel{ss}.B.field_count_total(ii) = size(find(task_sel{ss}.B.field_count == ii),2);
        elseif ii ==3
            task_sel{ss}.A.field_count_total(ii) = size(find(task_sel{ss}.A.field_count >= ii),2);
            task_sel{ss}.B.field_count_total(ii) = size(find(task_sel{ss}.B.field_count >= ii),2);
        end
    end
    
    
    
    %count the field width of task-selective neurons
    task_sel{ss}.A.width = session_vars{ss}.Place_cell{options.selectTrial(1)}.placeField.width(task_selective_ROIs{ss}.A.idx);
    task_sel{ss}.B.width = session_vars{ss}.Place_cell{options.selectTrial(2)}.placeField.width(task_selective_ROIs{ss}.B.idx);
    
    task_sel{ss}.A.field_select = select_fields{ss}{options.selectTrial(1)}(task_selective_ROIs{ss}.A.idx);
    task_sel{ss}.B.field_select = select_fields{ss}{options.selectTrial(2)}(task_selective_ROIs{ss}.B.idx);
    
end

for ss =sessionSelect
    %extract the field widths for task selective neurons for A and B
    %(cumulative)
    %A trials
    for pp=1:size(task_sel{ss}.A.width,2)
        task_sel{ss}.A.field_bin_width{pp} = task_sel{ss}.A.width{pp}(task_sel{ss}.A.field_select{pp});
    end
    %B trials
    for pp=1:size(task_sel{ss}.B.width,2)
        task_sel{ss}.B.field_bin_width{pp} = task_sel{ss}.B.width{pp}(task_sel{ss}.B.field_select{pp});
    end
    %convert to cm and vector from all ROIs
    %check if any fields detected (is field count empty) 
    if ~isempty(task_sel{ss}.A.field_count)
        task_sel{ss}.A.width_cm = cell2mat(task_sel{ss}.A.field_bin_width)*binToCm;
    else
        task_sel{ss}.A.width_cm = [];
    end
    
    if ~isempty(task_sel{ss}.B.field_count)
        task_sel{ss}.B.width_cm = cell2mat(task_sel{ss}.B.field_bin_width)*binToCm;
    else
        task_sel{ss}.B.width_cm = [];
    end
end

%% Parse place field properties for task-selective place cells
for ss =sessionSelect
    %remove neurons (logical idxs) with no id'd sig. place fields
    tuned{ss}.A.ts_PFadj = ~(pf_count_filtered_log{ss}(selectTrial(1),:) == 0) & tuned{ss}.A.ts;
    tuned{ss}.B.ts_PFadj = ~(pf_count_filtered_log{ss}(selectTrial(2),:) == 0) & tuned{ss}.B.ts;
    tuned{ss}.A.si_PFadj = ~(pf_count_filtered_log{ss}(selectTrial(1),:) == 0) & tuned{ss}.A.si;
    tuned{ss}.B.si_PFadj = ~(pf_count_filtered_log{ss}(selectTrial(2),:) == 0) & tuned{ss}.B.si;
    
    %ts/si tuned select min event filtered field counts
    all{ss}.A.ts.field_count = pf_count_filtered{ss}(selectTrial(1),tuned{ss}.A.ts_PFadj);
    all{ss}.B.ts.field_count = pf_count_filtered{ss}(selectTrial(2),tuned{ss}.B.ts_PFadj);
    all{ss}.A.si.field_count = pf_count_filtered{ss}(selectTrial(1),tuned{ss}.A.si_PFadj);
    all{ss}.B.si.field_count = pf_count_filtered{ss}(selectTrial(2),tuned{ss}.B.si_PFadj);
end

for ss =sessionSelect
    %count number of PFs cumulatively
    for ii=1:3
        if ii < 3
            all{ss}.A.ts.field_count_total(ii) = size(find(all{ss}.A.ts.field_count == ii),2);
            all{ss}.B.ts.field_count_total(ii) = size(find(all{ss}.B.ts.field_count == ii),2);
            all{ss}.A.si.field_count_total(ii) = size(find(all{ss}.A.si.field_count == ii),2);
            all{ss}.B.si.field_count_total(ii) = size(find(all{ss}.B.si.field_count == ii),2);
        elseif ii ==3
            all{ss}.A.ts.field_count_total(ii) = size(find(all{ss}.A.ts.field_count >= ii),2);
            all{ss}.B.ts.field_count_total(ii) = size(find(all{ss}.B.ts.field_count >= ii),2);
            all{ss}.A.si.field_count_total(ii) = size(find(all{ss}.A.si.field_count == ii),2);
            all{ss}.B.si.field_count_total(ii) = size(find(all{ss}.B.si.field_count == ii),2);
        end
    end
end

for ss =sessionSelect
    %count the field width of task-selective neurons
    %ts
    all{ss}.A.ts.width = session_vars{ss}.Place_cell{options.selectTrial(1)}.placeField.width(tuned{ss}.A.ts_PFadj);
    all{ss}.B.ts.width = session_vars{ss}.Place_cell{options.selectTrial(2)}.placeField.width(tuned{ss}.B.ts_PFadj);
    %si
    all{ss}.A.si.width = session_vars{ss}.Place_cell{options.selectTrial(1)}.placeField.width(tuned{ss}.A.si_PFadj);
    all{ss}.B.si.width = session_vars{ss}.Place_cell{options.selectTrial(2)}.placeField.width(tuned{ss}.B.si_PFadj);
    %ts
    all{ss}.A.ts.field_select = select_fields{ss}{options.selectTrial(1)}(tuned{ss}.A.ts_PFadj);
    all{ss}.B.ts.field_select = select_fields{ss}{options.selectTrial(2)}(tuned{ss}.B.ts_PFadj);
    %si
    all{ss}.A.si.field_select = select_fields{ss}{options.selectTrial(1)}(tuned{ss}.A.si_PFadj);
    all{ss}.B.si.field_select = select_fields{ss}{options.selectTrial(2)}(tuned{ss}.B.si_PFadj);
end

for ss =sessionSelect
    %extract the field widths for task selective neurons for A and B
    %(cumulative)
    %A trials
    %ts
    for pp=1:size(all{ss}.A.ts.width,2)
        all{ss}.A.ts.field_bin_width{pp} = all{ss}.A.ts.width{pp}(all{ss}.A.ts.field_select{pp});
    end
    %si
    for pp=1:size(all{ss}.A.si.width,2)
        all{ss}.A.si.field_bin_width{pp} = all{ss}.A.si.width{pp}(all{ss}.A.si.field_select{pp});
    end
    
    %B trials
    %ts
    for pp=1:size(all{ss}.B.ts.width,2)
        all{ss}.B.ts.field_bin_width{pp} =all{ss}.B.ts.width{pp}(all{ss}.B.ts.field_select{pp});
    end
    %si
    for pp=1:size(all{ss}.B.si.width,2)
        all{ss}.B.si.field_bin_width{pp} =all{ss}.B.si.width{pp}(all{ss}.B.si.field_select{pp});
    end
end

for ss =sessionSelect
    %convert to cm and vector from all ROIs
    all{ss}.A.ts.width_cm = cell2mat(all{ss}.A.ts.field_bin_width)*binToCm;
    all{ss}.B.ts.width_cm = cell2mat(all{ss}.B.ts.field_bin_width)*binToCm;
    all{ss}.A.si.width_cm = cell2mat(all{ss}.A.si.field_bin_width)*binToCm;
    all{ss}.B.si.width_cm = cell2mat(all{ss}.B.si.field_bin_width)*binToCm;
end

%% Plot PF count for SI/TS A/B neurons and histogram of PF widths

for ss=sessionSelect
    figure;
    subplot(1,2,1)
    hold on
    title('Place field counts for task-selective neurons - SI');
    xticks([1 2 3])
    xticklabels( {'1','2','3+'});
    ylim([0 400])
    bar([all{ss}.A.si.field_count_total',all{ss}.B.si.field_count_total'])
    subplot(1,2,2)
    hold on
    title('Place field counts for task-selective neurons - TS');
    xticks([1 2 3])
    xticklabels( {'1','2','3+'});
    ylim([0 400])
    bar([all{ss}.A.ts.field_count_total',all{ss}.B.ts.field_count_total'])
end

for ss=sessionSelect
    %plot histograms
    figure;
    subplot(2,2,1)
    hold on
    histogram(all{ss}.A.ts.width_cm ,0:5:70,'Normalization','probability')
    title('A tuned- TS')
    xlabel('Width [cm]')
    ylabel('Normalized density')
    xlim([15 70]);
    ylim([0 0.4]);
    
    subplot(2,2,2)
    hold on
    title('B tuned- TS')
    histogram(all{ss}.B.ts.width_cm ,0:5:70,'Normalization','probability')
    xlabel('Width [cm]')
    ylabel('Normalized density')
    xlim([15 70]);
    ylim([0 0.4]);
    
    subplot(2,2,3)
    hold on
    histogram(all{ss}.A.si.width_cm ,0:5:70,'Normalization','probability')
    title('A tuned- SI')
    xlabel('Width [cm]')
    ylabel('Normalized density')
    xlim([15 70]);
    ylim([0 0.4]);
    
    subplot(2,2,4)
    hold on
    title('B tuned- SI')
    histogram(all{ss}.B.si.width_cm ,0:5:70,'Normalization','probability')
    xlabel('Width [cm]')
    ylabel('Normalized density')
    xlim([15 70]);
    ylim([0 0.4]);
end

%% Export stuct

%place fields for all A and B tuned neurons by either TS or SI with place
%field and event filters in place
placeField_dist.all = all;
placeField_dist.task_sel = task_sel;

end
