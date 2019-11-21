function [mean_speed] = event_vs_speed(session_vars, task_selective_ROIs,ROI_idx_tuning_class,...
                                    select_fields,max_transient_peak,mean_bin_speed,lap_bin_split,options)

%% Define inputs variables

selectTrial = options.selectTrial;
sessionSelect = options.sessionSelect;

%% Define/load variables for each session

%% 

%run ones for each set of trials
split_run_ones = session_vars{1}.Behavior_split_lap.run_ones;

%all run events by lap
run_onset_matrix_each_lap = session_vars{1}.Events_split_lap.Run.run_onset_binary;

%order of correct and incorrect A/B trials
trialOrder = session_vars{1}.Behavior.performance.trialOrder;

%get lap idxs associated with the correct A and B trials
corr_laps_idx.A = find(trialOrder == 2);
corr_laps_idx.B = find(trialOrder ==3);

%get speed for events
%for correct A laps
%max place field filtered

trialType = 1;
[event_speeds.A, event_bins.A] = extract_event_speeds(corr_laps_idx.A, split_run_ones, run_onset_matrix_each_lap,lap_bin_split.A,session_vars,...
                    max_transient_peak,mean_bin_speed.A,trialType);
             
%for correct B laps                
trialType = 2;
[event_speeds.B, event_bins.B] = extract_event_speeds(corr_laps_idx.B, split_run_ones, run_onset_matrix_each_lap,lap_bin_split.B,session_vars,...
                    max_transient_peak,mean_bin_speed.B,trialType);                

%Get mean speed for each event on A laps and B laps
%last input should be from opposing laps - to calculate speed within that
%bin on opposing laps
%for A selective
[mean_speed.Asel,mean_speed.Asel_B_laps] = extract_mean_event_speed(task_selective_ROIs.A.idx, event_speeds.A, event_bins.A, mean_bin_speed.B);

[mean_speed.Bsel,mean_speed.Bsel_A_laps] = extract_mean_event_speed(task_selective_ROIs.B.idx, event_speeds.B, event_bins.B, mean_bin_speed.A);


%plot scatter for A selective and B selective neurons
figure
hold on
axis square
xlim([0 25])
ylim([0 25])
xlabel('Speed A laps [cm/s]');
ylabel('Speed B laps [cm/s]');
plot([0 25],[0 25],'k--')
%A sel
scatter(mean_speed.Asel,mean_speed.Asel_B_laps,10,'filled','MarkerFaceColor', 'b')
%B sel
scatter(mean_speed.Bsel_A_laps,mean_speed.Bsel,10,'filled','MarkerFaceColor', 'r')


%% Original development code used as input to function above
%{
for ll = 1:size(corr_laps_idx.A,1)
    bin_aligned_onsets.A{ll} = run_onset_matrix_each_lap{corr_laps_idx.A(ll)}(logical(split_run_ones{corr_laps_idx.A(ll)}),:);
end

%number of ROIs in session
nbROI = size(bin_aligned_onsets.A{1, 1},2);

%overlay bin map on top of onsets for each ROI
nb_laps.A = size(corr_laps_idx.A,1);


%repmat each bin assign to number of ROIs
for ll=1:size(lap_bin_split.A,2)
    lap_bin_split_ROI_ex.A{ll} = repmat(lap_bin_split.A{ll},1,nbROI);
end

%bin assign for each ROI
for ll=1:size(lap_bin_split.A,2)
    event_run_bins_mat.A{ll} = lap_bin_split_ROI_ex.A{ll}.*bin_aligned_onsets.A{ll};
end

%for each ROI, get bin activations on each lap
for rr=1:nbROI
    %for each laps
    for ll=1:size(event_run_bins_mat.A,2)
        event_run_bin_extract{rr}{ll} = event_run_bins_mat.A{ll}((event_run_bins_mat.A{ll}(:,rr) ~=0),rr); 
    end
    
end

%place field edges for correct A laps
place_field_bin_edges.A = session_vars{1}.Place_cell{1}.placeField.edge;

%get edges associated with the max transient peak
for rr=1:nbROI
    if~isnan(max_transient_peak{1}{1}(rr))
    max_place_field_bin_edges.A{rr} = place_field_bin_edges.A{rr}(max_transient_peak{1}{1}(rr),:);
    else
        max_place_field_bin_edges.A{rr} = nan;
    end
end

%extract the events associated with max transient place field
for rr=1:nbROI
    %if place field exists
    if ~isnan(max_place_field_bin_edges.A{rr})
        
        if (max_place_field_bin_edges.A{rr}(1) < max_place_field_bin_edges.A{rr}(2))
            %define edges of place field
            bin_start = max_place_field_bin_edges.A{rr}(1);
            bin_end = max_place_field_bin_edges.A{rr}(2);
            
            for ll=1:nb_laps.A
                event_idx_temp = find(event_run_bin_extract{rr}{ll} >= bin_start &  event_run_bin_extract{rr}{ll} <= bin_end);
                %extract the event within the field on each lap
                event_run_bin_extract_pf_sel{rr}{ll} = event_run_bin_extract{rr}{ll}(event_idx_temp);
            end
        else %if split field along track end
            %split into 2 finds
            %far edge (near end of track) (from x - 100)
            bin_start = max_place_field_bin_edges.A{rr}(1);
            %near edge (toward start of track) (from 1 -x)
            bin_end = max_place_field_bin_edges.A{rr}(2);
            %iterate through each lap
            for ll=1:nb_laps.A
                end_idx = find(event_run_bin_extract{rr}{ll} >= bin_start &  event_run_bin_extract{rr}{ll} <= 100);
                start_idx = find(event_run_bin_extract{rr}{ll} >= 1 &  event_run_bin_extract{rr}{ll} <= bin_end);
                
                %combine the indices above
                event_idx_temp = sort([start_idx,end_idx]);
                %extract the event within the field on each lap
                event_run_bin_extract_pf_sel{rr}{ll} = event_run_bin_extract{rr}{ll}(event_idx_temp);
            end
        end
    else
        event_run_bin_extract_pf_sel{rr} = [];
    end
end

%extract bin speed associated with each event
for rr=1:nbROI
    if ~isempty(event_run_bin_extract_pf_sel{rr})
        for ll=1:nb_laps.A
            event_run_pf_speed{rr}{ll} = mean_bin_speed.A(ll,event_run_bin_extract_pf_sel{rr}{ll});
        end
    else
         event_run_pf_speed{rr}{ll} = [];
    end
end

%extract task-selective speeds
Asel_speed = event_run_pf_speed(task_selective_ROIs.A.idx);

%convert each ROI speed to mat
for rr=1:size(Asel_speed,2)
    Asel_speed_mat{rr} = cell2mat(Asel_speed{rr});
    
end

%get mean speed of event for each ROI
figure
hold on
plot(cellfun(@mean,Asel_speed_mat))

%%%%%%%%%%% Development code for second script above %%%%%%%%%%%%%%%%%%

%get speeds of each event in field for A/B selective ROIs
event_speeds.Asel = event_speeds.A(task_selective_ROIs.A.idx);

%get mean across all laps for selective ROI
for rr=1:size(event_speeds.Asel,2)
    mean_speed.Asel(rr) = mean(cell2mat(event_speeds.Asel{rr}));
end

%get range each event in field for A selective
event_bins.Asel = event_bins.A(task_selective_ROIs.A.idx);

%extract the range of bins on which A events occur
for rr=1:size(event_bins.Asel,2)
    %set non empty zero vectors to empty
    set_empty = cellfun(@isempty,event_bins.Asel{rr},'UniformOutput',true);
    event_bins.Asel{rr}(set_empty) = {[]};
    unique_bin.Asel{rr} = unique(cell2mat(event_bins.Asel{rr}'));
end

%mean speeds in equivalanet bins across B laps
for rr=1:size(unique_bin.Asel,2)
    mean_speed.Asel_B_laps(rr) = mean(mean(mean_bin_speed.B(:,unique_bin.Asel{rr})));
end



%}

end

