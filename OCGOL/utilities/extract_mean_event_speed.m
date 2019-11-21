function [mean_speed_sel,mean_speed_sel_other] = extract_mean_event_speed(task_selective_ROIs,event_speeds,event_bins,mean_bin_speed)

%get speeds of each event in field for A/B selective ROIs
event_speeds_ex = event_speeds(task_selective_ROIs);

%get mean across all laps for selective ROI
for rr=1:size(event_speeds_ex,2)
    mean_speed_sel(rr) = mean(cell2mat(event_speeds_ex{rr}));
end

%get range each event in field for A selective
event_bins_ex = event_bins(task_selective_ROIs);

%extract the range of bins on which A events occur
for rr=1:size(event_bins_ex,2)
    %set non empty zero vectors to empty
    set_empty = cellfun(@isempty,event_bins_ex{rr},'UniformOutput',true);
    event_bins_ex{rr}(set_empty) = {[]};
    unique_bin_ex{rr} = unique(cell2mat(event_bins_ex{rr}'));
end

%mean speeds in equivalanet bins across B laps
for rr=1:size(unique_bin_ex,2)
    mean_speed_sel_other(rr) = mean(mean(mean_bin_speed(:,unique_bin_ex{rr})));
end



end

