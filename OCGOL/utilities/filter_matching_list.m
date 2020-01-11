function [match_list_run_event_filtered] = filter_matching_list(match_list_run_event_filtered,neuron_idx_less_than_5)

%% Filter based on run events
%match_list_run_event_filtered

for ss=1:size(match_list_run_event_filtered,2)
    %get the indexes that are matched
    [~,idx_list{ss},~] = intersect(match_list_run_event_filtered(:,ss),neuron_idx_less_than_5{ss},'stable');
    %set these to nan
    match_list_run_event_filtered(idx_list{ss},ss) = nan;
    
end

end

