function [registered] = filter_matching_components(registered,tunedLogical,select_fields)

%% Define variables

%filtered ROI list across days
matching_list.original = registered.multi.assigned_filtered;


%% Select SI A and B tuned for now for each session

%generate filtered match list according to chosen tuning criterioa
[matching_list] = filter_match_tuning_crit(matching_list,tunedLogical,select_fields);


%% Parse the tuning criterion filtered matching neurons

[matching_list] = filter_match_min_event_PF(matching_list,select_fields);


%% Output

registered.multi.matching_list_filtered =  matching_list;
    
end

