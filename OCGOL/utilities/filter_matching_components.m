function [registered] = filter_matching_components(registered,tunedLogical,select_fields,options)

%% Define variables

%filtered ROI list across days
matching_list.original = registered.multi.assigned_filtered;


%% Select SI A and B tuned for now for each session
%QC check
%generate filtered match list according to chosen tuning criteria
[matching_list] = filter_match_tuning_crit(matching_list,tunedLogical,select_fields,options);


%% Parse the tuning criterion filtered matching neurons
%makes sure that each neuron has at least 1 PF and 5 distinct events in
%field
[matching_list] = filter_match_min_event_PF(matching_list,select_fields,options);


%% Output

registered.multi.matching_list_filtered =  matching_list;
    
end

