function [outputArg1,outputArg2] = filter_matching_components(registered,tunedLogical,select_fields)

%% Define variables

%filtered ROI list across days
matching_list.original = registered.multi.assigned_filtered;


%% Select SI A and B tuned for now for each session
%copy first
matching_list.si_AB_filt = matching_list.original;

%for each session
for ss=1:size(select_fields,2)
    %get idxs
    si_idx_temp = find(tunedLogical(ss).si.AandB_tuned ==1);
    %get logical with values that are si tuned
    keep_idx_log = ismember(matching_list.si_AB_filt(:,ss),si_idx_temp);
    %set the negative of the log to nan (no match based on tuning criterion
    matching_list.si_AB_filt(~keep_idx_log,ss) = nan;

end

%% Select SI A and B fitl with at least 1 PF and 5 distinct events in field 

%convert logical field in binary
%for each session
for ss=1:size(select_fields,2)
    for tt=1:2
        %for each ROI
        for rr=1:size(select_fields{ss}{tt},2)
            %find at least 1 PF with 5 distinct events
            select_fields_binary{ss}{tt}{rr} = ~isempty(find(select_fields{ss}{tt}{rr} == 1));
        end
        %convert to matrix - each session and trial type
        select_fields_binary{ss}{tt} = cell2mat(select_fields_binary{ss}{tt});
    end
end

%do an & comparison for both trial types on each day - RESUME HERE

%%
    
end

