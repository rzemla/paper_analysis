function [matching_list] = filter_match_min_event_PF(matching_list,select_fields)


%% Generate logical vector filter to ensure at least 1 PF and 5 distinct events in field 

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

%do an & comparison for both trial types on each day for A and B tuned
%for each session (both trials)
%generate for A trials and B trials as well (naming)
for ss=1:size(select_fields,2)
    select_fieldsAB{ss} = select_fields_binary{ss}{1} & select_fields_binary{ss}{2};
    select_fieldsA{ss} = select_fields_binary{ss}{1};
    select_fieldsB{ss} = select_fields_binary{ss}{2};
end

%%
%copy match matices from previous filter
matching_list.si_AB_filt_event_filt = matching_list.si_AB_filt;

matching_list.si_Aall_filt_event_filt = matching_list.si_Aall_filt;
matching_list.si_Ball_filt_event_filt = matching_list.si_Ball_filt;

matching_list.ts_Aall_filt_event_filt = matching_list.ts_Aall_filt;
matching_list.ts_Ball_filt_event_filt = matching_list.ts_Ball_filt;

% for SI A
%for each session
for ss=1:size(select_fields,2)
    %get idxs
    event_idx_temp = find(select_fieldsA{ss} ==1);
    %get logical with values that are si tuned
    keep_idx_log = ismember(matching_list.si_Aall_filt_event_filt(:,ss),event_idx_temp);
    %set the negative of the log to nan (no match based on tuning criterion
    matching_list.si_Aall_filt_event_filt(~keep_idx_log,ss) = nan;
end

% for SI B
%for each session
for ss=1:size(select_fields,2)
    %get idxs
    event_idx_temp = find(select_fieldsB{ss} ==1);
    %get logical with values that are si tuned
    keep_idx_log = ismember(matching_list.si_Ball_filt_event_filt(:,ss),event_idx_temp);
    %set the negative of the log to nan (no match based on tuning criterion
    matching_list.si_Ball_filt_event_filt(~keep_idx_log,ss) = nan;
end


 % for SI AB
%for each session
for ss=1:size(select_fields,2)
    %get idxs
    event_idx_temp = find(select_fieldsAB{ss} ==1);
    %get logical with values that are si tuned
    keep_idx_log = ismember(matching_list.si_AB_filt_event_filt(:,ss),event_idx_temp);
    %set the negative of the log to nan (no match based on tuning criterion
    matching_list.si_AB_filt_event_filt(~keep_idx_log,ss) = nan;
end

% for TS A
%for each session
for ss=1:size(select_fields,2)
    %get idxs
    event_idx_temp = find(select_fieldsA{ss} ==1);
    %get logical with values that are si tuned
    keep_idx_log = ismember(matching_list.ts_Aall_filt_event_filt(:,ss),event_idx_temp);
    %set the negative of the log to nan (no match based on tuning criterion
    matching_list.ts_Aall_filt_event_filt(~keep_idx_log,ss) = nan;
end

% for TS B
%for each session
for ss=1:size(select_fields,2)
    %get idxs
    event_idx_temp = find(select_fieldsB{ss} ==1);
    %get logical with values that are si tuned
    keep_idx_log = ismember(matching_list.ts_Ball_filt_event_filt(:,ss),event_idx_temp);
    %set the negative of the log to nan (no match based on tuning criterion
    matching_list.ts_Ball_filt_event_filt(~keep_idx_log,ss) = nan;
end



end

