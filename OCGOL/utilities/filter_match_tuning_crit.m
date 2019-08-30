function [matching_list] = filter_match_tuning_crit(matching_list,tunedLogical,select_fields)

%copy first
matching_list.si_AB_filt = matching_list.original;
%si all A tuned; all B tuned
matching_list.si_Aall_filt = matching_list.original;
matching_list.si_Ball_filt = matching_list.original;
%si all A tuned; all B tuned
matching_list.ts_Aall_filt = matching_list.original;
matching_list.ts_Ball_filt = matching_list.original;

%for SI AB Filt
%for each session
for ss=1:size(select_fields,2)
    %get idxs
    idx_temp = find(tunedLogical(ss).si.AandB_tuned ==1);
    %get logical with values that are si tuned
    keep_idx_log = ismember(matching_list.si_AB_filt(:,ss),idx_temp);
    %set the negative of the log to nan (no match based on tuning criterion
    matching_list.si_AB_filt(~keep_idx_log,ss) = nan;

end

%for SI A filt (any A regardless of B tuning)
for ss=1:size(select_fields,2)
    %get idxs
    idx_temp = find(tunedLogical(ss).si.Atuned   == 1);
    %get logical with values that are si tuned
    keep_idx_log = ismember(matching_list.si_Aall_filt(:,ss),idx_temp);
    %set the negative of the log to nan (no match based on tuning criterion
    matching_list.si_Aall_filt(~keep_idx_log,ss) = nan;

end

%for SI B filt (any B regardless of A tuning)
for ss=1:size(select_fields,2)
    %get idxs
    idx_temp = find(tunedLogical(ss).si.Btuned   == 1);
    %get logical with values that are si tuned
    keep_idx_log = ismember(matching_list.si_Ball_filt(:,ss),idx_temp);
    %set the negative of the log to nan (no match based on tuning criterion
    matching_list.si_Ball_filt(~keep_idx_log,ss) = nan;

end

%for TS A filt (any A regardless of B tuning)
for ss=1:size(select_fields,2)
    %get idxs
    idx_temp = find(tunedLogical(ss).ts.Atuned   == 1);
    %get logical with values that are si tuned
    keep_idx_log = ismember(matching_list.ts_Aall_filt(:,ss),idx_temp);
    %set the negative of the log to nan (no match based on tuning criterion
    matching_list.ts_Aall_filt(~keep_idx_log,ss) = nan;

end

%for TS B filt (any B regardless of A tuning)
for ss=1:size(select_fields,2)
    %get idxs
    idx_temp = find(tunedLogical(ss).ts.Btuned   == 1);
    %get logical with values that are si tuned
    keep_idx_log = ismember(matching_list.ts_Ball_filt(:,ss),idx_temp);
    %set the negative of the log to nan (no match based on tuning criterion
    matching_list.ts_Ball_filt(~keep_idx_log,ss) = nan;

end


end

