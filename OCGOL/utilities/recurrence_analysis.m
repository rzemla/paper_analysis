function [recurr,frac_active,recurr_ex,frac_active_ex] = recurrence_analysis(registered,removedROI_clean,session_vars,tunedLogical,select_fields,options)

%assign select trial
selectTrial = options.selectTrial;

%set the number of sessions based on number of sessions imageed
nb_ses = size(options.sessionSelect,2);

%% Get all active neurons that are matched and post-match manually filtered

%all active assigned
all_active_match = registered.multi.assigned_filtered;

%convert to binary
all_active_match_bin = ~isnan(all_active_match);

%% Add all neurons from each session that have no matching partner to end of match matrix
%take into account soma cleaning from d1 where non-somatic point
%componenets were also included, but not in later sessions)

%day 1 cleaned components
%check if any i
d1_cleaned_ROIs = find(removedROI_clean ==1);

%check if any of the removed components have a match on d1/ses 1
if isempty(intersect(all_active_match(:,1),d1_cleaned_ROIs))
    disp('None of cleaned neurons from day 1 are in match')
else
    %remove those ROI from the removed components
    disp('Adjust removed neurons');
    nb_removed = size(d1_cleaned_ROIs,1);
    nb_overlap = size(intersect(all_active_match(:,1),d1_cleaned_ROIs),1);
    
    %remove the matching neurons (overwrite vector)
    d1_cleaned_ROIs = setdiff(d1_cleaned_ROIs,intersect(all_active_match(:,1),d1_cleaned_ROIs));
    
end

%get nb ROIs from each session
for ii=1:nb_ses
    %number ROIs
    nbROI(ii) = size(session_vars{ii}.Imaging.trace_restricted,2);
    %construct list of sequential ROIs from each ses from comp below
    ROI_list_ses{ii} = 1:nbROI(ii);
end

%get ROIs that do not have any matches to the other days
for ii=1:nb_ses
    % get all neurons that have no matches
    no_match_ROIs{ii} = setdiff(ROI_list_ses{ii},all_active_match(:,ii));
    %remove cleaned up soma from D1
    if ii==1
       no_match_ROIs{ii} = setdiff(no_match_ROIs{ii},d1_cleaned_ROIs);
    end
    
end

%% Extent the match matrix to include never matching neurons
%create extension matrix
unmatched_ROI_nb = cellfun(@(x) size(x,2),no_match_ROIs,'UniformOutput',true);
%blank matrix
extension_matrix = nan(sum(unmatched_ROI_nb),nb_ses);

%define fill indices
end_fill_idx = cumsum(unmatched_ROI_nb);
start_fill_idx = [1,end_fill_idx(1:end-1)+1];

%fill the blank matrix with 
for ii=1:nb_ses
    extension_matrix(start_fill_idx(ii):end_fill_idx(ii),ii) = no_match_ROIs{ii};
end

%extend active match matrix with no match ROIS
all_active_match_ex = [all_active_match; extension_matrix];

%convert to binary matrix
all_active_match_bin_ex = ~isnan(all_active_match_ex);

%% Convert count for each place cell into in logical vector (used in combination with TS and SI vectors)
%find A,sig T.S tuned and 5 min event with defined PF
%convert select_fields to logical - include or not include
for ss=1:nb_ses
    for tt=options.selectTrial   %for A vs.B
        %replace empty cells withs with 0
        empty_idx = find(cellfun(@isempty,select_fields{ss}{tt})==1);
        %equivalenet number of 0 cells to fill empty cells with
        zero_fill_cell(1,1:size(empty_idx,2)) = {0};
        select_fields{ss}{tt}(empty_idx) = zero_fill_cell;
        
        clear sig_field_sum empty_idx zero_fill_cell
    end
    
end

%check which fields are significant
for ss=1:nb_ses
    for tt=options.selectTrial   %for A vs.B
        %get logical output for neurons that have sig place fields
        select_field_log{ss}{tt} = logical(cell2mat(cellfun(@(x)sum(x,2), select_fields{ss}{tt},'UniformOutput',false)));
    end
end
%no_match_ROI
%extension_matrix

%% Extract ROI idxs that are tuned by some criterion and also have min5 events/sig place field

for ss=1:nb_ses
    %A
    tuned_ROI_sig_idxs{ss}.ts.allA = find((select_field_log{ss}{selectTrial(1)} & tunedLogical(ss).ts.Atuned)==1);
    tuned_ROI_sig_idxs{ss}.si.allA = find((select_field_log{ss}{selectTrial(1)} & tunedLogical(ss).si.Atuned)==1);
    
    %B 
    tuned_ROI_sig_idxs{ss}.ts.allB = find((select_field_log{ss}{selectTrial(2)} & tunedLogical(ss).ts.Btuned)==1);
    tuned_ROI_sig_idxs{ss}.si.allB = find((select_field_log{ss}{selectTrial(2)} & tunedLogical(ss).si.Btuned)==1);
end

%% Construct extended match_bin for TS/SI A/B
%copy matrices before 
ex_matrix.ts.A = extension_matrix;
ex_matrix.ts.B = extension_matrix;
ex_matrix.si.A = extension_matrix;
ex_matrix.si.B = extension_matrix;

for ss=1:nb_ses
    %A ts
    [~,select_idx_on,~] = intersect(extension_matrix(:,ss),tuned_ROI_sig_idxs{ss}.ts.allA);
    [~,select_idx_off] = setdiff(extension_matrix(:,ss),tuned_ROI_sig_idxs{ss}.ts.allA);
    %set on ROIs (1)
    ex_matrix.ts.A(select_idx_on,ss) = 1;
    %set off ROIs (0)
    ex_matrix.ts.A(select_idx_off,ss) = 0;
    
    %B ts
    [~,select_idx_on,~] = intersect(extension_matrix(:,ss),tuned_ROI_sig_idxs{ss}.ts.allB);
    [~,select_idx_off] = setdiff(extension_matrix(:,ss),tuned_ROI_sig_idxs{ss}.ts.allB);
    %set on ROIs (1)
    ex_matrix.ts.B(select_idx_on,ss) = 1;
    %set off ROIs (0)
    ex_matrix.ts.B(select_idx_off,ss) = 0; 
    
    %A si
    [~,select_idx_on,~] = intersect(extension_matrix(:,ss),tuned_ROI_sig_idxs{ss}.si.allA);
    [~,select_idx_off] = setdiff(extension_matrix(:,ss),tuned_ROI_sig_idxs{ss}.si.allA);
    %set on ROIs (1)
    ex_matrix.si.A(select_idx_on,ss) = 1;
    %set off ROIs (0)
    ex_matrix.si.A(select_idx_off,ss) = 0;
    
    %B si
    [~,select_idx_on,~] = intersect(extension_matrix(:,ss),tuned_ROI_sig_idxs{ss}.si.allB);
    [~,select_idx_off] = setdiff(extension_matrix(:,ss),tuned_ROI_sig_idxs{ss}.si.allB);
    %set on ROIs (1)
    ex_matrix.si.B(select_idx_on,ss) = 1;
    %set off ROIs (0)
    ex_matrix.si.B(select_idx_off,ss) = 0;  
end

%% SI and TS matching matrices for neurons with at least 1 match to another session
%SI binary match matrices for A and B trials
%SI
match_bin.si.A = ~isnan(registered.multi.matching_list_filtered.si_Aall_filt_event_filt);
match_bin.si.B = ~isnan(registered.multi.matching_list_filtered.si_Ball_filt_event_filt);
match_bin.si.AorB = match_bin.si.A | match_bin.si.B;

%TS
match_bin.ts.A = ~isnan(registered.multi.matching_list_filtered.ts_Aall_filt_event_filt);
match_bin.ts.B = ~isnan(registered.multi.matching_list_filtered.ts_Ball_filt_event_filt);
match_bin.ts.AorB = match_bin.ts.A | match_bin.ts.B;

%% Add the extended portion of non-match neuron tuned by SI/TS criteria

%SI A/B
match_bin_ex.si.A = [match_bin.si.A; ex_matrix.si.A];
match_bin_ex.si.B = [match_bin.si.B; ex_matrix.si.B];

%TS A/B
match_bin_ex.ts.A = [match_bin.ts.A; ex_matrix.ts.A];
match_bin_ex.ts.B = [match_bin.ts.B; ex_matrix.ts.B];


%% Calculate recurrence for neurons with at least 1 match to another session

%calculate recurrence
%for given 2 sessions, find neurons that are tuned in both divided by
%neurons that are tuned in either or both sessions
% (1/1)/(1/0 or 0/1 or 1/1)

%for session against the other
for ii=1:nb_ses
    for jj=1:nb_ses
        %TS A / B
        recurr.ts.A(ii,jj) =length(find(sum(match_bin.ts.A(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin.ts.A(:,[ii jj]),2) == 1)) + length(find(sum(match_bin.ts.A(:,[ii jj]),2) == 2)));
        recurr.ts.B(ii,jj) =length(find(sum(match_bin.ts.B(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin.ts.B(:,[ii jj]),2) == 1)) + length(find(sum(match_bin.ts.B(:,[ii jj]),2) == 2)));
        %TS A or B
        recurr.ts.AorB(ii,jj) =length(find(sum(match_bin.ts.AorB(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin.ts.AorB(:,[ii jj]),2) == 1)) + length(find(sum(match_bin.ts.AorB(:,[ii jj]),2) == 2)));
        
        %SI A / B
        recurr.si.A(ii,jj) =length(find(sum(match_bin.si.A(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin.si.A(:,[ii jj]),2) == 1)) + length(find(sum(match_bin.si.A(:,[ii jj]),2) == 2)));
        recurr.si.B(ii,jj) =length(find(sum(match_bin.si.B(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin.si.B(:,[ii jj]),2) == 1)) + length(find(sum(match_bin.si.B(:,[ii jj]),2) == 2)));
        %SI A or B
        recurr.si.AorB(ii,jj) =length(find(sum(match_bin.si.AorB(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin.si.AorB(:,[ii jj]),2) == 1)) + length(find(sum(match_bin.si.AorB(:,[ii jj]),2) == 2)));
    end
end

%% Fraction of active neurons matched between sessions (for neurons with at least 1 match to another session)
for ii=1:nb_ses
    for jj=1:nb_ses
        %relative to first session
        frac_active.first_ses(ii,jj) = length(find(sum(all_active_match_bin(:,[ii jj]),2) == 2))./length(find(sum(all_active_match_bin(:,ii),2) ==1));
        %relative to second session
        frac_active.second_ses(ii,jj) = length(find(sum(all_active_match_bin(:,[ii jj]),2) == 2))./length(find(sum(all_active_match_bin(:,jj),2) ==1));       
    end
end

%% Recurrence for neurons that do have a match partner to another session
for ii=1:nb_ses
    for jj=1:nb_ses
        %TS A / B
        recurr_ex.ts.A(ii,jj) =length(find(sum(match_bin_ex.ts.A(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin_ex.ts.A(:,[ii jj]),2) == 1)) + length(find(sum(match_bin_ex.ts.A(:,[ii jj]),2) == 2)));
        recurr_ex.ts.B(ii,jj) =length(find(sum(match_bin_ex.ts.B(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin_ex.ts.B(:,[ii jj]),2) == 1)) + length(find(sum(match_bin_ex.ts.B(:,[ii jj]),2) == 2)));
        %TS A or B
        %recurr_ex.ts.AorB(ii,jj) =length(find(sum(match_bin_ex.ts.AorB(:,[ii jj]),2) == 2))/...
        %    (length(find(sum(match_bin_ex.ts.AorB(:,[ii jj]),2) == 1)) + length(find(sum(match_bin_ex.ts.AorB(:,[ii jj]),2) == 2)));
        
        %SI A / B
        recurr_ex.si.A(ii,jj) =length(find(sum(match_bin_ex.si.A(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin_ex.si.A(:,[ii jj]),2) == 1)) + length(find(sum(match_bin_ex.si.A(:,[ii jj]),2) == 2)));
        recurr_ex.si.B(ii,jj) =length(find(sum(match_bin_ex.si.B(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin_ex.si.B(:,[ii jj]),2) == 1)) + length(find(sum(match_bin_ex.si.B(:,[ii jj]),2) == 2)));
        %SI A or B
        %recurr_ex.si.AorB(ii,jj) =length(find(sum(match_bin_ex.si.AorB(:,[ii jj]),2) == 2))/...
        %    (length(find(sum(match_bin_ex.si.AorB(:,[ii jj]),2) == 1)) + length(find(sum(match_bin_ex.si.AorB(:,[ii jj]),2) == 2)));
    end
end

%% Fraction of active neurons matched between sessions (including neurons not matching to other sessions)
for ii=1:nb_ses
    for jj=1:nb_ses
        %relative to first session
        frac_active_ex.first_ses(ii,jj) = length(find(sum(all_active_match_bin_ex(:,[ii jj]),2) == 2))./length(find(sum(all_active_match_bin_ex(:,ii),2) ==1));
        %relative to second session
        frac_active_ex.second_ses(ii,jj) = length(find(sum(all_active_match_bin_ex(:,[ii jj]),2) == 2))./length(find(sum(all_active_match_bin_ex(:,jj),2) ==1));       
    end
end

%% Plot recurrence rate across time (relative to D1) for neurons with at least 1 match
figure
hold on
title('T.S. tuned match rel D1 - at least 1 match')
plot(recurr.ts.A(1,:),'b')
plot(recurr.ts.B(1,:),'r')
%plot fraction active
ylim([0 1])
plot(frac_active.first_ses(1,:),'g')

figure
hold on
title('S.I. tuned match rel D1 - at least 1 match')
plot(recurr.si.A(1,:),'b')
plot(recurr.si.B(1,:),'r')
ylim([0 1])
plot(frac_active.first_ses(1,:),'g')

%A or B
figure
hold on
title('S.I. A or B tuned match rel D1 - at least 1 match')
plot(recurr.si.AorB(1,:),'m')
ylim([0 1])
plot(frac_active.first_ses(1,:),'g')

figure
hold on
title('T.S. A or B tuned match rel D1 - at least 1 match')
plot(recurr.ts.AorB(1,:),'m')
ylim([0 1])
plot(frac_active.first_ses(1,:),'g')

%% Plot recurrence rate across time (relative to D1) for neurons with at least 1 match
figure
hold on
title('T.S. tuned match rel D1 - all regardless of match')
plot(recurr_ex.ts.A(1,:),'b')
plot(recurr_ex.ts.B(1,:),'r')
%plot fraction active
ylim([0 1])
plot(frac_active_ex.first_ses(1,:),'g')

figure
hold on
title('S.I. tuned match rel D1 - all regardless of match')
plot(recurr_ex.si.A(1,:),'b')
plot(recurr_ex.si.B(1,:),'r')
ylim([0 1])
plot(frac_active_ex.first_ses(1,:),'g')

%A or B
% figure
% hold on
% title('S.I. A or B tuned match rel D1 - all regardless of match')
% plot(recurr_ex.si.AorB(1,:),'m')
% ylim([0 1])
% plot(frac_active_ex.first_ses(1,:),'g')

% figure
% hold on
% title('T.S. A or B tuned match rel D1 - all regardless of match')
% plot(recurr-ex.ts.AorB(1,:),'m')
% ylim([0 1])
% plot(frac_active.first_ses(1,:),'g')

end

