function [outputArg1,outputArg2] = si_score_comparison(session_vars, registered)

%combine S.I all with T.S all into one matrix
% a_all_log = ~isnan(registered.multi.matching_list_filtered.si_Aall_filt);
% b_all_log = ~isnan(registered.multi.matching_list_filtered.si_Ball_filt);

% comb_ab = nan(353,7)
% comb_ab(a_all_log) = registered.multi.matching_list_filtered.si_Aall_filt(a_all_log)
% comb_ab(b_all_log) = registered.multi.matching_list_filtered.si_Ball_filt(b_all_log)


%override matching choices
comb_ab = registered.multi.matching_list_filtered.si_Aall_filt_event_filt
ses_comp = [1,2];

%day2day matches
match_day = find(sum(~isnan(comb_ab(:,ses_comp)),2) ==2);
match_ROIs_only = comb_ab(match_day,ses_comp);

%get A si scores
%first session
scores_match_all{ses_comp(1),ses_comp(2)}.A(:,1) = session_vars{ses_comp(1)}.Place_cell{4}.Spatial_Info.Spatial_Info(8,match_ROIs_only(:,1));
%second session
scores_match_all{ses_comp(1),ses_comp(2)}.A(:,2) = session_vars{ses_comp(2)}.Place_cell{4}.Spatial_Info.Spatial_Info(8,match_ROIs_only(:,2));

scores_match_all{ses_comp(1),ses_comp(2)}.B(:,1) = session_vars{ses_comp(1)}.Place_cell{5}.Spatial_Info.Spatial_Info(8,match_ROIs_only(:,1));
%second session
scores_match_all{ses_comp(1),ses_comp(2)}.B(:,2) = session_vars{ses_comp(2)}.Place_cell{5}.Spatial_Info.Spatial_Info(8,match_ROIs_only(:,2));

%% Plot
figure
subplot(1,2,1)
hold on
title('A')
axis square
xlabel('Ses 1')
ylabel('Ses 2')
xlim([ 0 0.2])
ylim([0 0.2])
scatter(scores_match_all{ses_comp(1),ses_comp(2)}.A(:,1),scores_match_all{ses_comp(1),ses_comp(2)}.A(:,2))
plot([0,0.2],[0 0.2],'k--')

subplot(1,2,2)
hold on
title('B')
axis square
xlabel('Ses 1')
ylabel('Ses 2')
xlim([ 0 0.2])
ylim([0 0.2])
scatter(scores_match_all{ses_comp(1),ses_comp(2)}.B(:,1),scores_match_all{ses_comp(1),ses_comp(2)}.B(:,2))
plot([0,0.2],[0 0.2],'k--')

%%

figure
subplot(1,2,1)
hold on
%title('B')
axis square
xlabel('Ses 1 - A')
ylabel('Ses 1 - B')
xlim([ 0 0.2])
ylim([0 0.2])
scatter(scores_match_all{ses_comp(1),ses_comp(2)}.A(:,1),scores_match_all{ses_comp(1),ses_comp(2)}.B(:,1))
plot([0,0.2],[0 0.2],'k--')

subplot(1,2,2)
hold on
%title('B')
axis square
xlabel('Ses 2 - A')
ylabel('Ses 2 - B')
xlim([ 0 0.2])
ylim([0 0.2])
scatter(scores_match_all{ses_comp(1),ses_comp(2)}.A(:,2),scores_match_all{ses_comp(1),ses_comp(2)}.B(:,2))
plot([0,0.2],[0 0.2],'k--')

end

