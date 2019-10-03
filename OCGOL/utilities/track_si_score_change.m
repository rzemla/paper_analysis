function [outputArg1,outputArg2] = track_si_score_change(session_vars,registered,options)


%% Parameters 

sessionSelect = options.sessionSelect;
selectTrial = options.selectTrial;

%match list - all - take in only first 6 sessions/no punish
matchROI_all = registered.multi.assigned_filtered(:,1:6);

%ROIs for A and B trials (TS)
matchROI_ts_event.A = registered.multi.matching_list_filtered.ts_Aall_filt_event_filt(:,1:6);
matchROI_ts_event.B = registered.multi.matching_list_filtered.ts_Ball_filt_event_filt(:,1:6);

%ROIs for A and B trials (SI)
matchROI_si_event.A = registered.multi.matching_list_filtered.si_Aall_filt_event_filt(:,1:6);
matchROI_si_event.B = registered.multi.matching_list_filtered.si_Ball_filt_event_filt(:,1:6);

%% Load in SI scores (all neurons)
for ss=options.sessionSelect
    %100 bin spatial info score - A
    SI_score.A{ss} = session_vars{ss}.Place_cell{selectTrial(1)}.Spatial_Info.Spatial_Info(8,:);
    %100 bin spatial info score - B
    SI_score.B{ss} = session_vars{ss}.Place_cell{selectTrial(2)}.Spatial_Info.Spatial_Info(8,:);
end


%% For all neurons matching day-2-day, get S.I. scores for A and B trials

%get cell matrix with matching neurons first
%neighboring sessions (get matching indices
for ii=1:5
    matching_ROI_neighbor_idx{ii,ii+1} = matchROI_all(sum(~isnan(matchROI_all(:,[ii, ii+1])),2) == 2,[ii, ii+1]);
end

%relative to day 1
for ii=2:6
    matching_ROI_d1_idx{1,ii} = matchROI_all(sum(~isnan(matchROI_all(:,[1, ii])),2) == 2,[1, ii]);
end

%extract the S.I. score for ROIs matching relative to D1
for ii=2:6
    SI_scores_rel_d1.A{ii}(:,1) = SI_score.A{1}(matching_ROI_d1_idx{ii}(:,1));
    SI_scores_rel_d1.A{ii}(:,2) = SI_score.A{ii}(matching_ROI_d1_idx{ii}(:,2));
    
    SI_scores_rel_d1.B{ii}(:,1) = SI_score.B{1}(matching_ROI_d1_idx{ii}(:,1));
    SI_scores_rel_d1.B{ii}(:,2) = SI_score.B{ii}(matching_ROI_d1_idx{ii}(:,2));
end

%% Plot scatters - all neurons - matching A vs. B - D1 vs. D3

figure
subplot(1,2,1)
hold on
title('Neurons on D1 that match D3 - A  vs. B D1')
axis square
xlim([0 0.2])
ylim([0 0.2])
scatter(SI_scores_rel_d1.A{1, 3}(:,1), SI_scores_rel_d1.B{1, 3}(:,1),14,'filled')

subplot(1,2,2)
hold on
title('Neurons on D1 that match D3 - A  vs. B D3')
axis square
xlim([0 0.2])
ylim([0 0.2])
scatter(SI_scores_rel_d1.A{1, 3}(:,2), SI_scores_rel_d1.B{1, 3}(:,2),14,'filled')


%% Plot scatters - all neurons - matching A vs. A; B vs. B

day_sel = 6;

figure
subplot(1,2,1)
hold on
title('Neurons on D1 that match D1 -A')
axis square
xlim([0 0.2])
ylim([0 0.2])
scatter(SI_scores_rel_d1.A{1, day_sel}(:,1), SI_scores_rel_d1.A{1, day_sel}(:,2),14,'filled')
%plot unity line
plot([0 0.2],[0 0.2], 'k--');

subplot(1,2,2)
hold on
title('Neurons on D1 that match D1 - B')
axis square
xlim([0 0.2])
ylim([0 0.2])
scatter(SI_scores_rel_d1.B{1, day_sel}(:,1), SI_scores_rel_d1.B{1, day_sel}(:,2),14,'filled')
%plot unity line
plot([0 0.2],[0 0.2], 'k--');

%% Scatter spatial information scores

figure
for ss=sessionSelect
    subplot(1,6,ss)
    hold on
    xlim([0 0.2])
    ylim([0 0.2])
    axis square
    scatter(SI_score.A{ss}, SI_score.B{ss},14,'filled')
end

end

