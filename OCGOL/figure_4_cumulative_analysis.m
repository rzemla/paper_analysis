%% Load  data from learn and recall animals

%skip loading of session var data
options.skipSesVars = 1;

tic;
[CNMF_learn,reg_learn,reg_recall,reg_recall_long,session_vars_learn,...
        session_vars_recall,session_vars_recall_long,...
        short_term_learn,short_term_recall,long_term_recall] = figure4_load_data(options);
toc;

%% Learning performance line plot

cumulative_performance_plot(short_term_learn,short_term_recall)

%% Plot A/B/AB/neither distributions for learning/recall across sessions

plot_fraction_tuned(short_term_learn.tuned_frac,short_term_recall.tuned_frac)

 
%% Combine STC matches across time relative to D1 and neighboring days (all animals into 1)
%also orient maps by day (T.S.) only tuned for now for learning and recall

%short term learning vs recall
combine_STC_plot_multi_animal(short_term_learn.TC_corr_match,short_term_recall.TC_corr_match)


%vs long term data (supplement - not enough animals)
combine_STC_plot_multi_animal_short_vs_long_term(short_term_learn.TC_corr_match,long_term_recall.TC_corr_match)

%% Tuning curve correlation for A&B tuned neurons across days
%this runs for T.S. neurons
%uses output from tc_corr_matching_neurons
%correlation done on curves normalized to each other

AandB_corr_rel_d1(short_term_learn.TC_corr_match,short_term_recall.TC_corr_match)

%% Centroid difference across sessions

centroid_diff_sessions(short_term_learn,short_term_recall,reg_learn, reg_recall)

%% Recurrence analysis

%short and long recall
recurrence_cum_analysis(short_term_recall,long_term_recall)

%learning only
recurrence_cum_analysis_learn(short_term_learn)

%% Examine spatial trajectories for across time

%trajectory_analysis(TC_corr_match_learning,TC_corr_match_recall)


%% Show outlines of neightboring day matching componenets (1 vs. 3) - learning
%modify this to show discarded matches

show_component_match(CNMF_learn,reg_learn)

%% Mean TC histograms across time during learning

%combine mean TC scores for each session from each animal

for ss=1:6
        meanTC_SCE_combined{ss} = [];
    for aa=1:3
        meanTC_SCE_combined{ss} = [meanTC_SCE_combined{ss}, SCE_learning{aa}.SCE{ss}.meanTC];
    end
end

figure
hold on
for ss=1:6
    ecdf(meanTC_SCE_combined{ss})
pause
end

%% SCE participation vs. TC score

matching_ROI_matrix = reg_learn{1, 1}.registered.multi.assigned_filtered;

%create matching A neuron SCE participation matrix
SCE_A_ROI_engage = zeros(size(matching_ROI_matrix,1), size(matching_ROI_matrix,2));

for ss=1:6
    assign_counts = SCE_part_A{ss}(matching_ROI_matrix(~isnan(matching_ROI_matrix(:,ss)),ss));
    SCE_A_ROI_engage(~isnan(matching_ROI_matrix(:,ss)),ss) = assign_counts
    SCE_A_ROI_engage(isnan(matching_ROI_matrix(:,ss)),ss) = nan;
end


%for each session 
for ss=1:6
    SCE_part_A{ss} = sum(SCE_learning{1, 1}.SCE{ss}.sce_activity.A,2)
end

figure
hold on
for ss=1:6
    plot(SCE_part_A{ss})
pause
end
legend();

%% Plot learning and recall TC correlation relative to day 1 on same plot

figure;
subplot(1,2,1)
hold on
title('TC A')
ylim([0 1])
%recall data
for ss = 1:size(PV_TC_corr_recall,2)
    plot([2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanTC_rel_d1.ts.A,'g-')
end

%learning data
for ss = 1:size(PV_TC_corr_learning,2)
    plot([2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanTC_rel_d1.ts.A,'r--')
end

subplot(1,2,2)
hold on
title('TC B')
ylim([0 1])
%recall data
for ss = 1:size(PV_TC_corr_recall,2)
    plot([2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanTC_rel_d1.ts.B,'g-')
end

%learning data
for ss = 1:size(PV_TC_corr_learning,2)
    plot([2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanTC_rel_d1.ts.B,'r--')
end

%% Extract PV/TC correlation data relative to D1

%recall - collect rel PV data
for ss = 1:size(PV_TC_corr_recall,2)
    %collect data from all animals into matrix
    %A PV
   PV_rel_d1.recall.A.all(ss,:) = PV_TC_corr_recall(ss).PV_TC_corr.meanPV_rel_d1.A;
   %B PV
   PV_rel_d1.recall.B.all(ss,:) = PV_TC_corr_recall(ss).PV_TC_corr.meanPV_rel_d1.B;
   %A TC
   TC_rel_d1.recall.A.all(ss,:) = PV_TC_corr_recall(ss).PV_TC_corr.meanTC_rel_d1.ts.A;
   %B TC
   TC_rel_d1.recall.B.all(ss,:) = PV_TC_corr_recall(ss).PV_TC_corr.meanTC_rel_d1.ts.B;   
end

%learning - collect rel PV data
for ss = 1:size(PV_TC_corr_learning,2)
    %collect data from all animals into matrix
    %A
    PV_rel_d1.learn.A.all(ss,:) = PV_TC_corr_learning(ss).PV_TC_corr.meanPV_rel_d1.A;
    %B
    PV_rel_d1.learn.B.all(ss,:) = PV_TC_corr_learning(ss).PV_TC_corr.meanPV_rel_d1.B;
    %A TC
    TC_rel_d1.learn.A.all(ss,:) = PV_TC_corr_learning(ss).PV_TC_corr.meanTC_rel_d1.ts.A;
    %B TC
    TC_rel_d1.learn.B.all(ss,:) = PV_TC_corr_learning(ss).PV_TC_corr.meanTC_rel_d1.ts.B;
end

%get mean
PV_rel_d1.learn.A.mean = nanmean(PV_rel_d1.learn.A.all,1);
PV_rel_d1.learn.B.mean = nanmean(PV_rel_d1.learn.B.all,1);

PV_rel_d1.recall.A.mean = nanmean(PV_rel_d1.recall.A.all,1);
PV_rel_d1.recall.B.mean = nanmean(PV_rel_d1.recall.B.all,1);

TC_rel_d1.learn.A.mean = nanmean(TC_rel_d1.learn.A.all,1);
TC_rel_d1.learn.B.mean = nanmean(TC_rel_d1.learn.B.all,1);

TC_rel_d1.recall.A.mean = nanmean(TC_rel_d1.recall.A.all,1);
TC_rel_d1.recall.B.mean = nanmean(TC_rel_d1.recall.B.all,1);

%get sem - learn
PV_rel_d1.learn.A.sem = nanstd(PV_rel_d1.learn.A.all,0,1)./sqrt(size(PV_TC_corr_learning,2));
PV_rel_d1.learn.B.sem = nanstd(PV_rel_d1.learn.B.all,0,1)./sqrt(size(PV_TC_corr_learning,2));
%recall
PV_rel_d1.recall.A.sem = nanstd(PV_rel_d1.recall.A.all,0,1)./sqrt(size(PV_TC_corr_recall,2));
PV_rel_d1.recall.B.sem = nanstd(PV_rel_d1.recall.B.all,0,1)./sqrt(size(PV_TC_corr_recall,2));

TC_rel_d1.learn.A.sem = nanstd(TC_rel_d1.learn.A.all,0,1)./sqrt(size(PV_TC_corr_learning,2));
TC_rel_d1.learn.B.sem = nanstd(TC_rel_d1.learn.B.all,0,1)./sqrt(size(PV_TC_corr_learning,2));
%recall
TC_rel_d1.recall.A.sem = nanstd(TC_rel_d1.recall.A.all,0,1)./sqrt(size(PV_TC_corr_recall,2));
TC_rel_d1.recall.B.sem = nanstd(TC_rel_d1.recall.B.all,0,1)./sqrt(size(PV_TC_corr_recall,2));

%% Plot learning and recall PV correlation relative to day 1 on same plot

%color vectors - blue, red, magenta
color_vec(1,:) = [139, 0, 139]/255;
color_vec(2,:) = [65,105,225]/255;
color_vec(3,:) = [ 220,20,60]/255;


f = figure('Position', [2011 321 1154 556]);
set(f,'color','w');
subplot(1,2,1)
hold on
axis square
title('PV correlation')
ylim([0 1])
xlim([0 9])
xticks(1:2:8)
ylabel('Correlation coef.')
xlabel('Sessions since first session')
set(gca,'Linewidth',2)
set(gca,'FontSize', 20)
%recall data
for ss = 1:size(PV_TC_corr_recall,2)
    %A
    rA = errorbar([1 2 5 6 7 8],PV_rel_d1.recall.A.mean,PV_rel_d1.recall.A.sem,'Color', color_vec(2,:), 'LineStyle', '-','LineWidth',1.5)
    %B
    rB = errorbar([1 2 5 6 7 8],PV_rel_d1.recall.B.mean,PV_rel_d1.recall.B.sem,'Color', color_vec(3,:), 'LineStyle', '-','LineWidth',1.5)
        
end

%learning data
for ss = 1:size(PV_TC_corr_learning,2)
    %A
    lA = errorbar([1:5],PV_rel_d1.learn.A.mean,PV_rel_d1.learn.A.sem,'Color', color_vec(2,:), 'LineStyle', '--','LineWidth',1.5)
    %B
    lB = errorbar([1:5],PV_rel_d1.learn.B.mean,PV_rel_d1.learn.B.sem,'Color', color_vec(3,:), 'LineStyle', '--','LineWidth',1.5)
end

legend([lA lB rA rB], {'Learning A','Learning B','Recall A','Recall B'},'Location','northeast')

%TC
subplot(1,2,2)
hold on
axis square
title('TC correlation')
ylim([0 1])
xlim([0 9])
xticks(1:2:8)
%ylabel('Correlation coef.')
xlabel('Sessions since first session')
set(gca,'Linewidth',2)
set(gca,'FontSize', 20)
for ss = 1:size(PV_TC_corr_recall,2)
    %A
    rA = errorbar([1 2 5 6 7 8],TC_rel_d1.recall.A.mean,TC_rel_d1.recall.A.sem,'Color', color_vec(2,:), 'LineStyle', '-','LineWidth',1.5)
    %B
    rB = errorbar([1 2 5 6 7 8],TC_rel_d1.recall.B.mean,TC_rel_d1.recall.B.sem,'Color', color_vec(3,:), 'LineStyle', '-','LineWidth',1.5)
        
end

%learning data
for ss = 1:size(PV_TC_corr_learning,2)
    %A
    lA = errorbar([1:5],TC_rel_d1.learn.A.mean,TC_rel_d1.learn.A.sem,'Color', color_vec(2,:), 'LineStyle', '--','LineWidth',1.5)
    %B
    lB = errorbar([1:5],TC_rel_d1.learn.B.mean,TC_rel_d1.learn.B.sem,'Color', color_vec(3,:), 'LineStyle', '--','LineWidth',1.5)
end
legend([lA lB rA rB], {'Learning A','Learning B','Recall A','Recall B'},'Location','northeast')

%save performance figure
disp('Saving PV/TC correlation relative to D1 figure ')
export_fig(f ,fullfile('G:\Figure_4_figures','PV_TC_rel_D1.png'),'-r300')


%% Same day PV

figure;
hold on
title('Same day PV - at least 1 match to another session')
ylim([0 1])
%recall data
for ss = 1:size(PV_TC_corr_recall,2)
    plot([1 2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanPV_same_day.AB,'g-')
end

%learning data
for ss = 1:size(PV_TC_corr_learning,2)
    plot([1 2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanPV_same_day.AB,'r--')
end

%% Construct matrix with same day PV corr vs total performance
%combine across animals
%learning
for aa = 1:size(PV_TC_corr_learning,2)
    %mean PV over entire lap
    PV_corr_vs_perf.learn.meanPV(:,aa) = nanmean(PV_TC_corr_learning(aa).PV_TC_corr.PVcorr_all_same_day_diag,2)';
    %all performance
    PV_corr_vs_perf.learn.fracPerf(:,aa) = perf_learning{aa}.ses_perf(1,:)';
end

%recall
for aa = 1:size(PV_TC_corr_recall,2)
    %mean PV over entire lap
    PV_corr_vs_perf.recall.meanPV(:,aa) = nanmean(PV_TC_corr_recall(aa).PV_TC_corr.PVcorr_all_same_day_diag,2)';
    %all performance
    PV_corr_vs_perf.recall.fracPerf(:,aa) = perf_recall{aa}.ses_perf(1,:)';
end

%% Experiment with linear regression and confidence intervals

[r,p] = corrcoef(PV_corr_vs_perf.learn.fracPerf(:),PV_corr_vs_perf.learn.meanPV(:))

[r,p] = corrcoef(PV_corr_vs_perf.recall.fracPerf(:),PV_corr_vs_perf.recall.meanPV(:))

%[b,bint,~,~,stats] = regress(perf_learning{ss}.ses_perf(1,:)',[ones(6,1),nanmean(PV_TC_corr_learning(ss).PV_TC_corr.PVcorr_all_same_day_diag,2)])


%% Plot mean PV correlation (same day) as fraction of total correct scatter

figure
subplot(1,2,1)
hold on
xlim([0 1])
ylim([0 1])
axis('square')
title('Learning');
xlabel('Population correlation between trials');
ylabel('Fraction of total trials correct');
%for each animal
for ss=1:size(PV_TC_corr_learning,2)
    scatter(nanmean(PV_TC_corr_learning(ss).PV_TC_corr.PVcorr_all_same_day_diag,2)' ,perf_learning{ss}.ses_perf(1,:),'b')
end

subplot(1,2,2)
hold on
xlim([0 1])
ylim([0 1])
axis('square')
title('Recall');
xlabel('Population correlation between trials');
ylabel('Fraction of total trials correct');
for ss=1:size(PV_TC_corr_learning,2)
    scatter(nanmean(PV_TC_corr_recall(ss).PV_TC_corr.PVcorr_all_same_day_diag,2)' ,perf_recall{ss}.ses_perf(1,:),'b')
end

%% PV raster (whole track/diagonal) across learning

for ss=1:size(PV_TC_corr_learning,2)
    PV_diag_same_all_ses.learn(:,:,ss) = PV_TC_corr_learning(ss).PV_TC_corr.PVcorr_all_same_day_diag;
end

for ss=1:size(PV_TC_corr_recall,2)
    PV_diag_same_all_ses.recall(:,:,ss) = PV_TC_corr_recall(ss).PV_TC_corr.PVcorr_all_same_day_diag;
end
%means
PV_diag_same_mean.learn = mean(PV_diag_same_all_ses.learn,3);
PV_diag_same_mean.recall = mean(PV_diag_same_all_ses.recall,3);

f = figure('Position',[2170 240 1400 700]);
%set figure background color to white
set(f,'color','w');
subaxis(1,3,1,'SpacingHorizontal', 0.05,...
    'SpacingVertical',0.001,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1)
imagesc(PV_diag_same_mean.learn)
hold on
%plot odor, rewards locations
plot([10 10],[0 7],'k--', 'LineWidth',1.5)
plot([30 30],[0 7],'k--', 'LineWidth',1.5)
plot([70 70],[0 7],'k--', 'LineWidth',1.5)
caxis([0 1])
xticks([1 100])
title('Learning')
xticklabels({'0', '2'})
xlabel('Position [m]')
yticks(1:6)
ylabel('Learning session')
set(gca, 'TickLength', [0 0]);
set(gca,'FontSize',18)
set(gca,'linewidth',2)
colormap('jet')
%
subaxis(1,3,2,'SpacingHorizontal', 0.05,...
    'SpacingVertical',0.001,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1)
imagesc(PV_diag_same_mean.recall)
hold on
%plot odor, rewards locations
plot([10 10],[0 8],'k--', 'LineWidth',1.5)
plot([30 30],[0 8],'k--', 'LineWidth',1.5)
plot([70 70],[0 8],'k--', 'LineWidth',1.5)

caxis([0 1])
xticks([1 100])
xticklabels({'0', '2'})
title('Recall')
xlabel('Position [m]')
%ylabel('Recall session')
set(gca, 'TickLength', [0 0]);
set(gca,'FontSize',18)
set(gca,'linewidth',2)
colormap('jet')

subaxis(1,3,3,'SpacingHorizontal', 0.05,...
    'SpacingVertical',0.001,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1)
hold on
caxis([0 1])
colormap('jet')
c = colorbar('location','Westoutside',...
  'XTick',[-1 -0.5 0 0.5 1]);
axis off
c.Label.String = 'Correlation coef.';
c.FontSize = 20;
c.AxisLocation  = 'in';
set(gca,'FontSize',18)
set(gca,'linewidth',2)

disp('Saving match ROIs STC ')
export_fig(f ,fullfile('G:\Figure_4_figures','PV_across_track_raster.png'),'-r300')


%% SCE rate scatter plots and cdfs
%create blank entries for SCE related data for now
if 0
    SCE_rate_plots(session_vars_learn,session_vars_recall,SCE_learning,SCE_recall, perf_learning,perf_recall)
end
