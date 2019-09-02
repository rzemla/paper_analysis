%% Load in PV and TC correlation data from recall animals

%cross session directories (recall experiments)
cross_dirs_recall = {'G:\OCGOL_stability_recall\I46\crossSession',...
    'G:\OCGOL_stability_recall\I45_RT\crossSession',...
    'G:\OCGOL_stability_recall\I42L_1\crossSession',...
    'G:\OCGOL_stability_recall\I42R_1\crossSession'};

%cross session directories (recall experiments)
cross_dirs_learning = {'G:\OCGOL_learning_short_term\I56_RTLS\crossSession',...
    'G:\OCGOL_learning_short_term\I57_RTLS\crossSession',...
    'G:\OCGOL_learning_short_term\I57_LT\crossSession'};

%read in recall data
for ss=1:size(cross_dirs_recall,2)
    %TC/PV correlations
    PV_TC_corr_recall(ss) = load(fullfile(cross_dirs_recall{ss},'PV_TC_corr.mat'));
    %performance data
    perf_recall{ss} = load(fullfile(cross_dirs_recall{ss},'ses_perf.mat'));
    
end

%read in recall data
for ss=1:size(cross_dirs_learning,2)
    %TC/PV correlations
 PV_TC_corr_learning(ss) = load(fullfile(cross_dirs_learning{ss},'PV_TC_corr.mat'));
     %performance data
    perf_learning{ss} = load(fullfile(cross_dirs_learning{ss},'ses_perf.mat'));
end

%% Plot learning and recall TC correlation relative to day 1 on same plot

figure;
subplot(1,2,1)
hold on
title('TC A')
ylim([0 1])
%recall data
for ss = 1:size(cross_dirs_recall,2)
    plot([2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanTC_rel_d1.ts.A,'g-')
end

%learning data
for ss = 1:size(cross_dirs_learning,2)
    plot([2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanTC_rel_d1.ts.A,'r--')
end

subplot(1,2,2)
hold on
title('TC B')
ylim([0 1])
%recall data
for ss = 1:size(cross_dirs_recall,2)
    plot([2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanTC_rel_d1.ts.B,'g-')
end

%learning data
for ss = 1:size(cross_dirs_learning,2)
    plot([2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanTC_rel_d1.ts.B,'r--')
end

%% Plot learning and recall PV correlation relative to day 1 on same plot

figure;
subplot(1,2,1)
hold on
title('PV A')
ylim([0 1])
%recall data
for ss = 1:size(cross_dirs_recall,2)
    plot([2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanPV_rel_d1.A,'g-')
end

%learning data
for ss = 1:size(cross_dirs_learning,2)
    plot([2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanPV_rel_d1.A,'r--')
end

subplot(1,2,2)
hold on
title('PV B')
ylim([0 1])
%recall data
for ss = 1:size(cross_dirs_recall,2)
    plot([2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanPV_rel_d1.B,'g-')
end

%learning data
for ss = 1:size(cross_dirs_learning,2)
    plot([2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanPV_rel_d1.B,'r--')
end


%% Same day PV

figure;
hold on
title('Same day PV - at least 1 match to another session')
ylim([0 1])
%recall data
for ss = 1:size(cross_dirs_recall,2)
    plot([1 2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanPV_same_day.AB,'g-')
end

%learning data
for ss = 1:size(cross_dirs_learning,2)
    plot([1 2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanPV_same_day.AB,'r--')
end

%% Construct matrix with same day PV corr vs total performance
%combine across animals
%learning
for aa = 1:size(cross_dirs_learning,2)
    %mean PV over entire lap
    PV_corr_vs_perf.learn.meanPV(:,aa) = nanmean(PV_TC_corr_learning(aa).PV_TC_corr.PVcorr_all_same_day_diag,2)';
    %all performance
    PV_corr_vs_perf.learn.fracPerf(:,aa) = perf_learning{aa}.ses_perf(1,:)';
end

%recall
for aa = 1:size(cross_dirs_recall,2)
    %mean PV over entire lap
    PV_corr_vs_perf.recall.meanPV(:,aa) = nanmean(PV_TC_corr_recall(aa).PV_TC_corr.PVcorr_all_same_day_diag,2)';
    %all performance
    PV_corr_vs_perf.recall.fracPerf(:,aa) = perf_recall{aa}.ses_perf(1,:)';
end

%% Performance plot
%make one 3-D matrix with fractional performance from all animals
%recall
for aa = 1:size(cross_dirs_recall,2)
    perf_recall_comb(:,:,aa) = perf_recall{aa}.ses_perf;
end

%learning
for aa = 1:size(cross_dirs_learning,2)
    perf_learning_comb(:,:,aa) = perf_learning{aa}.ses_perf;
end

%get means
mean_recall_perf = mean(perf_recall_comb,3);
mean_learning_perf = mean(perf_learning_comb,3);
%get stds
std_recall_perf = std(perf_recall_comb,0,3);
std_learning_perf = std(perf_learning_comb,0,3);
%get sems
sem_recall_perf = std_recall_perf./sqrt(size(cross_dirs_recall,2));
sem_learning_perf = std_learning_perf./size(cross_dirs_learning,2);

%plot
%recall - dash
%learning - solid
figure;
hold on;
%xaxis(da
plot([1 2 3 6 7 8 9],mean_recall_perf(1,:),'m-')

plot([1:6],mean_learning_perf(1,:),'LineStyle','--','LineWidth',2,'Color', [139, 0, 139]/255)
%errorbar 
errorbar(xData,frac_tuned_each_mean.si',frac_tuned_each_sem.si,'k.')

%shaded plot parameters
% s = shadedErrorBar(1:6,squeeze(perf_learning_comb(1,:,:))',{@mean,@std},'lineprops','-','transparent',true,'patchSaturation',0.20);
% set(s.edge,'LineWidth',1,'LineStyle','-','Color',[[139, 0, 139]/255, 0.2]) %last value add transparency value
% s.mainLine.LineWidth = 2;
% s.mainLine.Color = [139, 0, 139]/255;
% s.patch.FaceColor = [139, 0, 139]/255;

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
for ss=1:size(cross_dirs_learning,2)
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
for ss=1:size(cross_dirs_learning,2)
    scatter(nanmean(PV_TC_corr_recall(ss).PV_TC_corr.PVcorr_all_same_day_diag,2)' ,perf_recall{ss}.ses_perf(1,:),'b')
end


for ss=1:size(cross_dirs_learning,2)
    scatter(PV_TC_corr_learning(ss).PV_TC_corr.meanPV_rel_d1.B,perf_learning{ss}.ses_perf(1,:),'r')
end

