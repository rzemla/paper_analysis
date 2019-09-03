%% Load  data from learn and recall animals
[CNMF_learn,reg_learn,reg_recall,PV_TC_corr_recall, perf_recall,PV_TC_corr_learning,perf_learning] = figure4_load_data();

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

%% Plot learning and recall PV correlation relative to day 1 on same plot

figure;
subplot(1,2,1)
hold on
title('PV A')
ylim([0 1])
%recall data
for ss = 1:size(PV_TC_corr_recall,2)
    plot([2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanPV_rel_d1.A,'g-')
end

%learning data
for ss = 1:size(PV_TC_corr_learning,2)
    plot([2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanPV_rel_d1.A,'r--')
end

subplot(1,2,2)
hold on
title('PV B')
ylim([0 1])
%recall data
for ss = 1:size(PV_TC_corr_recall,2)
    plot([2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanPV_rel_d1.B,'g-')
end

%learning data
for ss = 1:size(PV_TC_corr_learning,2)
    plot([2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanPV_rel_d1.B,'r--')
end


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

%% Performance plot
%make one 3-D matrix with fractional performance from all animals
%recall
for aa = 1:size(PV_TC_corr_recall,2)
    perf_recall_comb(:,:,aa) = perf_recall{aa}.ses_perf;
end

%learning
for aa = 1:size(PV_TC_corr_learning,2)
    perf_learning_comb(:,:,aa) = perf_learning{aa}.ses_perf;
end

%get means
mean_recall_perf = mean(perf_recall_comb,3);
mean_learning_perf = mean(perf_learning_comb,3);
%get stds
std_recall_perf = std(perf_recall_comb,0,3);
std_learning_perf = std(perf_learning_comb,0,3);
%get sems
sem_recall_perf = std_recall_perf./sqrt(size(PV_TC_corr_recall,2));
sem_learning_perf = std_learning_perf./size(PV_TC_corr_learning,2);

%plot
%recall - dash
%learning - solid
f = figure('Position',[2210 350 510 470]);
set(f,'color','w');
%combined - change color for each type
for tt=1
    %subplot(1,3,tt)
    hold on
    xlim([0.5 6.5])
    ylim([0 1.1])
    set(gca,'linewidth',2)
    set(gca,'FontSize',16)
yticks(0:0.2:1)
ylabel('Fraction of correct trials')
xticks((1:6))
xticklabels({'5A5B','5A5B','3A3B','3A3B','Random','Random'})
xtickangle(45)
%plot performance line (85%)
plot([0.5 6.5],[0.85 0.85],'k--','LineWidth',2);
%plot([1 2 3 6 7 8 9],mean_recall_perf(1,:),'m-')

    if tt ==1
        %all
        color_vec = [139, 0, 139]/255;
    elseif tt ==2
        %A
        color_vec = [65,105,225]/255;
    elseif tt==3
        %B
        color_vec = [ 220,20,60]/255;
    end
    %errorbar + mean  - learning
    errorbar([1:6],mean_learning_perf(tt,:),sem_learning_perf(tt,:),'Color', color_vec, 'LineStyle', '-','LineWidth',1.5)
    %errorbar + mean - recall
    %errorbar([1 2 3 6 7 8 9],mean_recall_perf(tt,:),sem_recall_perf(tt,:),'Color', color_vec, 'LineStyle', '-','LineWidth',1.5)
end

%save performance figure
disp('Saving performance figure ')
export_fig(f ,fullfile('G:\Figure_4_figures','learning_performance.png'),'-r300')

%plot([1:6],mean_learning_perf(1,:),'LineStyle','--','LineWidth',2,'Color', [139, 0, 139]/255)

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


%% Show outlines of neightboring day matching componenets (1 vs. 2) - learning
%which animal
learn_animal = 2;
% which two sessions to display
ses_comp = [1 3];

%remove soma parsed components from Coor_kp - coordinate outlines from CNMF
%ses_1
soma_keep{1} = CNMF_learn.removeROI_learn{learn_animal}{ses_comp(1)}.compSelect;
%ses 2
soma_keep{2} = CNMF_learn.removeROI_learn{learn_animal}{ses_comp(2)}.compSelect;
%coordinate outlines from both sessions
coor_keep{1} = CNMF_learn.CNMF_vars_learn{learn_animal}{ses_comp(1)}.Coor_kp(soma_keep{1});
coor_keep{2} = CNMF_learn.CNMF_vars_learn{learn_animal}{ses_comp(2)}.Coor_kp(soma_keep{2});

%extract filtered matching matrix for the two days
extract_match_idx = find(sum(~isnan(reg_learn{learn_animal}.registered.multi.assigned_filtered(:,ses_comp)),2) ==2);
%select session match matrix
vis_match_matrix = reg_learn{learn_animal}.registered.multi.assigned_filtered(extract_match_idx,ses_comp);

%plot BW outline of the component
f= figure('Position',[2100 150 1460 540]);
%set figure background color to white
set(f,'color','w');

subaxis(1,3,1, 'SpacingHorizontal', 0.0015,...
    'SpacingVertical',0.001,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
%subplot(1,3,1)
imagesc(CNMF_learn.templates_learn{learn_animal}{ses_comp(1)}.template);
hold on
axes(gca);
axis square
xticks(gca,[])
yticks(gca,[])
grayMap = brighten(gray,0.2);
colormap(gca,grayMap)

%plot all selected ROIs as green
for ROI = vis_match_matrix(:,1)'
    %plot componenet outline
    
    %plot(coor_keep{1}{ROI}(1,:),coor_keep{1}{ROI}(2,:),'g', 'LineWidth',1);
    f1= fill(coor_keep{1}{ROI}(1,:),coor_keep{1}{ROI}(2,:),'c', 'LineWidth',1,'EdgeColor','none');
    alpha(f1,0.5)
end

%subplot(1,3,2)
subaxis(1,3,2, 'SpacingHorizontal', 0.0015,...
    'SpacingVertical',0.001,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(CNMF_learn.templates_learn{learn_animal}{ses_comp(2)}.template);
hold on
axes(gca);
axis square
xticks(gca,[])
yticks(gca,[])
grayMap = brighten(gray,0.2);
colormap(gca,grayMap)
%plot all selected ROIs as green
for ROI = vis_match_matrix(:,2)'
    %plot componenet outline
    %plot(coor_keep{2}{ROI}(1,:),coor_keep{2}{ROI}(2,:),'m', 'LineWidth',1);
        f2= fill(coor_keep{2}{ROI}(1,:),coor_keep{2}{ROI}(2,:),'y', 'LineWidth',1,'EdgeColor','none');
    alpha(f2,0.5)
end

%combined
subaxis(1,3,3, 'SpacingHorizontal', 0.0015,...
    'SpacingVertical',0.001,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);

imagesc(CNMF_learn.templates_learn{learn_animal}{ses_comp(1)}.template);
hold on
axes(gca);
axis square
xticks(gca,[])
yticks(gca,[])
grayMap = brighten(gray,0.2);
colormap(gca,grayMap)

%plot all selected ROIs as green
for ROI = vis_match_matrix(:,1)'
    %plot componenet outline
    
    %plot(coor_keep{1}{ROI}(1,:),coor_keep{1}{ROI}(2,:),'g--', 'LineWidth',1);
     f1= fill(coor_keep{1}{ROI}(1,:),coor_keep{1}{ROI}(2,:),'y', 'LineWidth',1,'EdgeColor','none');
    alpha(f1,0.5)
end

for ROI = vis_match_matrix(:,2)'
    %plot componenet outline
    %plot(coor_keep{2}{ROI}(1,:),coor_keep{2}{ROI}(2,:),'m-', 'LineWidth',1);
f2= fill(coor_keep{2}{ROI}(1,:),coor_keep{2}{ROI}(2,:),'m','EdgeColor','none')
 alpha(f2,0.5)
end

%save figure 4a component match example learning
disp('Saving match ROIs STC ')
export_fig(f ,fullfile('G:\Figure_4_figures','componenent_matching.png'),'-r300')
