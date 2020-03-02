function [diff_angle] = centroid_diff_sessions_update(short_term_learn,short_term_recall,reg_learn, reg_recall, excl_day_combined_day_nan,options)


%% Define number of animals for each experiment

%use performance struct as input
%short term learn
nb_st_learn = size(short_term_learn.perf,2);

%short term recall
nb_st_recall = size(short_term_recall.perf,2);

%% Set up variables

%learning
tuned_log_learning = short_term_learn.tuned_log;
pf_vector_max_learning = short_term_learn.pf_vector_max;

%recall
tuned_log_recall = short_term_recall.tuned_log;
pf_vector_max_recall = short_term_recall.pf_vector_max;

%% Make assignments based on tuning state - all cross tuning sessions
%Calculate difference relative to D1 (for A, for B, for A vs.B - TS tuned/ min 5 events)

%get TS tuned - select only TS positive neurons in each reg column for A/B and A
%and B neurons

%run for ts

%LEARNING SETS
%for each animal
%all A and B trial regardless if correct
options.selectTrial = [4,5];
for aa =1:nb_st_learn %for each animal
    [theta_learn{aa}] = A_B_angle_diff_across_sessions(reg_learn,tuned_log_learning,pf_vector_max_learning,aa,options);
end

%RECALL SETS
%only correct A and B trials
options.selectTrial = [1,2];
for aa =1:nb_st_recall
    [theta_recall{aa}] = A_B_angle_diff_across_sessions(reg_recall,tuned_log_recall,pf_vector_max_recall,aa,options);
end

%LT RECALL here (30 days)


%% Filter and rearrange by days for short term learn and recall using updated approach

%ST LEARN rearrangement and filter
%for each animal
for aa=1:nb_st_learn
    %rearrange based on days
    for dd=1:9
        %index corresponding to day (#1 is st_learn)
        session_day_idx = find(excl_day_combined_day_nan{aa,1}(2,:) == dd);
        if ~isempty(session_day_idx)
            %get combined, A, and B performance (animal x day)
            %A
            theta_learn_day_filt_rel_d1.A{aa,dd} = theta_learn{aa}.A{1, session_day_idx};
            %B
            theta_learn_day_filt_rel_d1.B{aa,dd} = theta_learn{aa}.B{1, session_day_idx};
        else
            theta_learn_day_filt_rel_d1.A{aa,dd} = nan;
            theta_learn_day_filt_rel_d1.B{aa,dd} = nan;
        end
    end
end

%ST RECALL rearrangement and filter
%for each animal
for aa=1:nb_st_recall
    %rearrange based on days
    for dd=1:9
        %index corresponding to day (#2 is st_recall)
        session_day_idx = find(excl_day_combined_day_nan{aa,2}(2,:) == dd);
        if ~isempty(session_day_idx)
            %A
            theta_recall_day_filt_rel_d1.A{aa,dd} = theta_recall{aa}.A{1, session_day_idx};
            %B
            theta_recall_day_filt_rel_d1.B{aa,dd} = theta_recall{aa}.B{1, session_day_idx};
        else
            theta_recall_day_filt_rel_d1.A{aa,dd} = nan;
            theta_recall_day_filt_rel_d1.B{aa,dd} = nan;
        end
    end
end

%% Get mean for each animal --> get mean of means and sem

%st learn
%A
diff_angle.st_learn.animal.mean.A = cellfun(@nanmean,theta_learn_day_filt_rel_d1.A);
%mean of means
diff_angle.st_learn.animal.mean_mean.A = nanmean(diff_angle.st_learn.animal.mean.A,1);
%sem of means
diff_angle.st_learn.animal.mean_sem.A = nanstd(diff_angle.st_learn.animal.mean.A,0,1)./sqrt(sum(~isnan(diff_angle.st_learn.animal.mean.A),1));
%raw
diff_angle.st_learn.animal.raw.A = diff_angle.st_learn.animal.mean.A;

%B
diff_angle.st_learn.animal.mean.B = cellfun(@nanmean,theta_learn_day_filt_rel_d1.B);
%mean of means
diff_angle.st_learn.animal.mean_mean.B = nanmean(diff_angle.st_learn.animal.mean.B,1);
%sem of means
diff_angle.st_learn.animal.mean_sem.B = nanstd(diff_angle.st_learn.animal.mean.B,0,1)./sqrt(sum(~isnan(diff_angle.st_learn.animal.mean.B),1));
%raw
diff_angle.st_learn.animal.raw.B = diff_angle.st_learn.animal.mean.B;

%st recall
%A
diff_angle.st_recall.animal.mean.A = cellfun(@nanmean,theta_recall_day_filt_rel_d1.A);
%mean of means
diff_angle.st_recall.animal.mean_mean.A = nanmean(diff_angle.st_recall.animal.mean.A,1);
%sem of means
diff_angle.st_recall.animal.mean_sem.A = nanstd(diff_angle.st_recall.animal.mean.A,0,1)./sqrt(sum(~isnan(diff_angle.st_recall.animal.mean.A),1));
%raw
diff_angle.st_recall.animal.raw.A = diff_angle.st_recall.animal.mean.A;


%B
diff_angle.st_recall.animal.mean.B = cellfun(@nanmean,theta_recall_day_filt_rel_d1.B);
%mean of means
diff_angle.st_recall.animal.mean_mean.B = nanmean(diff_angle.st_recall.animal.mean.B,1);
%sem of means
diff_angle.st_recall.animal.mean_sem.B = nanstd(diff_angle.st_recall.animal.mean.B,0,1)./sqrt(sum(~isnan(diff_angle.st_recall.animal.mean.B),1));
%raw
diff_angle.st_recall.animal.raw.B = diff_angle.st_recall.animal.mean.B;

%% Merge all neurons from all animals - get mean and sem

%for each day, merge cells into 1 - all neurons pooled from all animals on
%that day
for dd=1:9
    %A st learn
    pooled_diff_angle.st_learn.A{dd} = cell2mat(theta_learn_day_filt_rel_d1.A(:,dd)');
    %B st learn
    pooled_diff_angle.st_learn.B{dd} = cell2mat(theta_learn_day_filt_rel_d1.B(:,dd)');
    
    %A st recall
    pooled_diff_angle.st_recall.A{dd} = cell2mat(theta_recall_day_filt_rel_d1.A(:,dd)');
    %B st recall
    pooled_diff_angle.st_recall.B{dd} = cell2mat(theta_recall_day_filt_rel_d1.B(:,dd)');
end

%% Calculate mean and sem for pooled neurons

%ST learn
%A mean
diff_angle.st_learn.pooled.mean.A = cellfun(@nanmean, pooled_diff_angle.st_learn.A);
%A sem
diff_angle.st_learn.pooled.sem.A = cellfun(@(x) nanstd(x,0,2),pooled_diff_angle.st_learn.A)./...
sqrt(cellfun(@sum,cellfun(@(x) ~isnan(x),pooled_diff_angle.st_learn.A,'UniformOutput',false)));
%raw
diff_angle.st_learn.pooled.raw.A = pooled_diff_angle.st_learn.A;

%B mean
diff_angle.st_learn.pooled.mean.B = cellfun(@nanmean, pooled_diff_angle.st_learn.B);
%B sem
diff_angle.st_learn.pooled.sem.B = cellfun(@(x) nanstd(x,0,2),pooled_diff_angle.st_learn.B)./...
sqrt(cellfun(@sum,cellfun(@(x) ~isnan(x),pooled_diff_angle.st_learn.B,'UniformOutput',false)));
%raw
diff_angle.st_learn.pooled.raw.B = pooled_diff_angle.st_learn.B;


%ST recall
%A mean
diff_angle.st_recall.pooled.mean.A = cellfun(@nanmean, pooled_diff_angle.st_recall.A);
%A sem
diff_angle.st_recall.pooled.sem.A = cellfun(@(x) nanstd(x,0,2),pooled_diff_angle.st_recall.A)./...
sqrt(cellfun(@sum,cellfun(@(x) ~isnan(x),pooled_diff_angle.st_recall.A,'UniformOutput',false)));
%raw
diff_angle.st_recall.pooled.raw.A = pooled_diff_angle.st_recall.A;

%B mean
diff_angle.st_recall.pooled.mean.B = cellfun(@nanmean, pooled_diff_angle.st_recall.B);
%B sem
diff_angle.st_recall.pooled.sem.B = cellfun(@(x) nanstd(x,0,2),pooled_diff_angle.st_recall.B)./...
sqrt(cellfun(@sum,cellfun(@(x) ~isnan(x),pooled_diff_angle.st_recall.B,'UniformOutput',false)));
%raw
diff_angle.st_recall.pooled.raw.B = pooled_diff_angle.st_recall.B;

%% QC check by animal
figure

subplot(1,2,1)
hold on
title('Animal - A')
%st learn
errorbar((1:7),diff_angle.st_learn.animal.mean_mean.A(1:7),diff_angle.st_learn.animal.mean_sem.A(1:7),'--')

%st recall
errorbar(([1:3,6:9]),diff_angle.st_recall.animal.mean_mean.A([1:3,6:9]),diff_angle.st_recall.animal.mean_sem.A([1:3,6:9]),'-')

subplot(1,2,2)
hold on
title('Animal - B')
%st learn
errorbar((1:7),diff_angle.st_learn.animal.mean_mean.B(1:7),diff_angle.st_learn.animal.mean_sem.B(1:7),'--')

%st recall
errorbar(([1:3,6:9]),diff_angle.st_recall.animal.mean_mean.B([1:3,6:9]),diff_angle.st_recall.animal.mean_sem.B([1:3,6:9]),'-')


%% Look at distribution of centroid differences in histograms - consistent with bar charts
%Relative to reward zones

%{

%A
subplot_indices = [1:8; 9:16];
%learning
figure
for ii=2:9
    subplot(2,8,subplot_indices(1,ii-1))
hold on
ylim([0 0.7])
xlim([0 pi])
title('Learning - A');
ylabel('Normalized density')
xlabel('Angular diff. [rad]')
    histogram(merge_theta_learn_days_all.A{ii},'Normalization','probability','BinEdges',0:(pi/9):pi)
end
%recall
for ii=2:9
    subplot(2,8,subplot_indices(2,ii-1))
    hold on
    ylim([0 0.7])
    xlim([0 pi])
    title('Recall - A')
    ylabel('Normalized density')
    xlabel('Angular diff. [rad]')
    histogram(merge_theta_recall_days_all.A{ii},'Normalization','probability','BinEdges',0:(pi/9):pi)
end

%% B
subplot_indices = [1:8; 9:16];
%learning
figure
for ii=2:9
    subplot(2,8,subplot_indices(1,ii-1))
hold on
ylim([0 0.7])
xlim([0 pi])
title('Learning - B');
ylabel('Normalized density')
xlabel('Angular diff. [rad]')
    histogram(merge_theta_learn_days_all.B{ii},'Normalization','probability','BinEdges',0:(pi/9):pi)
end
%recall
for ii=2:9
    subplot(2,8,subplot_indices(2,ii-1))
    hold on
    ylim([0 0.7])
    xlim([0 pi])
    title('Recall - B')
    ylabel('Normalized density')
    xlabel('Angular diff. [rad]')
    histogram(merge_theta_recall_days_all.B{ii},'Normalization','probability','BinEdges',0:(pi/9):pi)
end

%% Get mean,sem for learning and recall relative to D1 A vs. A centroid distance and B vs. B centroid diff
% A
%get mean
mean_learn.A = cellfun(@nanmean, merge_theta_learn_days_all.A(2:end));
mean_recall.A = cellfun(@nanmean, merge_theta_recall_days_all.A(2:end));

%number of neurons
nb_learn.A = cellfun(@(x) size(x,2), merge_theta_learn_days_all.A);
nb_recall.A = cellfun(@(x) size(x,2), merge_theta_recall_days_all.A);

%sem
sem_learn.A = cellfun(@nanstd, merge_theta_learn_days_all.A(2:end))./sqrt(nb_learn.A(2:end));
sem_recall.A = cellfun(@nanstd, merge_theta_recall_days_all.A(2:end))./sqrt(nb_recall.A(2:end));

% B 
%mean
mean_learn.B = cellfun(@nanmean, merge_theta_learn_days_all.B(2:end));
mean_recall.B = cellfun(@nanmean, merge_theta_recall_days_all.B(2:end));

%number of neurons
nb_learn.B = cellfun(@(x) size(x,2), merge_theta_learn_days_all.B);
nb_recall.B = cellfun(@(x) size(x,2), merge_theta_recall_days_all.B);

%sem
sem_learn.B = cellfun(@nanstd, merge_theta_learn_days_all.B(2:end))./sqrt(nb_learn.B(2:end));
sem_recall.B = cellfun(@nanstd, merge_theta_recall_days_all.B(2:end))./sqrt(nb_recall.B(2:end));

%% Statistics - Mann Whitney comp for each session - learning vs. recall 

for dd=2:9
    %learn vs. recall Mann Whitney U - A
    [pA(dd),~,~] = ranksum(merge_theta_learn_days_all.A{dd},merge_theta_recall_days_all.A{dd});
    
    %learn vs. recall Mann Whitney U - B
    [pB(dd),~,~] = ranksum(merge_theta_learn_days_all.B{dd},merge_theta_recall_days_all.B{dd});
end

%% Plot bar charts for learn/recall A vs learn/recall B comparison

%rad to cm - ~196 = 2*pi 
cm2rad = (2*pi)./196;
%rad ticks that correspond to 0 25 50 cm
rad_ticks = [0 25 50].*cm2rad;

figure('Position',[2055 79 480 724])
subplot(2,1,1)
hold on
ylim([0 1.7])
yticks(rad_ticks)
yticklabels({'0','25','50'});
xticks(1:8)
xlabel('Days since first session')
ylabel('Centroid diff. [cm]')

title('Centroid difference - A')
hb = bar([mean_learn.A(1:end)',  mean_recall.A(1:end)'],'BarWidth', 1);
pause(0.1)
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    if ib ==1
    errorbar(xData,mean_learn.A,sem_learn.A,'k.','linewidth',1.5)
    %color
    hb(ib).FaceColor =  [65,105,225]./255;
    %transparency
    hb(ib).FaceAlpha =  0.3;    
    
    elseif ib==2
      errorbar(xData,mean_recall.A,sem_recall.A,'k.','linewidth',1.5) 
      %color
      hb(ib).FaceColor =  [65,105,225]./255;
      %transparency
      hb(ib).FaceAlpha =  1;
    end
end

set(gca,'Fontsize',16)
set(gca,'LineWidth',1.5)

legend([hb(1),hb(2)],{'Learning','Recall'},'location','northwest');

subplot(2,1,2)
hold on
ylim([0 1.7])
yticks(rad_ticks)
yticklabels({'0','25','50'});
xticks(1:8)
xlabel('Days since first session')
ylabel('Centroid diff. [cm]')
title('Centroid difference - B')
hb = bar([mean_learn.B(1:end)',  mean_recall.B(1:end)'],'BarWidth', 1);

pause(0.1)
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    if ib ==1
    errorbar(xData,mean_learn.B,sem_learn.B,'k.','linewidth',1.5)
    %color
    hb(ib).FaceColor =  [220,20,60]./255;
    %transparency
    hb(ib).FaceAlpha =  0.3;    
    
    elseif ib==2
      errorbar(xData,mean_recall.B,sem_recall.B,'k.','linewidth',1.5) 
      %color
      hb(ib).FaceColor =  [220,20,60]./255;
      %transparency
      hb(ib).FaceAlpha =  1;
    end
end

set(gca,'Fontsize',16)
set(gca,'LineWidth',1.5)

legend([hb(1),hb(2)],{'Learning','Recall'},'location','northwest');

%}

end
