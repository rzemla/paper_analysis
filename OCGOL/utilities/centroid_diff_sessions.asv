function [outputArg1,outputArg2] = centroid_diff_sessions(learn_comb_data,recall_comb_data,reg_learn, reg_recall)
%% Get all A/allB A&B


%% Set up variables

%learning
tuned_log_learning = learn_comb_data.tuned_log_learning;
pf_vector_max_learning = learn_comb_data.pf_vector_max_learning;

%recall
tuned_log_recall = recall_comb_data.tuned_log_recall;
pf_vector_max_recall = recall_comb_data.pf_vector_max_recall;


%% Compare the distance to start of reward zone A and reward zone B for learning/recall
%LEARNING SETS
%for each animal
options.selectTrial = [4,5];
options.sessionSelect = 1:6;

for aa=1:3
    [theta_reward_zone_learn{aa}] = A_B_angle_diff_across_sessions_reward_zones(reg_learn,tuned_log_learning,pf_vector_max_learning,aa,options);
end


%RECALL SETS
options.selectTrial = [1,2];
options.sessionSelect = 1:7;
for aa =1:4
    [theta_reward_zone_recall{aa}] = A_B_angle_diff_across_sessions_reward_zones(reg_recall,tuned_log_recall,pf_vector_max_recall,aa,options);
end


%% Parse session in distance into day distance based on imaging schedule for each animal

%merge animals according to day distance for learning (diff from rew zone
%for each 
%LEARNING MERGE FOR DAY ALIGNMENT
[merge_theta_rew_zones_learn_days] = arrange_learning_session_days_rew_dist(theta_reward_zone_learn)

%RECALL MERGE FOR DAY ALIGNMENT
[merge_theta_rew_zones_recall_days] = arrange_recall_session_rew_dist_days(theta_reward_zone_recall)

%merge across all animals for learning
for dd=1:6
    merge_theta_learn_rew_zone_days_all.A{dd} = cell2mat(merge_theta_rew_zones_learn_days.A.rewA(:,dd));
    merge_theta_learn_rew_zone_days_all.B{dd} = cell2mat(merge_theta_rew_zones_learn_days.B.rewB(:,dd));
end

%merge across all animals for recall
for dd=1:9
    merge_theta_recall_rew_zone_days_all.A{dd} = cell2mat(merge_theta_rew_zones_recall_days.A.rewA(:,dd));
    merge_theta_recall_rew_zone_days_all.B{dd} = cell2mat(merge_theta_rew_zones_recall_days.B.rewB(:,dd));
end


%% Get mean for learn and recall of centroid difference relative to reward zone
% A
%get mean
mean_rew_zone.learn.A = cellfun(@nanmean, merge_theta_learn_rew_zone_days_all.A(2:end));
mean_rew_zone.recall.A = cellfun(@nanmean, merge_theta_recall_rew_zone_days_all.A(2:6));


%number of neurons
nb_rew_zone.learn.A = cellfun(@(x) size(x,1), merge_theta_learn_rew_zone_days_all.A);
nb_rew_zone.recall.A = cellfun(@(x) size(x,1), merge_theta_recall_rew_zone_days_all.A);

%sem
sem_rew_zone.learn.A = cellfun(@nanstd, merge_theta_learn_rew_zone_days_all.A(2:end))./sqrt(nb_rew_zone.learn.A(2:end));
sem_rew_zone.recall.A = cellfun(@nanstd,  merge_theta_recall_rew_zone_days_all.A(2:6))./sqrt(nb_rew_zone.recall.A(2:6));

% B
%get mean
mean_rew_zone.learn.B = cellfun(@nanmean, merge_theta_learn_rew_zone_days_all.B(2:end));
mean_rew_zone.recall.B = cellfun(@nanmean, merge_theta_recall_rew_zone_days_all.B(2:6));

%number of neurons
nb_rew_zone.learn.B = cellfun(@(x) size(x,1), merge_theta_learn_rew_zone_days_all.B);
nb_rew_zone.recall.B = cellfun(@(x) size(x,1), merge_theta_recall_rew_zone_days_all.B);

%sem
sem_rew_zone.learn.B = cellfun(@nanstd, merge_theta_learn_rew_zone_days_all.B(2:end))./sqrt(nb_rew_zone.learn.B(2:end));
sem_rew_zone.recall.B = cellfun(@nanstd,  merge_theta_recall_rew_zone_days_all.B(2:6))./sqrt(nb_rew_zone.recall.B(2:6));

%% Make assignments based on tuning state 
%Calculate difference relative to D1 (for A, for B, for A vs.B - TS tuned/ min 5 events)

%get TS tuned - select only TS positive neurons in each reg column for A/B and A
%and B neurons

%LEARNING SETS
%for each animal
options.selectTrial = [4,5];
options.sessionSelect = 1:6;
for aa =1:3 %for each animal
    [theta_learn{aa}] = A_B_angle_diff_across_sessions(reg_learn,tuned_log_learning,pf_vector_max_learning,aa,options);
end

%RECALL SETS
options.selectTrial = [1,2];
options.sessionSelect = 1:7;
for aa =1:4
    [theta_recall{aa}] = A_B_angle_diff_across_sessions(reg_recall,tuned_log_recall,pf_vector_max_recall,aa,options);
end

%% Parse sessions by day duration relative to day 1 (since not all sessions were recording as sequential days) 
%arrange the angles to match relative to d1..d9 for recall
[merge_theta_learn_days] = arrange_learning_session_days(theta_learn);

%arrange the angles to match relative to d1..d9 for recall
[merge_theta_recall_days] = arrange_recall_session_days(theta_recall);


%merge across all animals for learning
for dd=1:6
    merge_theta_learn_days_all.A{dd} = cell2mat(merge_theta_learn_days.A(:,dd)');
    merge_theta_learn_days_all.B{dd} = cell2mat(merge_theta_learn_days.B(:,dd)');
end

%merge across all animals for recall
for dd=1:9
    merge_theta_recall_days_all.A{dd} = cell2mat(merge_theta_recall_days.A(:,dd)');
    merge_theta_recall_days_all.B{dd} = cell2mat(merge_theta_recall_days.B(:,dd)');
end


%% Get mean,sem for learning and recall relative to D1 A vs. A centroid distance and B vs. B centroid diff
% A
%get mean
mean_learn.A = cellfun(@nanmean, merge_theta_learn_days_all.A(2:end));
mean_recall.A = cellfun(@nanmean, merge_theta_recall_days_all.A(2:6));

%number of neurons
nb_learn.A = cellfun(@(x) size(x,2), merge_theta_learn_days_all.A);
nb_recall.A = cellfun(@(x) size(x,2), merge_theta_recall_days_all.A);

%sem
sem_learn.A = cellfun(@nanstd, merge_theta_learn_days_all.A(2:end))./sqrt(nb_learn.A(2:end));
sem_recall.A = cellfun(@nanstd, merge_theta_recall_days_all.A(2:6))./sqrt(nb_recall.A(2:6));

% B 
%mean
mean_learn.B = cellfun(@nanmean, merge_theta_learn_days_all.B(2:end));
mean_recall.B = cellfun(@nanmean, merge_theta_recall_days_all.B(2:6));

%number of neurons
nb_learn.B = cellfun(@(x) size(x,2), merge_theta_learn_days_all.B);
nb_recall.B = cellfun(@(x) size(x,2), merge_theta_recall_days_all.B);

%sem
sem_learn.B = cellfun(@nanstd, merge_theta_learn_days_all.B(2:end))./sqrt(nb_learn.B(2:end));
sem_recall.B = cellfun(@nanstd, merge_theta_recall_days_all.B(2:6))./sqrt(nb_recall.B(2:6));

%% Mann Whitney comp for each session - learning vs. recall 

for dd=2:6
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
xticks(1:5)
xlabel('Days since first session')
ylabel('Centroid diff. [cm]')

title('Centroid difference - A')
hb = bar([mean_learn.A(1:5)',  mean_recall.A(1:5)'],'BarWidth', 1);
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
xticks(1:5)
xlabel('Days since first session')
ylabel('Centroid diff. [cm]')
title('Centroid difference - B')
hb = bar([mean_learn.B(1:5)',  mean_recall.B(1:5)'],'BarWidth', 1);

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


%% Plot bar charts for change in distance to A vs. B reward zone for learning vs. recall

%rad to cm - ~196 = 2*pi 
cm2rad = (2*pi)./196;
%rad ticks that correspond to 0 25 50 cm
rad_ticks = [0 25 50].*cm2rad;

figure('Position',[2055 79 480 724])
subplot(2,1,1)
hold on
ylim([0 1.7])
%yticks(rad_ticks)
%yticklabels({'0','25','50'});
xticks(1:5)
xlabel('Days since first session')
ylabel('Centroid diff. change relative to rew zone ')

title('Centroid difference - A')
hb = bar([mean_rew_zone.learn.A(1:5)',  mean_rew_zone.recall.A(1:5)'],'BarWidth', 1);
pause(0.1)
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    if ib ==1
    errorbar(xData,mean_rew_zone.learn.A,sem_rew_zone.learn.A  ,'k.','linewidth',1.5)
    %color
    hb(ib).FaceColor =  [65,105,225]./255;
    %transparency
    hb(ib).FaceAlpha =  0.3;    
    
    elseif ib==2
      errorbar(xData,mean_rew_zone.recall.A,sem_rew_zone.recall.A,'k.','linewidth',1.5) 
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
xticks(1:5)
xlabel('Days since first session')
ylabel('Centroid diff. [cm]')
title('Centroid difference - B')
hb = bar([mean_rew_zone.learn.B(1:5)', mean_rew_zone.recall.B(1:5)'],'BarWidth', 1);

pause(0.1)
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    if ib ==1
    errorbar(xData,mean_rew_zone.learn.B,sem_rew_zone.learn.B,'k.','linewidth',1.5)
    %color
    hb(ib).FaceColor =  [220,20,60]./255;
    %transparency
    hb(ib).FaceAlpha =  0.3;    
    
    elseif ib==2
      errorbar(xData,mean_rew_zone.recall.B,sem_rew_zone.recall.B,'k.','linewidth',1.5) 
      %color
      hb(ib).FaceColor =  [220,20,60]./255;
      %transparency
      hb(ib).FaceAlpha =  1;
    end
end

set(gca,'Fontsize',16)
set(gca,'LineWidth',1.5)

legend([hb(1),hb(2)],{'Learning','Recall'},'location','northwest');

%% OLD
%{
%merge outputs into 1 cell
merge_theta_recall.A = [];
merge_theta_recall.B = [];

%for each animal
for aa=1:4
    merge_theta_recall.A = [merge_theta_recall.A ; theta_recall{aa}.A];
    merge_theta_recall.B = [merge_theta_recall.B ; theta_recall{aa}.B];
end

%merge across session
for ss=1:7
    theta_recall_merge.A{ss} = cell2mat(merge_theta_recall.A(:,ss)');
    theta_recall_merge.B{ss} = cell2mat(merge_theta_recall.B(:,ss)');
end

%     %merge outputs into 1 cell
%     merge_theta_learn.A = [];
%     merge_theta_learn.B = [];
%
%     %for each animal
%     for aa=1:3
%         merge_theta_learn.A = [merge_theta_learn.A ; theta_learn{aa}.A];
%         merge_theta_learn.B = [merge_theta_learn.B ; theta_learn{aa}.B];
%     end
%
%     %merge across session
%     for ss=1:6
%         theta_learn_merge.A{ss} = cell2mat(merge_theta_learn.A(:,ss)');
%         theta_learn_merge.B{ss} = cell2mat(merge_theta_learn.B(:,ss)');
%     end


%}
