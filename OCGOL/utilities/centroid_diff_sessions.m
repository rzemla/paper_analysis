function [outputArg1,outputArg2] = centroid_diff_sessions(short_term_learn,short_term_recall,reg_learn, reg_recall)
%% Get all A/B A&B


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

%LEARNING SETS
%for each animal
options.selectTrial = [4,5];
for aa =1:6 %for each animal
    [theta_learn{aa}] = A_B_angle_diff_across_sessions(reg_learn,tuned_log_learning,pf_vector_max_learning,aa,options);
end

%RECALL SETS
options.selectTrial = [1,2];
for aa =1:5
    [theta_recall{aa}] = A_B_angle_diff_across_sessions(reg_recall,tuned_log_recall,pf_vector_max_recall,aa,options);
end

%% Filter sessions relative to session 1 based on field shift/trial count - LEARN
if 1
%animal 1, session 3,4 - A/B both trials - means
    %both row and column cross matches
    theta_learn{1}.A(:,4) = {[]};
    theta_learn{1}.A(:,3) = {[]};
    theta_learn{1}.A(4,:) = {[]};
    theta_learn{1}.A(3,:) = {[]};
    
    theta_learn{1}.B(:,4) = {[]};
    theta_learn{1}.B(:,3) = {[]};
    theta_learn{1}.B(4,:) = {[]};
    theta_learn{1}.B(3,:) = {[]};
    
    %animal 3, session 2 - A/B both trials - means
    theta_learn{3}.A(:,2) = {[]};
    theta_learn{3}.A(2,:) = {[]};
    
    theta_learn{3}.B(:,2) = {[]};
    theta_learn{3}.B(2,:) = {[]};
    
    %animal 5, session 6 and 7 - A/B both trials - means
    theta_learn{5}.A(:,6) = {[]};
    theta_learn{5}.A(:,7) = {[]};
    theta_learn{5}.A(6,:) = {[]};
    theta_learn{5}.A(7,:) = {[]};
    
    theta_learn{5}.B(:,6) = {[]};
    theta_learn{5}.B(:,7) = {[]};
    theta_learn{5}.B(6,:) = {[]};
    theta_learn{5}.B(7,:) = {[]};
end

%% Parse sessions by day duration relative to day 1 (since not all sessions were recording as sequential days) 
%arrange the angles to match relative to d1..d9 for recall
[merge_theta_learn_days] = arrange_learning_session_days(theta_learn);

%arrange the angles to match relative to d1..d9 for recall
[merge_theta_recall_days] = arrange_recall_session_days(theta_recall);

%merge across all animals for learning
for dd=1:9
    merge_theta_learn_days_all.A{dd} = cell2mat(merge_theta_learn_days.A(:,dd)');
    merge_theta_learn_days_all.B{dd} = cell2mat(merge_theta_learn_days.B(:,dd)');
end

%merge across all animals for recall
for dd=1:9
    merge_theta_recall_days_all.A{dd} = cell2mat(merge_theta_recall_days.A(:,dd)');
    merge_theta_recall_days_all.B{dd} = cell2mat(merge_theta_recall_days.B(:,dd)');
end

%% Look at distribution of centroid differences in histograms - consistent with bar charts
%Relative to reward zones

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

%% Distance relative to reward zone - old code
%{

%% Distance relative to respective reward zones!!! - T.S. neurons for now (add option for SI)


%% Compare the distance to start of reward zone A and reward zone B for learning/recall
% LEARNING SETS
%ADD SELECTOR FOR TS AND SI TUNED HERE
%for each animal
options.selectTrial = [4,5];
%flag to indicate if this is a learning or recall set
%options.learnSet = 1;
%for each animal
for aa=1:6
    [theta_reward_zone_learn{aa}] = A_B_angle_diff_across_sessions_reward_zones(reg_learn,tuned_log_learning,pf_vector_max_learning,aa,options);
end

% RECALL SETS
options.selectTrial = [1,2];
%flag to indicate if this is a learning or recall set
%options.learnSet = 0;
for aa =1:5
    [theta_reward_zone_recall{aa}] = A_B_angle_diff_across_sessions_reward_zones(reg_recall,tuned_log_recall,pf_vector_max_recall,aa,options);
end


%% Parse session in distance into day distance based on imaging schedule for each animal

%merge animals according to day distance for learning (diff from rew zone
%for each 
%LEARNING MERGE FOR DAY ALIGNMENT
[merge_theta_rew_zones_learn_days] = arrange_learning_session_days_rew_dist(theta_reward_zone_learn);

%RECALL MERGE FOR DAY ALIGNMENT
[merge_theta_rew_zones_recall_days] = arrange_recall_session_rew_dist_days(theta_reward_zone_recall);

%merge across all animals for learning
for dd=1:9
    merge_theta_learn_rew_zone_days_all.A{dd} = cell2mat(merge_theta_rew_zones_learn_days.A.rewA(:,dd));
    merge_theta_learn_rew_zone_days_all.B{dd} = cell2mat(merge_theta_rew_zones_learn_days.B.rewB(:,dd));
end

%merge across all animals for recall
for dd=1:9
    merge_theta_recall_rew_zone_days_all.A{dd} = cell2mat(merge_theta_rew_zones_recall_days.A.rewA(:,dd));
    merge_theta_recall_rew_zone_days_all.B{dd} = cell2mat(merge_theta_rew_zones_recall_days.B.rewB(:,dd));
end

%% Look at distribution of centroid differences in histograms - consistent with bar charts
%Relative to reward zones

subplot_indices = [1:8; 9:16];
%learning
figure
for ii=2:9
    subplot(2,8,subplot_indices(1,ii-1))
hold on
ylim([0 0.7])
xlim([0 pi])
title('Learning - dist to rew zone');
ylabel('Normalized density')
xlabel('Angular diff. [rad]')
    histogram(merge_theta_learn_rew_zone_days_all.A{ii},'Normalization','probability','BinEdges',0:(pi/9):pi)
end
%recall
for ii=2:9
    subplot(2,8,subplot_indices(2,ii-1))
    hold on
    ylim([0 0.7])
    xlim([0 pi])
    title('Recall - dist to rew zone')
    ylabel('Normalized density')
    xlabel('Angular diff. [rad]')
    histogram(merge_theta_recall_rew_zone_days_all.A{ii},'Normalization','probability','BinEdges',0:(pi/9):pi)
end


%% Get mean for learn and recall of centroid difference relative to reward zone

%TODO = check for nan values to prevent skew of sem calculation

%

%%%%%% A %%%%%%%
%get mean
mean_rew_zone.learn.A = cellfun(@nanmean, merge_theta_learn_rew_zone_days_all.A(2:end));
mean_rew_zone.recall.A = cellfun(@nanmean, merge_theta_recall_rew_zone_days_all.A(2:end));

%number of neurons
nb_rew_zone.learn.A = cellfun(@(x) size(x,1), merge_theta_learn_rew_zone_days_all.A);
nb_rew_zone.recall.A = cellfun(@(x) size(x,1), merge_theta_recall_rew_zone_days_all.A);

%sem
sem_rew_zone.learn.A = cellfun(@nanstd, merge_theta_learn_rew_zone_days_all.A(2:end))./sqrt(nb_rew_zone.learn.A(2:end));
sem_rew_zone.recall.A = cellfun(@nanstd,  merge_theta_recall_rew_zone_days_all.A(2:end))./sqrt(nb_rew_zone.recall.A(2:end));

%%%%% B %%%%%
%get mean
mean_rew_zone.learn.B = cellfun(@nanmean, merge_theta_learn_rew_zone_days_all.B(2:end));
mean_rew_zone.recall.B = cellfun(@nanmean, merge_theta_recall_rew_zone_days_all.B(2:end));

%number of neurons
nb_rew_zone.learn.B = cellfun(@(x) size(x,1), merge_theta_learn_rew_zone_days_all.B);
nb_rew_zone.recall.B = cellfun(@(x) size(x,1), merge_theta_recall_rew_zone_days_all.B);

%sem
sem_rew_zone.learn.B = cellfun(@nanstd, merge_theta_learn_rew_zone_days_all.B(2:end))./sqrt(nb_rew_zone.learn.B(2:end));
sem_rew_zone.recall.B = cellfun(@nanstd,  merge_theta_recall_rew_zone_days_all.B(2:end))./sqrt(nb_rew_zone.recall.B(2:end));
%}

%% Plot bar charts for change in distance to A vs. B reward zone for learning vs. recall - REWARD ZONES (SKIP)
if 0
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
end
