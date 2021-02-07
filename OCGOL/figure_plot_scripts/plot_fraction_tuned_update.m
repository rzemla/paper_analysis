function [frac_tuned_across_days] = plot_fraction_tuned_update(tuned_frac_learning,tuned_frac_recall,tuned_frac_long_recall)

%% Get number of session in learning and recall

nb_learn_animals = size(tuned_frac_learning,2);
nb_recall_animals = size(tuned_frac_recall,2);
nb_long_recall_animals = size(tuned_frac_long_recall,2);

for aa=1:nb_learn_animals
    nb_ses_learn(aa) = size(tuned_frac_learning{aa}.tuned_fractions.tuned_count_filt,2);
end

for aa=1:nb_recall_animals
    nb_ses_recall(aa) = size(tuned_frac_recall{aa}.tuned_fractions.tuned_count_filt,2);
end



%% Get mean and sem for fractions learned and recall - combine fractions into one matrix

%preallocation with nans given variable number of sessions for each animal
%1st dim - row - session
%middle dim - column is the relative A,B,A&B, neither fractional count
%3rd dim - animal index
frac_learning_all_si = nan(max(nb_ses_learn),4,nb_learn_animals);
frac_learning_all_ts = nan(max(nb_ses_learn),4,nb_learn_animals);

%learning
for aa=1:nb_learn_animals
    frac_learning_all_si(1:nb_ses_learn(aa),:,aa) = tuned_frac_learning{aa}.tuned_fractions.fracTuned_si_filt;
    frac_learning_all_ts(1:nb_ses_learn(aa),:,aa) = tuned_frac_learning{aa}.tuned_fractions.fracTuned_ts_filt;
end

%% Filter out sessions with low trial count/field shift (conservative)

%animal 1, session 3,4 - A/B both trials - means
frac_learning_all_si([3,4],:,1) = nan;
frac_learning_all_ts([3,4],:,1) = nan;
%animal 3, session 2 - A/B both trials - means
frac_learning_all_si([2],:,3) = nan;
frac_learning_all_ts([2],:,3) = nan;
%animal 5, session 6 and 7 - A/B both trials - means
frac_learning_all_si([6,7],:,5) = nan;
frac_learning_all_ts([6,7],:,5) = nan;

%% Rearrange learning trials into day order - LEARN
%preallocate with nan
frac_learning_all_si_day = nan(max(nb_ses_learn),4,nb_learn_animals);
frac_learning_all_ts_day = nan(max(nb_ses_learn),4,nb_learn_animals);

for dd=1:9
    switch dd
        case 1
            for aa=1:6
                frac_learning_all_si_day(dd,:,aa) = frac_learning_all_si(dd,:,aa);
                frac_learning_all_ts_day(dd,:,aa) = frac_learning_all_ts(dd,:,aa);
            end            
        case 2
            for aa=1:6
                frac_learning_all_si_day(dd,:,aa) = frac_learning_all_si(dd,:,aa);
                frac_learning_all_ts_day(dd,:,aa) = frac_learning_all_ts(dd,:,aa);
            end
        case 3
            for aa=1:6
                frac_learning_all_si_day(dd,:,aa) = frac_learning_all_si(dd,:,aa);
                frac_learning_all_ts_day(dd,:,aa) = frac_learning_all_ts(dd,:,aa);
            end
        case 4
            for aa=1:6
                frac_learning_all_si_day(dd,:,aa) = frac_learning_all_si(dd,:,aa);
                frac_learning_all_ts_day(dd,:,aa) = frac_learning_all_ts(dd,:,aa);
            end
        case 5
            for aa=1:6
                if aa==1 || aa==2
                    frac_learning_all_si_day(dd,:,aa) = nan(1,4);
                    frac_learning_all_ts_day(dd,:,aa) = nan(1,4);
                elseif (aa>= 3 && aa<= 6)
                    frac_learning_all_si_day(dd,:,aa) = frac_learning_all_si(5,:,aa);
                    frac_learning_all_ts_day(dd,:,aa) = frac_learning_all_ts(5,:,aa);
                end
            end
        case 6
            for aa=1:6
                if aa==1 || aa==2
                    frac_learning_all_si_day(dd,:,aa) = frac_learning_all_si(5,:,aa);
                    frac_learning_all_ts_day(dd,:,aa) = frac_learning_all_ts(5,:,aa);
                elseif (aa>= 3 && aa<= 6)
                    frac_learning_all_si_day(dd,:,aa) = frac_learning_all_si(6,:,aa);
                    frac_learning_all_ts_day(dd,:,aa) = frac_learning_all_ts(6,:,aa);
                end
            end            
        case 7
            for aa=1:6
                if aa==1 || aa==2
                    frac_learning_all_si_day(dd,:,aa) = frac_learning_all_si(6,:,aa);
                    frac_learning_all_ts_day(dd,:,aa) = frac_learning_all_ts(6,:,aa);
                elseif aa==3 || aa==4
                    frac_learning_all_si_day(dd,:,aa) = frac_learning_all_si(7,:,aa);
                    frac_learning_all_ts_day(dd,:,aa) = frac_learning_all_ts(7,:,aa);
                elseif aa==5 || aa==6
                    frac_learning_all_si_day(dd,:,aa) = nan(1,4);
                    frac_learning_all_ts_day(dd,:,aa) = nan(1,4);
                end
            end
        case 8
            for aa=1:6
                if aa==1 || aa==2 || aa==5
                    frac_learning_all_si_day(dd,:,aa) = frac_learning_all_si(7,:,aa);
                    frac_learning_all_ts_day(dd,:,aa) = frac_learning_all_ts(7,:,aa);
                elseif aa==3 || aa==4
                    frac_learning_all_si_day(dd,:,aa) = frac_learning_all_si(8,:,aa);
                    frac_learning_all_ts_day(dd,:,aa) = frac_learning_all_ts(8,:,aa);
                elseif aa==6
                    frac_learning_all_si_day(dd,:,aa) = nan(1,4);
                    frac_learning_all_ts_day(dd,:,aa) = nan(1,4);
                end
            end 
        case 9
            for aa=1:6
                if aa==1 || aa==2 || aa==3 || aa==5 || aa==6
                    frac_learning_all_si_day(dd,:,aa) = nan(1,4);
                    frac_learning_all_ts_day(dd,:,aa) = nan(1,4);
                elseif aa==4
                    frac_learning_all_si_day(dd,:,aa) = frac_learning_all_si(9,:,aa);
                    frac_learning_all_ts_day(dd,:,aa) = frac_learning_all_ts(9,:,aa);
                end
            end            
    end
end


%% Export data to Prism for analysis only extract data until day 7 for learning

%1st dim - day
%2nd dim - A,B,A&B,neither, 
%3rd dim - animal

%SI
%output: animal x day
si_frac_export.learning.A = squeeze(frac_learning_all_si_day(1:7,1,:))';
si_frac_export.learning.B = squeeze(frac_learning_all_si_day(1:7,2,:))';
si_frac_export.learning.AB = squeeze(frac_learning_all_si_day(1:7,3,:))';

%TS
ts_frac_export.learning.A = squeeze(frac_learning_all_ts_day(1:7,1,:))';
ts_frac_export.learning.B = squeeze(frac_learning_all_ts_day(1:7,2,:))';
ts_frac_export.learning.AB = squeeze(frac_learning_all_ts_day(1:7,3,:))';


%% Export data to Prism for analysis only extract data until day 7 for short term recall

%same structure - no missing values - 
%order is d1, d2, d3, d6,d7,d8,d9 (linear following session)

%1st dim - day
%2nd dim - A,B,A&B,neither, 
%3rd dim - animal

%recall
for aa=1:size(tuned_frac_recall,2)
    frac_recall_all_si(:,:,aa) = tuned_frac_recall{aa}.tuned_fractions.fracTuned_si_filt;
    frac_recall_all_ts(:,:,aa) = tuned_frac_recall{aa}.tuned_fractions.fracTuned_ts_filt;
end

%frac_recall_all_session_si = frac_recall_all_si 
%ordered by session/day (1 2 3 6 7 8 9)

%SI
%output: animal x day
si_frac_export.recall.A = squeeze(frac_recall_all_si(1:7,1,:))';
si_frac_export.recall.B = squeeze(frac_recall_all_si(1:7,2,:))';
si_frac_export.recall.AB = squeeze(frac_recall_all_si(1:7,3,:))';

%TS 
ts_frac_export.recall.A = squeeze(frac_recall_all_ts(1:7,1,:))';
ts_frac_export.recall.B = squeeze(frac_recall_all_ts(1:7,2,:))';
ts_frac_export.recall.AB = squeeze(frac_recall_all_ts(1:7,3,:))';


%% Import long recall data

%session, trial type, animal
for aa=1:size(tuned_frac_long_recall,2)
    frac_long_recall_all_si(:,:,aa) = tuned_frac_long_recall{aa}.tuned_fractions.fracTuned_si_filt;
    frac_long_recall_all_ts(:,:,aa) = tuned_frac_long_recall{aa}.tuned_fractions.fracTuned_ts_filt;
end

%% Exclude sessions for long term recall data 
%Exclude sessions with less than or equal to 5 laps in A or B
%animal 3 , both A/B due to field shift, session 2 3 4 

frac_long_recall_all_si([2 3 4],:,3) = nan;
frac_long_recall_all_ts([2 3 4],:,3) = nan;

%% Export long term recall data to Prism

%SI
%output: animal x day
si_frac_export.long_recall.A = squeeze(frac_long_recall_all_si(1:6,1,:))';
si_frac_export.long_recall.B = squeeze(frac_long_recall_all_si(1:6,2,:))';
si_frac_export.long_recall.AB = squeeze(frac_long_recall_all_si(1:6,3,:))';

%TS 
ts_frac_export.long_recall.A = squeeze(frac_long_recall_all_ts(1:6,1,:))';
ts_frac_export.long_recall.B = squeeze(frac_long_recall_all_ts(1:6,2,:))';
ts_frac_export.long_recall.AB = squeeze(frac_long_recall_all_ts(1:6,3,:))';


%% Get mean and SEM for learning animal (SEM not used in plot)
%modify to reflect # of animals
%mean for learning
mean_learning_si = nanmean(frac_learning_all_si,3);
mean_learning_ts = nanmean(frac_learning_all_ts,3);

%for relative days
mean_learning_si_day = nanmean(frac_learning_all_si_day,3);
mean_learning_ts_day = nanmean(frac_learning_all_ts_day,3);

%sem for learning
sem_learning_si = nanstd(frac_learning_all_si,0,3)./sqrt(size(tuned_frac_learning,2));
sem_learning_ts = nanstd(frac_learning_all_ts,0,3)./sqrt(size(tuned_frac_learning,2));

%recall
for aa=1:size(tuned_frac_recall,2)
    frac_recall_all_si(:,:,aa) = tuned_frac_recall{aa}.tuned_fractions.fracTuned_si_filt;
    frac_recall_all_ts(:,:,aa) = tuned_frac_recall{aa}.tuned_fractions.fracTuned_ts_filt;
end

%mean for recall
mean_recall_si = nanmean(frac_recall_all_si,3);
mean_recall_ts = nanmean(frac_recall_all_ts,3);

%mean for long recall
mean_long_recall_si = nanmean(frac_long_recall_all_si,3);
mean_long_recall_ts = nanmean(frac_long_recall_all_ts,3);

%sem for recall
sem_recall_si = nanstd(frac_recall_all_si,0,3)./sqrt(size(tuned_frac_recall,2));
sem_recall_ts = nanstd(frac_recall_all_ts,0,3)./sqrt(size(tuned_frac_recall,2));


%% Export T.S. data for figure 4C -SI
bar_colorset = [65,105,225; 220,20,60; 139, 0, 139; 128 128 128]./255;

figure('Position' ,[2016 473 1016 307])
subplot(1,3,1)
hold on
axis square
ylim([0 1.2])
xlabel('Sessions')
ylabel('Fraction tuned')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:6)
title('Learning SI')
sh1 = bar(mean_learning_si,'stacked');
set(gca,'FontSize',18)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh1,2)
    sh1(ii).FaceColor = 'flat';
    sh1(ii).CData = bar_colorset(ii,:);
end

subplot(1,3,2)
hold on
axis square
ylim([0 1.2])
xlabel('Sessions')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
title('Recall SI')
sh2 = bar(mean_recall_si,'stacked');
set(gca,'FontSize',18)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh2,2)
    sh2(ii).FaceColor = 'flat';
    sh2(ii).CData = bar_colorset(ii,:);
end

%third subplot just to get the bar legend
subplot(1,3,3)
hold on
axis square
ylim([0 1.2])
xlabel('Sessions')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
title('Get legend -delete in illustrator')
sh2 = bar(mean_recall_si,'stacked');
set(gca,'FontSize',18)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh2,2)
    sh2(ii).FaceColor = 'flat';
    sh2(ii).CData = bar_colorset(ii,:);
end


legend({'A','B','A&B','Neither'})


%% Export T.S. data for figure 4C (Day)
bar_colorset = [65,105,225; 220,20,60; 139, 0, 139; 128 128 128]./255;

figure('Position' ,[2091 342 1021 348])
subplot(1,3,1)
hold on
axis square
ylim([0 1.2])
xlabel('Day')
ylabel('Fraction tuned')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
title('Learning - day - T.S.')
sh1 = bar(1:7,mean_learning_ts_day(1:7,:),'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh1,2)
    sh1(ii).FaceColor = 'flat';
    sh1(ii).CData = bar_colorset(ii,:);
end

subplot(1,3,2)
hold on
axis square
ylim([0 1.2])
xlabel('Day')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
xticklabels({'1','2','3','6','7','8','9'})
title('Recall - day - T.S.')
sh2 = bar(mean_recall_ts,'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh2,2)
    sh2(ii).FaceColor = 'flat';
    sh2(ii).CData = bar_colorset(ii,:);
end

%third subplot just to get the bar legend
subplot(1,3,3)
hold on
axis square
ylim([0 1.2])
xlabel('Day')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
xticklabels({'1','2','3','6','7','8','9'})
title('Get legend -delete in illustrator')
sh2 = bar(mean_recall_ts,'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh2,2)
    sh2(ii).FaceColor = 'flat';
    sh2(ii).CData = bar_colorset(ii,:);
end


legend({'A','B','A&B','Neither'})

%% Plot S.I. data for figure 4C (Day)
bar_colorset = [65,105,225; 220,20,60; 139, 0, 139; 128 128 128]./255;

figure('Position' ,[2091 342 1021 348])
subplot(1,3,1)
hold on
axis square
ylim([0 1.2])
xlabel('Day')
ylabel('Fraction tuned')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
title('Learning - day - S.I.')
sh1 = bar(1:7,mean_learning_si_day(1:7,:),'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh1,2)
    sh1(ii).FaceColor = 'flat';
    sh1(ii).CData = bar_colorset(ii,:);
end

subplot(1,3,2)
hold on
axis square
ylim([0 1.2])
xlabel('Day')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
xticklabels({'1','2','3','6','7','8','9'})
title('Recall - day - S.I.')
sh2 = bar(mean_recall_si,'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh2,2)
    sh2(ii).FaceColor = 'flat';
    sh2(ii).CData = bar_colorset(ii,:);
end

%third subplot just to get the bar legend
subplot(1,3,3)
hold on
axis square
ylim([0 1.2])
xlabel('Day')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
xticklabels({'1','2','3','6','7','8','9'})
title('Get legend -delete in illustrator')
sh2 = bar(mean_recall_si,'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh2,2)
    sh2(ii).FaceColor = 'flat';
    sh2(ii).CData = bar_colorset(ii,:);
end


legend({'A','B','A&B','Neither'})

%% Plot Long recall S.I. and T.S data for supplement/thesis
bar_colorset = [65,105,225; 220,20,60; 139, 0, 139; 128 128 128]./255;

figure('Position' ,[2091 342 1021 348])
subplot(1,3,1)
hold on
axis square
ylim([0 1.2])
xlabel('Day')
ylabel('Fraction tuned')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:6)
xticklabels({'1', '6', '16', '20', '25', '30'})
title('Long Recall - day - S.I.')
sh1 = bar(1:6,mean_long_recall_si(1:6,:),'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh1,2)
    sh1(ii).FaceColor = 'flat';
    sh1(ii).CData = bar_colorset(ii,:);
end

subplot(1,3,2)
hold on
axis square
ylim([0 1.2])
xlabel('Day')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:6)
xticklabels({'1', '6', '16', '20', '25', '30'})
title('Long Recall - day - T.S.')
sh2 = bar(1:6,mean_long_recall_ts(1:6,:),'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh2,2)
    sh2(ii).FaceColor = 'flat';
    sh2(ii).CData = bar_colorset(ii,:);
end

%third subplot just to get the bar legend
subplot(1,3,3)
hold on
axis square
title('Irrelevant - space filler')
ylim([0 1.2])
xlabel('Day')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
xticklabels({'1','2','3','6','7','8','9'})
title('Get legend -delete in illustrator')
sh2 = bar(1:6,mean_long_recall_ts(1:6,:),'stacked');
set(gca,'FontSize',14)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh2,2)
    sh2(ii).FaceColor = 'flat';
    sh2(ii).CData = bar_colorset(ii,:);
end


legend({'A','B','A&B','Neither'})

%% Export data for plotting

frac_tuned_across_days.mean_learning_si = mean_learning_si;
frac_tuned_across_days.mean_recall_si = mean_recall_si;
frac_tuned_across_days.mean_learning_ts_day = mean_learning_ts_day;
frac_tuned_across_days.mean_recall_ts = mean_recall_ts;
frac_tuned_across_days.mean_learning_si_day = mean_learning_si_day;
frac_tuned_across_days.mean_recall_si = mean_recall_si;
frac_tuned_across_days.mean_long_recall_si = mean_long_recall_si;
frac_tuned_across_days.mean_long_recall_ts = mean_long_recall_ts;


end

