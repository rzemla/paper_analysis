function [outputArg1,outputArg2] = plot_fraction_tuned(tuned_frac_learning,tuned_frac_recall)


%% Get mean and sem for fractions learned and recall

%combined all animals into 1 matrix

%learning
for aa=1:size(tuned_frac_learning,2)
    frac_learning_all_si(:,:,aa) = tuned_frac_learning{aa}.tuned_fractions.fracTuned_si_filt;
    frac_learning_all_ts(:,:,aa) = tuned_frac_learning{aa}.tuned_fractions.fracTuned_ts_filt;
end

%mean for learning
mean_learning_si = nanmean(frac_learning_all_si,3);
mean_learning_ts = nanmean(frac_learning_all_ts,3);
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
%sem for recall
sem_recall_si = nanstd(frac_recall_all_si,0,3)./sqrt(size(tuned_frac_recall,2));
sem_recall_ts = nanstd(frac_recall_all_ts,0,3)./sqrt(size(tuned_frac_recall,2));


%% Plot as bar chart for each session

figure('Position' ,[1964 643 1063 225])
subplot(1,2,1)
hold on
ylim([0 1])
title('Learning - S.I.')
bar(mean_learning_si,'stacked')
subplot(1,2,2)
hold on
ylim([0 1])
title('Recall - S.I.')
bar(mean_recall_si,'stacked')

%% Export T.S. data for figure 4C
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
title('Learning')
sh1 = bar(mean_learning_ts,'stacked');
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
title('Recall')
sh2 = bar(mean_recall_ts,'stacked');
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
sh2 = bar(mean_recall_ts,'stacked');
set(gca,'FontSize',18)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh2,2)
    sh2(ii).FaceColor = 'flat';
    sh2(ii).CData = bar_colorset(ii,:);
end


legend({'A','B','A&B','Neither'})

end

