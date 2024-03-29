function [outputArg1,outputArg2] = plot_fraction_tuned(tuned_frac_learning,tuned_frac_recall)

%% Get number of session in learning and recall

nb_learn_animals = size(tuned_frac_learning,2);
nb_recall_animals = size(tuned_frac_recall,2);

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
%sem for recall
sem_recall_si = nanstd(frac_recall_all_si,0,3)./sqrt(size(tuned_frac_recall,2));
sem_recall_ts = nanstd(frac_recall_all_ts,0,3)./sqrt(size(tuned_frac_recall,2));


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
title('Recall')
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
sh2 = bar(mean_recall_ts,'stacked');
set(gca,'FontSize',18)
set(gca,'LineWidth',2)

%set individual color bars
for ii = 1:size(sh2,2)
    sh2(ii).FaceColor = 'flat';
    sh2(ii).CData = bar_colorset(ii,:);
end


legend({'A','B','A&B','Neither'})

%% Export T.S. data for figure 4C
bar_colorset = [65,105,225; 220,20,60; 139, 0, 139; 128 128 128]./255;

figure('Position' ,[2091 342 1021 348])
subplot(1,3,1)
hold on
axis square
ylim([0 1.2])
xlabel('Day')
ylabel('Fraction tuned')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:9)
title('Learning - day')
sh1 = bar(mean_learning_ts_day,'stacked');
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
xlabel('Day')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
xticklabels({'1','2','3','6','7','8','9'})
title('Recall - day')
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
xlabel('Day')
yticks([0 0.2 0.4 0.6 0.8 1])
xticks(1:7)
xticklabels({'1','2','3','6','7','8','9'})
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

