function [outputArg1,outputArg2] = PV_TC_correlation_analysis(short_term_learn, short_term_recall, long_term_recall,excl_day_combined_day_nan)



%% PV correlation - Filter out low quality sessions and convert to days (learning,ST and LT recall)
%exp_type:
%(1) - short term learn
%(2) - short term recall
%(3) - long term recall

%day_range - vector
%short term learn and short term recall: [1:9];
%long term recall: [1 6 16 20 25 30];

%short term learn
exp_type =1;
day_range = 1:9;

[st_learn.mean_PV,st_learn.sem_PV] = filter_convert_day_return_mean_sem_PV(short_term_learn,excl_day_combined_day_nan,exp_type,day_range);

%short term recall (fills in day 4 and day 5)
exp_type =2;
day_range = 1:9;

[st_recall.mean_PV,st_recall.sem_PV] = filter_convert_day_return_mean_sem_PV(short_term_recall,excl_day_combined_day_nan,exp_type,day_range);

%long term recall
exp_type =3;
day_range = [1 6 16 20 25 30];

[lt_recall.mean_PV,lt_recall.sem_PV] = filter_convert_day_return_mean_sem_PV(long_term_recall,excl_day_combined_day_nan,exp_type,day_range);



%% Plot PV as function of time with sem

%short term plots
figure
hold on
ylim([0 1])
%plot(mean_PV.A(1:7),'b--')
title('PV correlation - all neurons')
errorbar(1:7, st_learn.mean_PV.A(1:7),st_learn.sem_PV.A(1:7))
errorbar(1:7, st_learn.mean_PV.B(1:7),st_learn.sem_PV.B(1:7))

errorbar(1:9, st_recall.mean_PV.A(1:9),st_recall.sem_PV.A(1:9))
errorbar(1:9, st_recall.mean_PV.B(1:9),st_recall.sem_PV.B(1:9))

%long term plots
figure
hold on
title('PV correlation - all neurons')
ylim([0 1])
plot_days = [1 6 16 20 25 30];
xticks(plot_days)
errorbar(plot_days, lt_recall.mean_PV.A(plot_days),lt_recall.sem_PV.A(plot_days))
errorbar(plot_days, lt_recall.mean_PV.B(plot_days),lt_recall.sem_PV.B(plot_days))

%% TC correlation - TS

%short term learn
exp_type =1;
day_range = 1:9;
tuning_type = 'ts';

[st_learn.ts.mean_TC,st_learn.ts.sem_TC] = filter_convert_day_return_mean_sem_TC_global(short_term_learn,excl_day_combined_day_nan,exp_type,day_range,tuning_type);

%short term recall (fills in day 4 and day 5)
exp_type =2;
day_range = 1:9;
tuning_type = 'ts';

[st_recall.ts.mean_TC,st_recall.ts.sem_TC] = filter_convert_day_return_mean_sem_TC_global(short_term_recall,excl_day_combined_day_nan,exp_type,day_range,tuning_type);


%long term recall
exp_type =3;
day_range = [1 6 16 20 25 30];
tuning_type = 'ts';
[lt_recall.ts.mean_TC,lt_recall.ts.sem_TC] = filter_convert_day_return_mean_sem_TC_global(long_term_recall,excl_day_combined_day_nan,exp_type,day_range,tuning_type);


%% Run TC correlation analysis on NORMALIZED spatial tuning curves (control for non-norm)
%modify filter_convert_day_return _mean_sem_TC_global
exp_type =1;
day_range = 1:9;
tuning_type = 'ts';

[st_learn.ts.norm.mean_TC,st_learn.ts.norm.sem_TC] = filter_convert_day_return_mean_sem_TC_global_normalized(short_term_learn,excl_day_combined_day_nan,exp_type,day_range,tuning_type);


%% TC correlation - SI

%short term learn
exp_type =1;
day_range = 1:9;
tuning_type = 'si';

[st_learn.si.mean_TC,st_learn.si.sem_TC] = filter_convert_day_return_mean_sem_TC_global(short_term_learn,excl_day_combined_day_nan,exp_type,day_range,tuning_type);

%short term recall (fills in day 4 and day 5)
exp_type =2;
day_range = 1:9;
tuning_type = 'si';

[st_recall.si.mean_TC,st_recall.si.sem_TC] = filter_convert_day_return_mean_sem_TC_global(short_term_recall,excl_day_combined_day_nan,exp_type,day_range,tuning_type);

%long term recall
exp_type =3;
day_range = [1 6 16 20 25 30];
tuning_type = 'si';
[lt_recall.si.mean_TC,lt_recall.si.sem_TC] = filter_convert_day_return_mean_sem_TC_global(long_term_recall,excl_day_combined_day_nan,exp_type,day_range,tuning_type);



%% Short term - Plot TC as function of time with sem - T.S - by animal

%short term plots
figure('Position', [2136 350 1031 420])
subplot(1,2,1)
hold on
title('TC - ts - by animal')
ylim([0 1])
errorbar(1:7, st_learn.ts.mean_TC.animal.A(1:7),st_learn.ts.sem_TC.animal.A(1:7))
errorbar(1:7, st_learn.ts.mean_TC.animal.B(1:7),st_learn.ts.sem_TC.animal.B(1:7))

errorbar(1:9, st_recall.ts.mean_TC.animal.A(1:9),st_recall.ts.sem_TC.animal.A(1:9))
errorbar(1:9, st_recall.ts.mean_TC.animal.B(1:9),st_recall.ts.sem_TC.animal.B(1:9))

%% Short term - Plot TC as function of time with sem - T.S - by neuron

%short term plots
subplot(1,2,2)
hold on
title('TC - ts - by neuron')
ylim([0 1])
errorbar(1:7, st_learn.ts.mean_TC.neuron.A(1:7),st_learn.ts.sem_TC.neuron.A(1:7))
errorbar(1:7, st_learn.ts.mean_TC.neuron.B(1:7),st_learn.ts.sem_TC.neuron.B(1:7))

errorbar(1:9, st_recall.ts.mean_TC.neuron.A(1:9),st_recall.ts.sem_TC.neuron.A(1:9))
errorbar(1:9, st_recall.ts.mean_TC.neuron.B(1:9),st_recall.ts.sem_TC.neuron.B(1:9))

%% Long term recall TC plots - TS - by animal

figure('Position', [2136 350 1031 420])
subplot(1,2,1)
hold on
title('TC correlation -  LT - TS - by animal')
ylim([0 1])
plot_days = [1 6 16 20 25 30];
xticks(plot_days)
errorbar(plot_days, lt_recall.ts.mean_TC.animal.A(plot_days),lt_recall.ts.sem_TC.animal.A(plot_days))
errorbar(plot_days, lt_recall.ts.mean_TC.animal.B(plot_days),lt_recall.ts.sem_TC.animal.B(plot_days))


%% Long term recall TC plots - TS - by neuron

subplot(1,2,2)
hold on
title('TC correlation -  LT - TS - by neuron')
ylim([0 1])
plot_days = [1 6 16 20 25 30];
xticks(plot_days)
errorbar(plot_days, lt_recall.ts.mean_TC.neuron.A(plot_days),lt_recall.ts.sem_TC.neuron.A(plot_days))
errorbar(plot_days, lt_recall.ts.mean_TC.neuron.B(plot_days),lt_recall.ts.sem_TC.neuron.B(plot_days))


%% Plot TC as function of time with sem - S.I. - by animal

figure('Position', [2136 350 1031 420])
subplot(1,2,1)
hold on
title('TC - si - by animal')
ylim([0 1])
errorbar(1:7, st_learn.si.mean_TC.animal.A(1:7),st_learn.si.sem_TC.animal.A(1:7))
errorbar(1:7, st_learn.si.mean_TC.animal.B(1:7),st_learn.si.sem_TC.animal.B(1:7))

errorbar(1:9, st_recall.si.mean_TC.animal.A(1:9),st_recall.si.sem_TC.animal.A(1:9))
errorbar(1:9, st_recall.si.mean_TC.animal.B(1:9),st_recall.si.sem_TC.animal.B(1:9))

%% Plot TC as function of time with sem - S.I. - by neuron

%short term plots
subplot(1,2,2)
hold on
title('TC - si - by neuron')
ylim([0 1])
errorbar(1:7, st_learn.si.mean_TC.neuron.A(1:7),st_learn.si.sem_TC.neuron.A(1:7))
errorbar(1:7, st_learn.si.mean_TC.neuron.B(1:7),st_learn.si.sem_TC.neuron.B(1:7))

errorbar(1:9, st_recall.si.mean_TC.neuron.A(1:9),st_recall.si.sem_TC.neuron.A(1:9))
errorbar(1:9, st_recall.si.mean_TC.neuron.B(1:9),st_recall.si.sem_TC.neuron.B(1:9))

%% DEVELOPMENT CODE USED TO GENERATE FUNCTION ABOVE

%% Number of animals in each category

% nb_st_learn = size(short_term_learn.PV_TC_corr,2);
% nb_st_recall = size(short_term_recall.PV_TC_corr,2);
% nb_lt_recall = size(long_term_recall.PV_TC_corr,2);


%% Load in non-normalized PV correlation matrices for each class of experiments

% for aa=1:nb_st_learn
%     %short term learn
%     st_learn.PVmat.A = short_term_learn.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.A;
%     st_learn.PVmat.B = short_term_learn.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.B;
% end
% 
% for aa=1:nb_st_recall
%     %short term recall
%     st_recall.PVmat.A = short_term_recall.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.A;
%     st_recall.PVmat.B = short_term_recall.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.B;
% end
% 
% for aa=1:nb_lt_recall
%     %long term recall
%     lt_recall.PVmat.A = long_term_recall.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.A;
%     lt_recall.PVmat.B = long_term_recall.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.B;    
% end



%relative to day 1
%for each animal
% for aa=1:nb_st_learn
%     %for each day
%     for dd=1:9
%         %look for day and get session index that corresponds to that day
%         ses2day_idx = find(excl_day_combined_day_nan{aa,1}(2,:) == dd);
%         
%         %if no assignmenet for that day
%         if ~isempty(ses2day_idx)
%             %animal x day PV matrix cell
%             %for A trials
%             day_PV_mat.st_learn.A{aa,dd} = short_term_learn.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.A{1,ses2day_idx};
%             %for B trials
%             day_PV_mat.st_learn.B{aa,dd} = short_term_learn.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.B{1,ses2day_idx};
%         else
%             day_PV_mat.st_learn.A{aa,dd} = [];
%             day_PV_mat.st_learn.B{aa,dd} = [];
%         end
%         
%     end
%     
% end
% 
% %% Get diagnonal and mean of PV for each day
% 
% %get diagnonal for A
% day_PV_diag.st_learn.A = cellfun(@diag,day_PV_mat.st_learn.A,'UniformOutput',false);
% %get mean of diagnonal for A
% day_PV_diag_mean.st_learn.A = cell2mat(cellfun(@nanmean,day_PV_diag.st_learn.A,'UniformOutput',false));
% %get mean across animals of the mean of the diagnonals
% day_PV_diag_mean_mean_st_learn.A = nanmean(day_PV_diag_mean.st_learn.A,1);

end

