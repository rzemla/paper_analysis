function [outputArg1,outputArg2] = PV_TC_correlation_plotter(corr_analysis)

%% Export the structures

st_learn = corr_analysis.st_learn;
st_recall = corr_analysis.st_recall;
lt_recall = corr_analysis.lt_recall;


%% PV and TC line plots for figures for Figure 4 and supplement

%WORK ON THIS



%% Below are control\QC plots 


%% Short term - Plot TC as function of time with sem - T.S - by animal

%short term plots
figure('Position', [2136 350 1031 420])
subplot(1,2,1)
hold on
title('TC - ts - by animal (normalized overlay)')
ylim([0 1])
errorbar(1:7, st_learn.ts.mean_TC.animal.A(1:7),st_learn.ts.sem_TC.animal.A(1:7),'b')
errorbar(1:7, st_learn.ts.mean_TC.animal.B(1:7),st_learn.ts.sem_TC.animal.B(1:7),'r')

errorbar(1:9, st_recall.ts.mean_TC.animal.A(1:9),st_recall.ts.sem_TC.animal.A(1:9),'b')
errorbar(1:9, st_recall.ts.mean_TC.animal.B(1:9),st_recall.ts.sem_TC.animal.B(1:9),'r')
%overlay
errorbar(1:7, st_learn.ts.norm.mean_TC.animal.A(1:7),st_learn.ts.norm.sem_TC.animal.A(1:7),'c--')
errorbar(1:7, st_learn.ts.norm.mean_TC.animal.B(1:7),st_learn.ts.norm.sem_TC.animal.B(1:7),'m--')

errorbar(1:9, st_recall.ts.norm.mean_TC.animal.A(1:9),st_recall.ts.norm.sem_TC.animal.A(1:9),'c--')
errorbar(1:9, st_recall.ts.norm.mean_TC.animal.B(1:9),st_recall.ts.norm.sem_TC.animal.B(1:9),'m--')


%overlay normalized data as a control
subplot(1,2,2)
hold on
ylim([0 1])
title('TC - ts - by animal - correlation on normalized STCs for control (no overlay)')
errorbar(1:7, st_learn.ts.norm.mean_TC.animal.A(1:7),st_learn.ts.norm.sem_TC.animal.A(1:7),'--')
errorbar(1:7, st_learn.ts.norm.mean_TC.animal.B(1:7),st_learn.ts.norm.sem_TC.animal.B(1:7),'--')

errorbar(1:9, st_recall.ts.norm.mean_TC.animal.A(1:9),st_recall.ts.norm.sem_TC.animal.A(1:9),'--')
errorbar(1:9, st_recall.ts.norm.mean_TC.animal.B(1:9),st_recall.ts.norm.sem_TC.animal.B(1:9),'--')


%% Short term - Plot TC as function of time with sem - T.S - by animal (compare with cross day correlations)

%short term plots
figure('Position', [1921 41 1920 963])
subplot(2,3,1)
hold on
axis square
xlabel('Day')
ylabel('Correlation')
title('TC - ts - by animal ')
ylim([0 1])
errorbar(1:8, st_learn.ts.mean_TC.animal.A(1:8),st_learn.ts.sem_TC.animal.A(1:8),'b--')
errorbar(1:8, st_learn.ts.mean_TC.animal.B(1:8),st_learn.ts.sem_TC.animal.B(1:8),'r--')

errorbar(1:9, st_recall.ts.mean_TC.animal.A(1:9),st_recall.ts.sem_TC.animal.A(1:9),'b')
errorbar(1:9, st_recall.ts.mean_TC.animal.B(1:9),st_recall.ts.sem_TC.animal.B(1:9),'r')
legend({'A learn','B learn', 'A recall sub','B recall sub'})

%overlay normalized data as a control
subplot(2,3,2)
hold on
axis square
xlabel('Day')
ylabel('Correlation')
ylim([0 1])
title('TC - ts - by animal - \newline using recall data generated from distance-day correlations')
errorbar(1:8, st_learn.ts.mean_TC.animal.A(1:8),st_learn.ts.sem_TC.animal.A(1:8),'b--')
errorbar(1:8, st_learn.ts.mean_TC.animal.B(1:8),st_learn.ts.sem_TC.animal.B(1:8),'r--')

errorbar(1:9, st_recall.ts.mean_TC.all_corr.A(1:9),st_recall.ts.sem_TC.all_corr.A(1:9),'b-')
errorbar(1:9, st_recall.ts.mean_TC.all_corr.B(1:9),st_recall.ts.sem_TC.all_corr.B(1:9),'r-')
legend({'A learn','B learn', 'A recall all-corr','B recall all-corr'})

subplot(2,3,3)
hold on
axis square
xlabel('Day')
ylabel('Correlation')
ylim([0 1])
title('TC - ts - by animal - \newline overlay D4/D5 substituion and all-day correlation data')
%D4/D5 substitution
errorbar(1:9, st_recall.ts.mean_TC.animal.A(1:9),st_recall.ts.sem_TC.animal.A(1:9),'b')
errorbar(1:9, st_recall.ts.mean_TC.animal.B(1:9),st_recall.ts.sem_TC.animal.B(1:9),'r')

%all day cross correlations
errorbar(1:9, st_recall.ts.mean_TC.all_corr.A(1:9),st_recall.ts.sem_TC.all_corr.A(1:9),'b--')
errorbar(1:9, st_recall.ts.mean_TC.all_corr.B(1:9),st_recall.ts.sem_TC.all_corr.B(1:9),'r--')
legend({'A','B', 'A all-corr','B all-corr'})

% Short term - Plot TC as function of time with sem - S.I. - by animal (compare with cross day correlations)

%short term plots
%figure('Position', [2136 350 1031 420])
subplot(2,3,4)
hold on
axis square
xlabel('Day')
ylabel('Correlation')
title('TC - si - by animal ')
ylim([0 1])
errorbar(1:8, st_learn.si.mean_TC.animal.A(1:8),st_learn.si.sem_TC.animal.A(1:8),'b--')
errorbar(1:8, st_learn.si.mean_TC.animal.B(1:8),st_learn.si.sem_TC.animal.B(1:8),'r--')

errorbar(1:9, st_recall.si.mean_TC.animal.A(1:9),st_recall.si.sem_TC.animal.A(1:9),'b')
errorbar(1:9, st_recall.si.mean_TC.animal.B(1:9),st_recall.si.sem_TC.animal.B(1:9),'r')
legend({'A learn','B learn', 'A recall sub','B recall sub'})

%overlay normalized data as a control
subplot(2,3,5)
hold on
axis square
xlabel('Day')
ylabel('Correlation')
ylim([0 1])
title('TC - si - by animal - \newline using recall data generated from distance-day correlations')
errorbar(1:8, st_learn.si.mean_TC.animal.A(1:8),st_learn.si.sem_TC.animal.A(1:8),'b--')
errorbar(1:8, st_learn.si.mean_TC.animal.B(1:8),st_learn.si.sem_TC.animal.B(1:8),'r--')

errorbar(1:9, st_recall.si.mean_TC.all_corr.A(1:9),st_recall.si.sem_TC.all_corr.A(1:9),'b-')
errorbar(1:9, st_recall.si.mean_TC.all_corr.B(1:9),st_recall.si.sem_TC.all_corr.B(1:9),'r-')
legend({'A learn','B learn', 'A recall all-corr','B recall all-corr'})

subplot(2,3,6)
hold on
axis square
xlabel('Day')
ylabel('Correlation')
ylim([0 1])
title('TC - si - by animal - \newline overlay D4/D5 substituion and all-day correlation data')
%D4/D5 substitution
errorbar(1:9, st_recall.si.mean_TC.animal.A(1:9),st_recall.si.sem_TC.animal.A(1:9),'b')
errorbar(1:9, st_recall.si.mean_TC.animal.B(1:9),st_recall.si.sem_TC.animal.B(1:9),'r')

%all day cross correlations
errorbar(1:9, st_recall.si.mean_TC.all_corr.A(1:9),st_recall.si.sem_TC.all_corr.A(1:9),'b--')
errorbar(1:9, st_recall.si.mean_TC.all_corr.B(1:9),st_recall.si.sem_TC.all_corr.B(1:9),'r--')
legend({'A','B', 'A all-corr','B all-corr'})

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

%% Plot PV as function of time with sem; 
%compare against plot with correlation generated from all day correlation
%for short term recall

%short term plots
figure
subplot(1,2,1)
hold on
ylim([0 1])
%plot(mean_PV.A(1:7),'b--')
title('PV correlation - by animal')
errorbar(1:7, st_learn.mean_PV.A(1:7),st_learn.sem_PV.A(1:7))
errorbar(1:7, st_learn.mean_PV.B(1:7),st_learn.sem_PV.B(1:7))

errorbar(1:9, st_recall.mean_PV.A(1:9),st_recall.sem_PV.A(1:9))
errorbar(1:9, st_recall.mean_PV.B(1:9),st_recall.sem_PV.B(1:9))

subplot(1,2,2)
hold on
ylim([0 1])
%plot(mean_PV.A(1:7),'b--')
title('PV correlation - by animal - using all-day correlations for recall')
errorbar(1:7, st_learn.mean_PV.A(1:7),st_learn.sem_PV.A(1:7))
errorbar(1:7, st_learn.mean_PV.B(1:7),st_learn.sem_PV.B(1:7))

errorbar(1:9, st_recall.mean_PV.all_corr.A(1:9),st_recall.sem_PV.all_corr.A(1:9))
errorbar(1:9, st_recall.mean_PV.all_corr.B(1:9),st_recall.sem_PV.all_corr.B(1:9))


%long term plots
figure
hold on
title('PV correlation - by animal')
ylim([0 1])
plot_days = [1 6 16 20 25 30];
xticks(plot_days)
errorbar(plot_days, lt_recall.mean_PV.A(plot_days),lt_recall.sem_PV.A(plot_days))
errorbar(plot_days, lt_recall.mean_PV.B(plot_days),lt_recall.sem_PV.B(plot_days))

%% Compare PV plot with short term recall generated for substitution of day 4 and day 5 vs. correlation against each point

figure('Position', [2093 485 1601 420])
subplot(1,3,1)
hold on
ylim([0 1])
title('PV - ST recall with D4 & D5 substitution')
xlabel('Day')
ylabel('Correlation')

%short term learning data
errorbar(1:8, st_learn.mean_PV.A(1:8),st_learn.sem_PV.A(1:8),'b--')
errorbar(1:8, st_learn.mean_PV.B(1:8),st_learn.sem_PV.B(1:8),'r--')

%st recall with substitution for D4 and D5
errorbar(1:9, st_recall.mean_PV.A(1:9),st_recall.sem_PV.A(1:9),'b')
errorbar(1:9, st_recall.mean_PV.B(1:9),st_recall.sem_PV.B(1:9),'r')
legend({'A learn','B learn', 'A recall sub','B recall sub'})

subplot(1,3,2)
hold on
ylim([0 1])
title('PV - ST recall with all day correlation distance days')
xlabel('Day')
ylabel('Correlation')

%short term learning data
errorbar(1:8, st_learn.mean_PV.A(1:8),st_learn.sem_PV.A(1:8),'b')
errorbar(1:8, st_learn.mean_PV.B(1:8),st_learn.sem_PV.B(1:8),'r')

%st recall with all cross day correlations
errorbar(1:9, st_recall.mean_PV.all_corr.A(1:9),st_recall.sem_PV.all_corr.A(1:9),'b--')
errorbar(1:9, st_recall.mean_PV.all_corr.B(1:9),st_recall.sem_PV.all_corr.B(1:9),'r--')
legend({'A learn','B learn', 'A recall all-corr','B recall all-corr'})

subplot(1,3,3)
hold on
ylim([0 1])
title('PV - ST recall with D4 & D5 substitution \newline vs. ST recall with all day correlation distance days')
xlabel('Day')
ylabel('Correlation')

%st recall with all cross day correlations
errorbar(1:9, st_recall.mean_PV.all_corr.A(1:9),st_recall.sem_PV.all_corr.A(1:9),'b--')
errorbar(1:9, st_recall.mean_PV.all_corr.B(1:9),st_recall.sem_PV.all_corr.B(1:9),'r--')

%st recall with substitution for D4 and D5
errorbar(1:9, st_recall.mean_PV.A(1:9),st_recall.sem_PV.A(1:9),'b')
errorbar(1:9, st_recall.mean_PV.B(1:9),st_recall.sem_PV.B(1:9),'r')
legend({'A','B', 'A all-corr','B all-corr'})


end

