%% import prism fraction of licks data on A and B trials

%% import data for figure 4 stats

dirpath = 'C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data';


%processed data

load(fullfile(dirpath,'source_data_fig4_5_and_sup.mat'))

%load(fullfile(dirpath,''))
%load(fullfile(dirpath,''))
%% long term recall

%% Print mean and sem data to word
% Start Word document that will contain the formatted stats data
WordFileName='long_term_neuron_sup_data_mean_sem.doc';
CurDir=pwd;
FileSpec = fullfile(CurDir,WordFileName);
%active X handle for manipulating document (ActXWord)
[ActXWord,WordHandle]=StartWord(FileSpec);

fprintf('Document will be saved in %s\n',FileSpec);

txt_input = 'PV corr 6 vs 30 recall A trials';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.mean_PV.A(6);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.sem_PV.A(6);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.mean_PV.A(30);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.sem_PV.A(30);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'PV corr 6 vs 30 recall B trials';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.mean_PV.B(6);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.sem_PV.B(6);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.mean_PV.B(30);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.sem_PV.B(30);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'TC SI corr 6 vs 30 recall A trials';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.si.mean_TC.animal.A(6);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.si.sem_TC.animal.A(6);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.si.mean_TC.animal.A(30);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.si.sem_TC.animal.A(30);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'TC SI corr 6 vs 30 recall B trials';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.si.mean_TC.animal.B(6);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.si.sem_TC.animal.B(6);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.si.mean_TC.animal.B(30);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.si.sem_TC.animal.B(30);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);


txt_input = 'TC ts corr 6 vs 30 recall A trials';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.ts.mean_TC.animal.A(6);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.ts.sem_TC.animal.A(6);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.ts.mean_TC.animal.A(30);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.ts.sem_TC.animal.A(30);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'TC ts corr 6 vs 30 recall B trials';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.ts.mean_TC.animal.B(6);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.ts.sem_TC.animal.B(6);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.ts.mean_TC.animal.B(30);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.ts.sem_TC.animal.B(30);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%% Neighboring day correlation long term
%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Neighboring PV corr 1 vs 6 vs 20 vs 25 recall A trials';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.mean_PV.A(1);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.sem_PV.A(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.mean_PV.A(4);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.sem_PV.A(4);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Neighboring PV corr 1 vs 6 vs 20 vs 25 recall B trials';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.mean_PV.B(1);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.sem_PV.B(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.mean_PV.B(4);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.sem_PV.B(4);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Neighboring TC SI corr 1 vs 6 vs 20 vs 25 recall A trials';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.si.mean_TC.animal.A(1);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.si.sem_TC.animal.A(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.si.mean_TC.animal.A(4);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.si.sem_TC.animal.A(4);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Neighboring TC SI corr 1 vs 6 vs 20 vs 25 recall B trials';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.si.mean_TC.animal.B(1);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.si.sem_TC.animal.B(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.si.mean_TC.animal.B(4);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.si.sem_TC.animal.B(4);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);


txt_input = 'Neighboring TC TS corr 1 vs 6 vs 20 vs 25 recall A trials';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.ts.mean_TC.animal.A(1);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.ts.sem_TC.animal.A(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.ts.mean_TC.animal.A(4);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.ts.sem_TC.animal.A(4);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Neighboring TC TS corr 1 vs 6 vs 20 vs 25 recall B trials';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.ts.mean_TC.animal.B(1);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.ts.sem_TC.animal.B(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.ts.mean_TC.animal.B(4);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.ts.sem_TC.animal.B(4);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);


%% Normalized data comparison long term

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Norm A vs. B corr score - SI Day 6';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.AB_corr.si.mean_TC.pooled.AB_corr_ratio_median(2);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.AB_corr.si.sem_TC.pooled.AB_corr_ratio_95ci(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Norm A vs. B corr score - SI Day 25';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.AB_corr.si.mean_TC.pooled.AB_corr_ratio_median(5);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.AB_corr.si.sem_TC.pooled.AB_corr_ratio_95ci(5);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Norm A vs. B corr score - TS Day 6';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.AB_corr.ts.mean_TC.pooled.AB_corr_ratio_median(2);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.AB_corr.ts.sem_TC.pooled.AB_corr_ratio_95ci(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Norm A vs. B corr score - TS Day 25';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_short_learn_recall.corr_analysis.lt_recall.AB_corr.ts.mean_TC.pooled.AB_corr_ratio_median(5);
sem_txt = source_data_short_learn_recall.corr_analysis.lt_recall.AB_corr.ts.sem_TC.pooled.AB_corr_ratio_95ci(5);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);


