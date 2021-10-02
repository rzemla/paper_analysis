%% import prism fraction of licks data on A and B trials

%% import data for figure 4 stats

dirpath = 'C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data';


%processed data
load(fullfile(dirpath,'source_data_ex10_sup.mat'))
load(fullfile(dirpath,'fig_1_speed_rew_zone.mat'))
load(fullfile(dirpath,'source_data_sup_2_3.mat'))

%load(fullfile(dirpath,''))
%load(fullfile(dirpath,''))
%% Sup 8 Common vs activity remapping index score and AUC/min

%% Print mean and sem data to word
% Start Word document that will contain the formatted stats data
WordFileName='activity_vs_remap_speed_PF_neuron_sup_data_mean_sem.doc';
CurDir=pwd;
FileSpec = fullfile(CurDir,WordFileName);
%active X handle for manipulating document (ActXWord)
[ActXWord,WordHandle]=StartWord(FileSpec);

fprintf('Document will be saved in %s\n',FileSpec);

txt_input = 'Common vs activity remapping index scores';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.remap_sup_plot_data.common_idx_values;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in =source_data_sup_2_3.remap_sup_plot_data.remap_idx_values;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Activity rate in RUN on A vs. B remapping neurons';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = source_data_sup_2_3.remap_sup_plot_data.activity_remap_mean_AUC_min(1);
sem_txt = source_data_sup_2_3.remap_sup_plot_data.activity_remap_sem_AUC_min(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = source_data_sup_2_3.remap_sup_plot_data.activity_remap_mean_AUC_min(2);
sem_txt = source_data_sup_2_3.remap_sup_plot_data.activity_remap_sem_AUC_min(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%% Partial field centroid dist
txt_input = 'Partial remapping place field centroids';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_AB_remap.partial_shift.partial_A_common_far_input(:,2)*1.96;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in =source_data_AB_remap.partial_shift.partial_B_common_far_input(:,2)*1.96;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%% Pre vs post reward zone mean speed difference

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Pre post speed A zone on A laps - RF';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_sem_out.A.mean(1,1);
sem_txt = mean_sem_out.A.sem(1,1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_sem_out.A.mean(1,2);
sem_txt = mean_sem_out.A.sem(1,2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Pre post speed A zone on A laps - Random';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_sem_out.A.mean(4,1);
sem_txt = mean_sem_out.A.sem(4,1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_sem_out.A.mean(4,2);
sem_txt = mean_sem_out.A.sem(4,2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);


%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Pre post speed A zone on B laps - 5A5B';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_sem_out.ArelB.mean(2,1);
sem_txt = mean_sem_out.ArelB.sem(2,1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_sem_out.ArelB.mean(2,2);
sem_txt = mean_sem_out.ArelB.sem(2,2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Pre post speed A zone on B laps - Random';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_sem_out.ArelB.mean(4,1);
sem_txt = mean_sem_out.ArelB.sem(4,1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_sem_out.ArelB.mean(4,2);
sem_txt = mean_sem_out.ArelB.sem(4,2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Pre post speed B zone on B laps - RF';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_sem_out.B.mean(1,1);
sem_txt = mean_sem_out.B.sem(1,1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_sem_out.B.mean(1,2);
sem_txt = mean_sem_out.B.sem(1,2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Pre post speed B zone on B laps - Random';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_sem_out.B.mean(4,1);
sem_txt = mean_sem_out.B.sem(4,1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_sem_out.B.mean(4,2);
sem_txt = mean_sem_out.B.sem(4,2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);


%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Pre post speed B zone on A laps - 5A5B';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_sem_out.BrelA.mean(2,1);
sem_txt = mean_sem_out.BrelA.sem(2,1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_sem_out.BrelA.mean(2,2);
sem_txt = mean_sem_out.BrelA.sem(2,2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'Pre post speed B zone on A laps - Random';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_sem_out.BrelA.mean(4,1);
sem_txt = mean_sem_out.BrelA.sem(4,1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_sem_out.BrelA.mean(4,2);
sem_txt = mean_sem_out.BrelA.sem(4,2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);
%newline
writeWordEnter(ActXWord,WordHandle,1);

%% Speed difference sup Fig 2

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'A-sel speed diff between laps';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);

var_in = source_data_sup_2_3.event_speed_plot.Asel_speed_diff;
mean_txt = nanmean(var_in);
sem_txt = nanstd(var_in, 0, 1)./sqrt(sum(~isnan(var_in)));

write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'B-sel speed diff between laps';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);

var_in = source_data_sup_2_3.event_speed_plot.Bsel_speed_diff;
mean_txt = nanmean(var_in);
sem_txt = nanstd(var_in, 0, 1)./sqrt(sum(~isnan(var_in)));

write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%% Place field count - single

txt_input = 'PF count single A vs. B';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_count.Asel_norm(:,1);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_count.Bsel_norm(:,1);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'PF count single A vs. AB-A';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_count.Asel_norm(:,1);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_count.AB.A_norm(:,1);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'PF count single A vs. AB-B';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_count.Asel_norm(:,1);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_count.AB.B_norm(:,1);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%% Place field count - double

txt_input = 'PF count double A vs. B';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_count.Asel_norm(:,2);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_count.Bsel_norm(:,2);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'PF count double A vs. AB-A';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_count.Asel_norm(:,2);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_count.AB.A_norm(:,2);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'PF count double A vs. AB-B';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_count.Asel_norm(:,2);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_count.AB.B_norm(:,2);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);


%% Place field count - triple

txt_input = 'PF count triple A vs. B';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_count.Asel_norm(:,3);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_count.Bsel_norm(:,3);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'PF count triple A vs. AB-A';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_count.Asel_norm(:,3);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);
txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_count.AB.A_norm(:,3);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'PF count triple A vs. AB-B';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_count.Asel_norm(:,3);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_count.AB.B_norm(:,3);
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);


%% Place field width


txt_input = 'PF width A vs. B';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_width_pool.Asel;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_width_pool.Bsel;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);






%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'PF width A vs. AB-A';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_width_pool.Asel;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_width_pool.AB.A;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'PF width A vs. AB-B';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_width_pool.Asel;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_width_pool.AB.B;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'PF width B vs. AB-A';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_width_pool.Bsel;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_width_pool.AB.A;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);


%newline
writeWordEnter(ActXWord,WordHandle,1);

txt_input = 'PF width B vs. AB-B';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_width_pool.Bsel;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_width_pool.AB.B;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);


txt_input = 'PF width AB-A vs. AB-B';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
var_in = source_data_sup_2_3.pf_width_pool.AB.A;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
var_in = source_data_sup_2_3.pf_width_pool.AB.B;
mean_txt = mean(var_in);
sem_txt =  std(var_in,0,1)./sqrt(numel(var_in));
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);



%% Place field width
source_data_sup_2_3.pf_width_pool.Asel
source_data_sup_2_3.pf_width_pool.Bsel

source_data_sup_2_3.pf_width_pool.AB.A
source_data_sup_2_3.pf_width_pool.AB.B
