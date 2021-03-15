%% Import all table formatted entries here for data generated for each figure

%import data for each figure
laptop_access = 0;

%laptop path directory
laptop_path_dir = 'C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data';
%desktop path directory
desktop_path_dir = 'G:\Google_drive\task_selective_place_paper\matlab_data';

%load table data
if laptop_access == 1
cd(laptop_path_dir)
else
    cd(desktop_path_dir)
end

%load in figure 4 and 5 data
ex_fig_data = load('ex_fig_table_data.mat');

%% Start Word document that will contain the formatted stats data

WordFileName='legend_stats_formatted_extended_fig.doc';
CurDir=pwd;
FileSpec = fullfile(CurDir,WordFileName);
%active X handle for manipulating document (ActXWord)
[ActXWord,WordHandle]=StartWord(FileSpec);

fprintf('Document will be saved in %s\n',FileSpec);


%% Ex Fig. 4b A-B place field speed difference for and A and B selective neurons
%description of statistics
txt_input = 'Ex. Fig. 4b - A-B in-field speed difference for A and B selective neurons';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input to 1-sample paired Wilcoxon test
dof = table2array(fig3_data.table_list.t_zone_shift(2,7));
test_stat = table2array(fig3_data.table_list.t_zone_shift(2,8));
p_val = table2array(fig3_data.table_list.t_zone_shift(2,9));
sample_n = table2array(fig3_data.table_list.t_zone_shift(2,5));

%description of comparison
comp_descrip = 'Zone II A vs. B lap shifts';
writeOneSampleWilcoxAnimal_0(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)


%% Figure 4e - Fraction of A and B-selective place cells during learning
%description of statistics
txt_input = 'Fig 4e - Fraction of A and B-selective place cells during learning';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%1-way LME input data - extract DOFs from test stats field
dof = table2array(fig4_5_data.table_list.t_1_rm_lme.learnA (1,7));
test_stat = table2array(fig4_5_data.table_list.t_1_rm_lme.learnA (1,8));
p_val = table2array(fig4_5_data.table_list.t_1_rm_lme.learnA (1,9));
sample_n = table2array(fig4_5_data.table_list.t_1_rm_lme.learnA (1,5));

comp_descrip = 'Fraction of A-trial tuned place cells during learning - T.S.';

%1-way RM linear mixed effects analysis
write_1wayLME(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.learnA_6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.learnA_6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.learnA_6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.learnA_6_7(1,5));

%description of comparison
comp_descrip = 'Day 1 vs. Day 6 A trial learning';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.learnA_6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.learnA_6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.learnA_6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.learnA_6_7(2,5));

%description of comparison
comp_descrip = 'Day 1 vs. Day 7 A trial learning';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);


%B trial group test
%1-way LME input data - extract DOFs from test stats field
dof = table2array(fig4_5_data.table_list.t_1_rm_lme.learnB(1,7));
test_stat = table2array(fig4_5_data.table_list.t_1_rm_lme.learnB(1,8));
p_val = table2array(fig4_5_data.table_list.t_1_rm_lme.learnB(1,9));
sample_n = table2array(fig4_5_data.table_list.t_1_rm_lme.learnB(1,5));

comp_descrip = 'Fraction of B-trial tuned place cells during learning - T.S.';

%1-way RM linear mixed effects analysis
write_1wayLME(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.learnB_6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.learnB_6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.learnB_6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.learnB_6_7(1,5));

%description of comparison
comp_descrip = 'Day 1 vs. Day 6 B trial learning';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.learnB_6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.learnB_6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.learnB_6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.learnB_6_7(2,5));

%description of comparison
comp_descrip = 'Day 1 vs. Day 7 B trial learning';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%close parenthesis
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);


%% Figure 4e - Fraction of A and B-selective place cells during recall
%description of statistics
txt_input = 'Fig 4e - Fraction of A and B-selective place cells during recall';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%1-way LME input data - extract DOFs from test stats field
dof = table2array(fig4_5_data.table_list.t_1_rm_lme.recallA(1,7));
test_stat = table2array(fig4_5_data.table_list.t_1_rm_lme.recallA(1,8));
p_val = table2array(fig4_5_data.table_list.t_1_rm_lme.recallA(1,9));
sample_n = table2array(fig4_5_data.table_list.t_1_rm_lme.recallA (1,5));

comp_descrip = 'Fraction of A-trial tuned place cells during recall - T.S.';

%1-way RM linear mixed effects analysis
write_1wayLME(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.recallA_6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.recallA_6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.recallA_6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.recallA_6_7(1,5));

%description of comparison
comp_descrip = 'Day 1 vs. Day 6 A trial recall';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.recallA_6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.recallA_6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.recallA_6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.recallA_6_7(2,5));

%description of comparison
comp_descrip = 'Day 1 vs. Day 7 A trial recall';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);


%B trial group test
%1-way LME input data - extract DOFs from test stats field
dof = table2array(fig4_5_data.table_list.t_1_rm_lme.recallB(1,7));
test_stat = table2array(fig4_5_data.table_list.t_1_rm_lme.recallB(1,8));
p_val = table2array(fig4_5_data.table_list.t_1_rm_lme.recallB(1,9));
sample_n = table2array(fig4_5_data.table_list.t_1_rm_lme.recallB(1,5));

comp_descrip = 'Fraction of B-trial tuned place cells during recall - T.S.';

%1-way RM linear mixed effects analysis
write_1wayLME(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.recallB_6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.recallB_6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.recallB_6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.recallB_6_7(1,5));

%description of comparison
comp_descrip = 'Day 1 vs. Day 6 B trial recall';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.recallB_6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.recallB_6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.recallB_6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.recallB_6_7(2,5));

%description of comparison
comp_descrip = 'Day 1 vs. Day 7 B trial recall';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%close parenthesis
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Fig 4f - PV correlation across days  - 2-way LME - Recall - A and B trials

%description of statistics
txt_input = 'Fig 4f - PV correlation relative to day 1 on A and B trials during recall';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%2-way LME input data - extract DOFs from test stats field
dof = table2array(fig4_5_data.table_list.t_2way_rm_lme.PV_learn_recallA(1,7));
test_stat = table2array(fig4_5_data.table_list.t_2way_rm_lme.PV_learn_recallA(1,8));
p_val = table2array(fig4_5_data.table_list.t_2way_rm_lme.PV_learn_recallA(1,9));
sample_n = table2array(fig4_5_data.table_list.t_2way_rm_lme.PV_learn_recallA(1,5));

comp_descrip = 'PV correlation on A trials';

%2-way RM linear mixed effects analysis
write_2wayLME(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%paired t-tests
%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.PV_recallA2v6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.PV_recallA2v6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.PV_recallA2v6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.PV_recallA2v6_7(1,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 6 A trial recall';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.PV_recallA2v6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.PV_recallA2v6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.PV_recallA2v6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.PV_recallA2v6_7(2,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 7 A trial recall';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%B trials
%2-way LME input data - extract DOFs from test stats field
dof = table2array(fig4_5_data.table_list.t_2way_rm_lme.PV_learn_recallB(1,7));
test_stat = table2array(fig4_5_data.table_list.t_2way_rm_lme.PV_learn_recallB(1,8));
p_val = table2array(fig4_5_data.table_list.t_2way_rm_lme.PV_learn_recallB(1,9));
sample_n = table2array(fig4_5_data.table_list.t_2way_rm_lme.PV_learn_recallB(1,5));

comp_descrip = 'PV correlation on B trials';

%2-way RM linear mixed effects analysis
write_2wayLME(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%paired t-tests
%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.PV_recallB2v6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.PV_recallB2v6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.PV_recallB2v6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.PV_recallB2v6_7(1,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 6 B trial recall';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.PV_recallB2v6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.PV_recallB2v6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.PV_recallB2v6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.PV_recallB2v6_7(2,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 7 B trial recall';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%close parenthesis
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Fig 4f - PV learning and learning vs recall t tests

%description of statistics
txt_input = 'Fig 4f - Paired and unpaired t-test for learning and learning vs recall';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.PV_learnA2v6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.PV_learnA2v6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.PV_learnA2v6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.PV_learnA2v6_7(1,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 6 A trial learning';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.PV_learnA2v6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.PV_learnA2v6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.PV_learnA2v6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.PV_learnA2v6_7(2,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 7 A trial learning';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

% B trials
%paired t-tests
%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.PV_learnB2v6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.PV_learnB2v6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.PV_learnB2v6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.PV_learnB2v6_7(1,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 6 B trial learning';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.PV_learnB2v6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.PV_learnB2v6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.PV_learnB2v6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.PV_learnB2v6_7(2,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 7 B trial learning';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into unpaired
dof = table2array(fig4_5_data.table_list.t_ttest.learn_recallA_6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.learn_recallA_6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.learn_recallA_6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.learn_recallA_6_7(1,5));

%description of comparison
comp_descrip = 'Day 6 A trials learning vs. recall';
writeUnPairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into unpaired
dof = table2array(fig4_5_data.table_list.t_ttest.learn_recallA_6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.learn_recallA_6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.learn_recallA_6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.learn_recallA_6_7(2,5));

%description of comparison
comp_descrip = 'Day 7 A trials learning vs. recall';
writeUnPairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into unpaired
dof = table2array(fig4_5_data.table_list.t_ttest.learn_recallB_6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.learn_recallB_6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.learn_recallB_6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.learn_recallB_6_7(1,5));

%description of comparison
comp_descrip = 'Day 6 B trials learning vs. recall';
writeUnPairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into unpaired
dof = table2array(fig4_5_data.table_list.t_ttest.learn_recallB_6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.learn_recallB_6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.learn_recallB_6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.learn_recallB_6_7(2,5));

%description of comparison
comp_descrip = 'Day 7 B trials learning vs. recall';
writeUnPairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%close parenthesis
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Fig 4f Tuning curve correction - TS - equivalent analysis - learning then recall

%description of statistics
txt_input = 'Fig 4f - TC correlation (TS) relative to day 1 on A and B trials during learning';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%2-way LME input data - extract DOFs from test stats field
dof = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_learn_recallA(1,7));
test_stat = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_learn_recallA(1,8));
p_val = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_learn_recallA(1,9));
sample_n = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_learn_recallA(1,5));

comp_descrip = 'TC correlation on A trials';

%2-way RM linear mixed effects analysis
write_2wayLME(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%separate with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%Learn A 2 vs 6 and 2 vs 7

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnA2v6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnA2v6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnA2v6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnA2v6_7(1,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 6 A trial learning';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnA2v6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnA2v6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnA2v6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnA2v6_7(2,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 7 A trial learning';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%separate with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%%%%

%B trials
%2-way LME input data - extract DOFs from test stats field
dof = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_learn_recallB(1,7));
test_stat = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_learn_recallB(1,8));
p_val = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_learn_recallB(1,9));
sample_n = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_learn_recallB(1,5));

comp_descrip = 'TC correlation on B trials';

%2-way RM linear mixed effects analysis
write_2wayLME(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%separate with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%Learn A 2 vs 6 and 2 vs 7

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnB2v6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnB2v6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnB2v6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnB2v6_7(1,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 6 B trial learning';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnB2v6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnB2v6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnB2v6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learnB2v6_7(2,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 7 B trial learning';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%%%%%%
%rest of t-test stats comparing recall days and Learning vs Recall

%semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallA2v6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallA2v6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallA2v6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallA2v6_7(1,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 6 A trial recall';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallA2v6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallA2v6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallA2v6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallA2v6_7(2,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 7 A trial recall';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

% B trials
%paired t-tests
%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallB2v6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallB2v6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallB2v6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallB2v6_7(1,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 6 B trial recall';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallB2v6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallB2v6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallB2v6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_recallB2v6_7(2,5));

%description of comparison
comp_descrip = 'Day 2 vs. Day 7 B trial recall';
writePairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into unpaired
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallA_6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallA_6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallA_6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallA_6_7(1,5));

%description of comparison
comp_descrip = 'Day 6 A trials learning vs. recall';
writeUnPairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into unpaired
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallA_6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallA_6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallA_6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallA_6_7(2,5));

%description of comparison
comp_descrip = 'Day 7 A trials learning vs. recall';
writeUnPairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into unpaired
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallB_6_7(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallB_6_7(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallB_6_7(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallB_6_7(1,5));

%description of comparison
comp_descrip = 'Day 6 B trials learning vs. recall';
writeUnPairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into unpaired
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallB_6_7(2,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallB_6_7(2,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallB_6_7(2,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_learn_recallB_6_7(2,5));

%description of comparison
comp_descrip = 'Day 7 B trials learning vs. recall';
writeUnPairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%close parenthesis
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Fig 4g - Neighboring PV correlation for A and B trials

%description of statistics
txt_input = 'Fig 4g - Neighboring PV correlation for A and B trials learning and recall';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%2-way LME input data - extract DOFs from test stats field
dof = table2array(fig4_5_data.table_list.t_2way_rm_lme.PVn_learn_recallA(1,7));
test_stat = table2array(fig4_5_data.table_list.t_2way_rm_lme.PVn_learn_recallA(1,8));
p_val = table2array(fig4_5_data.table_list.t_2way_rm_lme.PVn_learn_recallA(1,9));
sample_n = table2array(fig4_5_data.table_list.t_2way_rm_lme.PVn_learn_recallA(1,5));

comp_descrip = 'Neighboring session PV correlation on A trials';

%2-way RM linear mixed effects analysis
write_2wayLME(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%separate with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%2-way LME input data - extract DOFs from test stats field
dof = table2array(fig4_5_data.table_list.t_2way_rm_lme.PVn_learn_recallB(1,7));
test_stat = table2array(fig4_5_data.table_list.t_2way_rm_lme.PVn_learn_recallB(1,8));
p_val = table2array(fig4_5_data.table_list.t_2way_rm_lme.PVn_learn_recallB(1,9));
sample_n = table2array(fig4_5_data.table_list.t_2way_rm_lme.PVn_learn_recallB(1,5));

comp_descrip = 'Neighboring session PV correlation on B trials';

%2-way RM linear mixed effects analysis
write_2wayLME(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%separate with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%Learning 2 paired t-tests(no Holm-Sidak correction for these)
%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.PVn_learnA12v67(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.PVn_learnA12v67(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.PVn_learnA12v67(1,9));
sample_n = table2array(fig4_5_data.table_list.t_ttest.PVn_learnA12v67(1,5));

%description of comparison
comp_descrip = 'Days 1 vs. 2 Vs. Day 6 vs. 7 A trials learn';
writePairedTtestAnimal_no_pcor(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.PVn_learnB12v67(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.PVn_learnB12v67(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.PVn_learnB12v67(1,9));
sample_n = table2array(fig4_5_data.table_list.t_ttest.PVn_learnB12v67(1,5));

%description of comparison
comp_descrip = 'Days 1 vs. 2 Vs. Day 6 vs. 7 B trials learn';
writePairedTtestAnimal_no_pcor(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%Recall 2 paired t-tests

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.PVn_recallA12v67(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.PVn_recallA12v67(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.PVn_recallA12v67(1,9));
sample_n = table2array(fig4_5_data.table_list.t_ttest.PVn_recallA12v67(1,5));

%description of comparison
comp_descrip = 'Days 1 vs. 2 Vs. Day 6 vs. 7 A trials recall';
writePairedTtestAnimal_no_pcor(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

dof = table2array(fig4_5_data.table_list.t_ttest.PVn_recallB12v67(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.PVn_recallB12v67(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.PVn_recallB12v67(1,9));
sample_n = table2array(fig4_5_data.table_list.t_ttest.PVn_recallB12v67(1,5));

%description of comparison
comp_descrip = 'Days 1 vs. 2 Vs. Day 6 vs. 7 B trials recall';
writePairedTtestAnimal_no_pcor(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);


%% Fig 4g - Neighboring TC (TS) correlation for A and B trials

%description of statistics
txt_input = 'Fig 4g - Neighboring TC (TS) correlation for A and B trials learning and recall';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%2-way LME input data - extract DOFs from test stats field
dof = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_n_learn_recallA(1,7));
test_stat = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_n_learn_recallA(1,8));
p_val = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_n_learn_recallA(1,9));
sample_n = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_n_learn_recallA(1,5));

comp_descrip = 'Neighboring session TC correlation on A trials';

%2-way RM linear mixed effects analysis
write_2wayLME(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%separate with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%2-way LME input data - extract DOFs from test stats field
dof = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_n_learn_recallB(1,7));
test_stat = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_n_learn_recallB(1,8));
p_val = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_n_learn_recallB(1,9));
sample_n = table2array(fig4_5_data.table_list.t_2way_rm_lme.TC_ts_n_learn_recallB(1,5));

comp_descrip = 'Neighboring session TC correlation on B trials';

%2-way RM linear mixed effects analysis
write_2wayLME(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%separate with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%Learning 2 paired t-tests(no Holm-Sidak correction for these)
%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_learnA12v67(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_learnA12v67(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_learnA12v67(1,9));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_learnA12v67(1,5));

%description of comparison
comp_descrip = 'Days 1 vs. 2 Vs. Day 6 vs. 7 A trials learn';
writePairedTtestAnimal_no_pcor(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_learnB12v67(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_learnB12v67(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_learnB12v67(1,9));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_learnB12v67(1,5));

%description of comparison
comp_descrip = 'Days 1 vs. 2 Vs. Day 6 vs. 7 B trials learn';
writePairedTtestAnimal_no_pcor(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%Recall 2 paired t-tests

%input into paired t-test
dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_recallA12v67(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_recallA12v67(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_recallA12v67(1,9));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_recallA12v67(1,5));

%description of comparison
comp_descrip = 'Days 1 vs. 2 Vs. Day 6 vs. 7 A trials recall';
writePairedTtestAnimal_no_pcor(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

dof = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_recallB12v67(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_recallB12v67(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_recallB12v67(1,9));
sample_n = table2array(fig4_5_data.table_list.t_ttest.TC_ts_n_recallB12v67(1,5));

%description of comparison
comp_descrip = 'Days 1 vs. 2 Vs. Day 6 vs. 7 B trials recall';
writePairedTtestAnimal_no_pcor(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Fig 4h  Mean centroid difference (TS) relative to d1 - A and B trials

%description of statistics
txt_input = 'Fig 4h - Mean centroid difference (TS) relative to day 1 on A and B trials during learning and recall';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%2-way LME input data - extract DOFs from test stats field
dof = table2array(fig4_5_data.table_list.t_2way_rm_lme.cent_ts_learn_recallA(1,7));
test_stat = table2array(fig4_5_data.table_list.t_2way_rm_lme.cent_ts_learn_recallA(1,8));
p_val = table2array(fig4_5_data.table_list.t_2way_rm_lme.cent_ts_learn_recallA(1,9));
sample_n = table2array(fig4_5_data.table_list.t_2way_rm_lme.cent_ts_learn_recallA(1,5));

comp_descrip = 'Mean centroid difference relative to Day 1 A trials';

%2-way RM linear mixed effects analysis
write_2wayLME(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%2-way LME input data - extract DOFs from test stats field
dof = table2array(fig4_5_data.table_list.t_2way_rm_lme.cent_ts_learn_recallB(1,7));
test_stat = table2array(fig4_5_data.table_list.t_2way_rm_lme.cent_ts_learn_recallB(1,8));
p_val = table2array(fig4_5_data.table_list.t_2way_rm_lme.cent_ts_learn_recallB(1,9));
sample_n = table2array(fig4_5_data.table_list.t_2way_rm_lme.cent_ts_learn_recallB(1,5));

comp_descrip = 'Mean centroid difference relative to Day 1 B trials';

%2-way RM linear mixed effects analysis
write_2wayLME(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into unpaired
dof = table2array(fig4_5_data.table_list.t_ttest.cent_ts_learn_recallA_15_16(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.cent_ts_learn_recallA_15_16(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.cent_ts_learn_recallA_15_16(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.cent_ts_learn_recallA_15_16(1,5));

%description of comparison
comp_descrip = 'Day 5 A trials learning vs. recall';
writeUnPairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into unpaired
dof = table2array(fig4_5_data.table_list.t_ttest.cent_ts_learn_recallB_15_16(1,7));
test_stat = table2array(fig4_5_data.table_list.t_ttest.cent_ts_learn_recallB_15_16(1,8));
p_val = table2array(fig4_5_data.table_list.t_ttest.cent_ts_learn_recallB_15_16(1,10));
sample_n = table2array(fig4_5_data.table_list.t_ttest.cent_ts_learn_recallB_15_16(1,5));

%description of comparison
comp_descrip = 'Day 5 B trials learning vs. recall';
writeUnPairedTtestAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Fig. 5c - A&B tuned A vs B lap spatial tuning correlation relative to D1 during learning

txt_input = 'Fig 5c - Matching A&B tuned neurons A vs B lap correlation normalized to D1 during learning';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%Friedman test input data
dof = table2array(fig4_5_data.table_list.t_krusall_learn(1,7));
test_stat = table2array(fig4_5_data.table_list.t_krusall_learn(1,8));
p_val = table2array(fig4_5_data.table_list.t_krusall_learn(1,9));
sample_n = table2array(fig4_5_data.table_list.t_krusall_learn(1,5));

%description of comparison
comp_descrip = 'Day one normalized A vs. B lap correlation scores for matching neurons during learning';
%kruskal wallis test write
writeKruskalWallisTest(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%1-sample Wilcox test no multi comp p_val correction against 1
%input to 1-sample paired Wilcoxon test
dof = table2array(fig4_5_data.table_list.paired_wilcox_AB_learn(1,7));
test_stat = table2array(fig4_5_data.table_list.paired_wilcox_AB_learn(1,8));
p_val = table2array(fig4_5_data.table_list.paired_wilcox_AB_learn(1,9));
sample_n = table2array(fig4_5_data.table_list.paired_wilcox_AB_learn(1,5));

%description of comparison
comp_descrip = 'Day 2 learn';
writeOneSampleWilcoxAnimal_1(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%1-sample Wilcox test no multi comp p_val correction against 1
%input to 1-sample paired Wilcoxon test
dof = table2array(fig4_5_data.table_list.paired_wilcox_AB_learn(6,7));
test_stat = table2array(fig4_5_data.table_list.paired_wilcox_AB_learn(6,8));
p_val = table2array(fig4_5_data.table_list.paired_wilcox_AB_learn(6,9));
sample_n = table2array(fig4_5_data.table_list.paired_wilcox_AB_learn(6,5));

%description of comparison
comp_descrip = 'Day 7 learn';
writeOneSampleWilcoxAnimal_1(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Fig. 5d - A&B tuned A vs B lap spatial tuning correlation relative to D1 - recall

txt_input = 'Fig 5d - Matching A&B tuned neurons A vs B lap correlation normalized to D1 during recall';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%Kruskal-Wallis test input data
dof = table2array(fig4_5_data.table_list.t_krusall_recall(1,7));
test_stat = table2array(fig4_5_data.table_list.t_krusall_recall(1,8));
p_val = table2array(fig4_5_data.table_list.t_krusall_recall(1,9));
sample_n = table2array(fig4_5_data.table_list.t_krusall_recall(1,5));

%description of comparison
comp_descrip = 'Day one normalized A vs. B lap correlation scores for matching neurons during recall';
%kruskal wallis test write
writeKruskalWallisTest(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%1-sample Wilcox test no multi comp p_val correction against 1
%input to 1-sample paired Wilcoxon test
dof = table2array(fig4_5_data.table_list.paired_wilcox_AB_recall(1,7));
test_stat = table2array(fig4_5_data.table_list.paired_wilcox_AB_recall(1,8));
p_val = table2array(fig4_5_data.table_list.paired_wilcox_AB_recall(1,9));
sample_n = table2array(fig4_5_data.table_list.paired_wilcox_AB_recall(1,5));

%description of comparison
comp_descrip = 'Day 2 recall';
writeOneSampleWilcoxAnimal_1(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%1-sample Wilcox test no multi comp p_val correction against 1
%input to 1-sample paired Wilcoxon test
dof = table2array(fig4_5_data.table_list.paired_wilcox_AB_recall(6,7));
test_stat = table2array(fig4_5_data.table_list.paired_wilcox_AB_recall(6,8));
p_val = table2array(fig4_5_data.table_list.paired_wilcox_AB_recall(6,9));
sample_n = table2array(fig4_5_data.table_list.paired_wilcox_AB_recall(6,5));

%description of comparison
comp_descrip = 'Day 7 recall';
writeOneSampleWilcoxAnimal_1(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Close Word document
CloseWord(ActXWord,WordHandle,FileSpec);


%order of entry
%test name, test statistic, p value (include rounding and sigstar add),
%test statistic - 2 decimal places
%p-value 3 decimal places and sub to <0.001 for low p values
%nb of samples as a char entry

%% Original function customized
%WriteToWordFromMatlab_testing
