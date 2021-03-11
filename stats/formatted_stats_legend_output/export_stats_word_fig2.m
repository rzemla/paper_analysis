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

%load in figure 2 data
fig2_data = load('fig2_table_data.mat');

%% Make Word export stats script for each figure legend

%% Start Word document that will contain the formatted stats data

WordFileName='legend_stats_formatted.doc';
CurDir=pwd;
FileSpec = fullfile(CurDir,WordFileName);
%active X handle for manipulating document (ActXWord)
[ActXWord,WordHandle]=StartWord(FileSpec);

fprintf('Document will be saved in %s\n',FileSpec);

%% Figure 2c AUC/min comparison for A and B selective place cells - RUN

%description of statistics
txt_input = 'Fig 2c - AUC/min RUN epochs';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%A-selective place cells (1) - adjusted p values here
dof = table2array(fig2_data.table_list.t_auc_run(1,7));
test_stat = table2array(fig2_data.table_list.t_auc_run(1,8));
p_val = table2array(fig2_data.table_list.t_auc_run(1,10));
sample_n = table2array(fig2_data.table_list.t_auc_run(1,5));

%description of comparison
comp_descrip = 'A-selective A vs. B lap activity rate';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%B-selective place cells (1)
dof = table2array(fig2_data.table_list.t_auc_run(2,7));
test_stat = table2array(fig2_data.table_list.t_auc_run(2,8));
p_val = table2array(fig2_data.table_list.t_auc_run(2,10));
sample_n = table2array(fig2_data.table_list.t_auc_run(2,5));

%description of comparison
comp_descrip = 'B-selective A vs. B lap activity rate';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%A&B-selective (non-selective) place cells (1)
dof = table2array(fig2_data.table_list.t_auc_run(3,7));
test_stat = table2array(fig2_data.table_list.t_auc_run(3,8));
p_val = table2array(fig2_data.table_list.t_auc_run(3,10));
sample_n = table2array(fig2_data.table_list.t_auc_run(3,5));

%description of comparison
comp_descrip = 'Non-selective A vs. B lap activity rate';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%close parenthesis and 2 newlines
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Figure 2c AUC/min comparison - NO RUN epochs

%description of statistics
txt_input = 'Fig 2c - AUC/min NO RUN epochs';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%A-selective place cells (1) - adjusted p values here
dof = table2array(fig2_data.table_list.t_auc_norun(1,7));
test_stat = table2array(fig2_data.table_list.t_auc_norun(1,8));
p_val = table2array(fig2_data.table_list.t_auc_norun(1,10));
sample_n = table2array(fig2_data.table_list.t_auc_norun(1,5));

%description of comparison
comp_descrip = 'A-selective A vs. B lap activity rate';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%B-selective place cells (1)
dof = table2array(fig2_data.table_list.t_auc_norun(2,7));
test_stat = table2array(fig2_data.table_list.t_auc_norun(2,8));
p_val = table2array(fig2_data.table_list.t_auc_norun(2,10));
sample_n = table2array(fig2_data.table_list.t_auc_norun(2,5));

%description of comparison
comp_descrip = 'B-selective A vs. B lap activity rate';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%A&B-selective (non-selective) place cells (1)
dof = table2array(fig2_data.table_list.t_auc_norun(3,7));
test_stat = table2array(fig2_data.table_list.t_auc_norun(3,8));
p_val = table2array(fig2_data.table_list.t_auc_norun(3,10));
sample_n = table2array(fig2_data.table_list.t_auc_norun(3,5));

%description of comparison
comp_descrip = 'Non-selective A vs. B lap activity rate';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%close parenthesis and 2 newlines
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);


%% Figure 2d - Fraction of task-selective neurons by each tuning criterion

%description of statistics
txt_input = 'Fig 2d - SI/TS fraction tuned Friedman test';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%Fig 2d - SI criterion Friedman (non-adjusted p vals (no comp)
%get input data for Friedman test:
dof = table2array(fig2_data.table_list.t_frac_si(1,7));
test_stat = table2array(fig2_data.table_list.t_frac_si(1,8));
p_val = table2array(fig2_data.table_list.t_frac_si(1,9));
sample_n = table2array(fig2_data.table_list.t_frac_si(1,5));
%description of comparison
comp_descrip = 'S.I. criterion';
%Fig. 2d friedman test write
writeFriedmanTest(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolor
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%Fig 2d - TS criterion Friedman
%get input data for Friedman test:
dof = table2array(fig2_data.table_list.t_frac_ts(1,7));
test_stat = table2array(fig2_data.table_list.t_frac_ts(1,8));
p_val = table2array(fig2_data.table_list.t_frac_ts(1,9));
sample_n = table2array(fig2_data.table_list.t_frac_ts(1,5));
%description of comparison
comp_descrip = 'T.S. criterion';
%Fig. 2d friedman test write
writeFriedmanTest(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)
%close parenthesis
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Fig 2d A vs B comparison by SI and TS criterion

%description of statistics
txt_input = 'Fig 2d - A vs B by SI and TS criterion';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%A vs B by SI criterion
dof = table2array(fig2_data.table_list.t_frac_si(2,7));
test_stat = table2array(fig2_data.table_list.t_frac_si(2,8));
p_val = table2array(fig2_data.table_list.t_frac_si(2,10));
sample_n = table2array(fig2_data.table_list.t_frac_si(2,5));

%description of comparison
comp_descrip = 'A vs. B - S.I. criterion';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n);

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%A vs B by TS criterion
dof = table2array(fig2_data.table_list.t_frac_ts(2,7));
test_stat = table2array(fig2_data.table_list.t_frac_ts(2,8));
p_val = table2array(fig2_data.table_list.t_frac_ts(2,10));
sample_n = table2array(fig2_data.table_list.t_frac_ts(2,5));

%description of comparison
comp_descrip = 'A vs. B - T.S. criterion';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n);

%closed parenthesis
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Fig. 2d A vs A&B and B vs A&B for SI and TS criteria

%description of statistics
txt_input = 'Fig 2d - A vs A&B and B vs A&B by SI and TS criteria';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%A vs A&BB by SI criterion
dof = table2array(fig2_data.table_list.t_frac_si(3,7));
test_stat = table2array(fig2_data.table_list.t_frac_si(3,8));
p_val = table2array(fig2_data.table_list.t_frac_si(3,10));
sample_n = table2array(fig2_data.table_list.t_frac_si(3,5));

%description of comparison
comp_descrip = 'A vs. A&B - S.I. criterion';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n);

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%B vs A&B by SI criterion
dof = table2array(fig2_data.table_list.t_frac_si(4,7));
test_stat = table2array(fig2_data.table_list.t_frac_si(4,8));
p_val = table2array(fig2_data.table_list.t_frac_si(4,10));
sample_n = table2array(fig2_data.table_list.t_frac_si(4,5));

%description of comparison
comp_descrip = 'B vs. A&B - S.I. criterion';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n);

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%A vs A&B by TS criterion
dof = table2array(fig2_data.table_list.t_frac_ts(3,7));
test_stat = table2array(fig2_data.table_list.t_frac_ts(3,8));
p_val = table2array(fig2_data.table_list.t_frac_ts(3,10));
sample_n = table2array(fig2_data.table_list.t_frac_ts(3,5));

%description of comparison
comp_descrip = 'A vs. A&B - T.S. criterion';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n);

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%B vs A&B by TS criterion
dof = table2array(fig2_data.table_list.t_frac_ts(4,7));
test_stat = table2array(fig2_data.table_list.t_frac_ts(4,8));
p_val = table2array(fig2_data.table_list.t_frac_ts(4,10));
sample_n = table2array(fig2_data.table_list.t_frac_ts(4,5));

%description of comparison
comp_descrip = 'B vs. A&B - T.S. criterion';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n);

%closed parenthesis
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Figure 2f - A and B selective Rayleigh tests of circular uniformity - finish modifying this

%description of statistics
txt_input = 'Fig 2f - A and B selective tests of uniformity';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

dof = table2array(fig2_data.table_list.t_pf_dist(1,7));
test_stat = table2array(fig2_data.table_list.t_pf_dist(1,8));
p_val = table2array(fig2_data.table_list.t_pf_dist(1,9));
sample_n = table2array(fig2_data.table_list.t_pf_dist(1,5));
%description of comparison
comp_descrip = 'A-selective place field distribution';

writeRayleighTest(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

dof = table2array(fig2_data.table_list.t_pf_dist(2,7));
test_stat = table2array(fig2_data.table_list.t_pf_dist(2,8));
p_val = table2array(fig2_data.table_list.t_pf_dist(2,9));
sample_n = table2array(fig2_data.table_list.t_pf_dist(2,5));
%description of comparison
comp_descrip = 'B-selective place field distribution';

writeRayleighTest(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%close parenthesis
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Figure 2g - A vs B place cell centroid difference - 2sample KS test

%description of statistics
txt_input = 'Fig 2g - A vs B place field centroid distribution';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

dof = table2array(fig2_data.table_list.t_2ks_pf_dist(1,7));
test_stat = table2array(fig2_data.table_list.t_2ks_pf_dist(1,8));
p_val = table2array(fig2_data.table_list.t_2ks_pf_dist(1,9));
sample_n = table2array(fig2_data.table_list.t_2ks_pf_dist(1,5));


%description of comparison
comp_descrip = 'A vs. B centroid difference';
%modify this to run 2KS test
write2ksTest(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%% Figure 2h - TC correlation between task selective/non selective place cells

%description of statistics
txt_input = 'Fig 2h - TC correlation between A vs B laps fo task selective/nonselective place cells';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%A vs B sel
sample_n = table2array(fig2_data.table_list.t_tc_corr(1,5));
dof = table2array(fig2_data.table_list.t_tc_corr(1,7));
test_stat = table2array(fig2_data.table_list.t_tc_corr(1,8));
p_val = table2array(fig2_data.table_list.t_tc_corr(1,10));

%description of comparison
comp_descrip = 'A- vs. B- selective';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n);

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%A vs AB
sample_n = table2array(fig2_data.table_list.t_tc_corr(2,5));
dof = table2array(fig2_data.table_list.t_tc_corr(2,7));
test_stat = table2array(fig2_data.table_list.t_tc_corr(2,8));
p_val = table2array(fig2_data.table_list.t_tc_corr(2,10));

%description of comparison
comp_descrip = 'A-selective vs. A&B';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n);

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);


%B vs AB
sample_n = table2array(fig2_data.table_list.t_tc_corr(3,5));
dof = table2array(fig2_data.table_list.t_tc_corr(3,7));
test_stat = table2array(fig2_data.table_list.t_tc_corr(3,8));
p_val = table2array(fig2_data.table_list.t_tc_corr(3,10));

%description of comparison
comp_descrip = 'B-selective vs. A&B';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n);


%% Import and embed the number of neurons in the above stats (last)


%% Insert two enters(newlines)
% writeWordEnter(ActXWord,WordHandle,2);
% 
% txt_input = 'Write text here';
% writeDefaultWordText(ActXWord,WordHandle,txt_input);

%% Close Word document
%CloseWord(ActXWord,WordHandle,FileSpec);


%% 
%chi character code: 1D712 - unicode hex value;
%convert hex to dec: 120594

%order of entry
%test name, test statistic, p value (include rounding and sigstar add),
%test statistic - 2 decimal places
%p-value 3 decimal places and sub to <0.001 for low p values
%nb of samples as a char entry

%% Original function customized
%WriteToWordFromMatlab_testing
