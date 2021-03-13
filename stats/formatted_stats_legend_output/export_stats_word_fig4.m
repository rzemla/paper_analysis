%% Import all table formatted entries here for data generated for each figure

%import data for each figure
laptop_access = 1;

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
fig4_5_data = load('fig4_5_table_data.mat');

%% Start Word document that will contain the formatted stats data

WordFileName='legend_stats_formatted_fig4_5.doc';
CurDir=pwd;
FileSpec = fullfile(CurDir,WordFileName);
%active X handle for manipulating document (ActXWord)
[ActXWord,WordHandle]=StartWord(FileSpec);

fprintf('Document will be saved in %s\n',FileSpec);

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

%CONITNUE HERE

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



%% 3d - paired Wilcoxon test comparisons of each class

%input into paired Wilcoxon test
dof = table2array(fig3_data.table_list.t_frac_remap(1,7));
test_stat = table2array(fig3_data.table_list.t_frac_remap(1,8));
p_val = table2array(fig3_data.table_list.t_frac_remap(1,10));
sample_n = table2array(fig3_data.table_list.t_frac_remap(1,5));

%description of comparison
comp_descrip = 'Common vs. activity';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolor
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired Wilcoxon test
dof = table2array(fig3_data.table_list.t_frac_remap(2,7));
test_stat = table2array(fig3_data.table_list.t_frac_remap(2,8));
p_val = table2array(fig3_data.table_list.t_frac_remap(2,10));
sample_n = table2array(fig3_data.table_list.t_frac_remap(2,5));

%description of comparison
comp_descrip = 'Common vs. global';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolor
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input into paired Wilcoxon test
dof = table2array(fig3_data.table_list.t_frac_remap(3,7));
test_stat = table2array(fig3_data.table_list.t_frac_remap(3,8));
p_val = table2array(fig3_data.table_list.t_frac_remap(3,10));
sample_n = table2array(fig3_data.table_list.t_frac_remap(3,5));

%description of comparison
comp_descrip = 'Common vs. partial';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%close parenthesis
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% 3e - Distribution of place field centroids among common remapping neurons

%description of statistics
txt_input = 'Fig 3e - Distribution of place field centroid location in common remapping neurons';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%rayleigh input test parameters
dof = table2array(fig3_data.table_list.t_pf_dist_com(1,7));
test_stat = table2array(fig3_data.table_list.t_pf_dist_com(1,8));
p_val = table2array(fig3_data.table_list.t_pf_dist_com(1,9));
sample_n = table2array(fig3_data.table_list.t_pf_dist_com(1,5));
%description of comparison
comp_descrip = 'Place field distribution for common neurons';

writeRayleighTest(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%close parenthesis
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% 3f - Zone shift for A vs B place fields for global remapping neurons

%description of statistics
txt_input = 'Fig 3f - Zone shift for A vs B place fields for global remapping neurons';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input to 1-sample paired Wilcoxon test
dof = table2array(fig3_data.table_list.t_zone_shift(1,7));
test_stat = table2array(fig3_data.table_list.t_zone_shift(1,8));
p_val = table2array(fig3_data.table_list.t_zone_shift(1,9));
sample_n = table2array(fig3_data.table_list.t_zone_shift(1,5));

%description of comparison
comp_descrip = 'Zone I A vs. B lap shifts';
writeOneSampleWilcoxAnimal_0(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolor
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input to 1-sample paired Wilcoxon test
dof = table2array(fig3_data.table_list.t_zone_shift(2,7));
test_stat = table2array(fig3_data.table_list.t_zone_shift(2,8));
p_val = table2array(fig3_data.table_list.t_zone_shift(2,9));
sample_n = table2array(fig3_data.table_list.t_zone_shift(2,5));

%description of comparison
comp_descrip = 'Zone II A vs. B lap shifts';
writeOneSampleWilcoxAnimal_0(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolor
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%input to 1-sample paired Wilcoxon test
dof = table2array(fig3_data.table_list.t_zone_shift(3,7));
test_stat = table2array(fig3_data.table_list.t_zone_shift(3,8));
p_val = table2array(fig3_data.table_list.t_zone_shift(3,9));
sample_n = table2array(fig3_data.table_list.t_zone_shift(3,5));

%description of comparison
comp_descrip = 'Zone III A vs. B lap shifts';
writeOneSampleWilcoxAnimal_0(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%close parenthesis
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% 3g - Global remapping neuron shift between reward zones - significant trends between zones

%description of statistics
txt_input = 'Fig 3g - Global remapping neurons switching fields between A and B laps between zones';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);


%insert Friedman test for all zone comparisons for global remappers
%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%Friedman test input data
dof = table2array(fig3_data.table_list.t_frac_zones_friedman(1,7));
test_stat = table2array(fig3_data.table_list.t_frac_zones_friedman(1,8));
p_val = table2array(fig3_data.table_list.t_frac_zones_friedman(1,9));
sample_n = table2array(fig3_data.table_list.t_frac_zones_friedman(1,5));

%description of comparison
comp_descrip = 'Global remapping neurons place field zone shift';
%Fig. 3f friedman test write
writeFriedmanTest(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%split entry with semicolon
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);


%insert individual paired Wilcoxon corrections
%II-III and I-II - trend, are statistically significant, but after
%correcting for comparisons across all zones, they no longer meet
%statistical significance - put non-correct p value in parentheses

%AII to BIII vs. AIII to BII
dof = table2array(fig3_data.table_list.t_frac_zone_remap(5,7));
test_stat = table2array(fig3_data.table_list.t_frac_zone_remap(5,8));
p_val = table2array(fig3_data.table_list.t_frac_zone_remap(5,10));
p_val_non_cor = table2array(fig3_data.table_list.t_frac_zone_remap(5,9));
sample_n = table2array(fig3_data.table_list.t_frac_zone_remap(5,5));

%variant function that places non-corrected p val in parenthesis (fig 3)
comp_descrip = 'AII to BIII vs. AIII to BII shift';
writePairedWilcoxAnimal_non_cor_pval(ActXWord,WordHandle,comp_descrip,test_stat,p_val,p_val_non_cor, dof, sample_n)

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%AI to BII vs. AII to BI
dof = table2array(fig3_data.table_list.t_frac_zone_remap(4,7));
test_stat = table2array(fig3_data.table_list.t_frac_zone_remap(4,8));
p_val = table2array(fig3_data.table_list.t_frac_zone_remap(4,10));
p_val_non_cor = table2array(fig3_data.table_list.t_frac_zone_remap(4,9));
sample_n = table2array(fig3_data.table_list.t_frac_zone_remap(4,5));

%variant function that places non-corrected p val in parenthesis (fig 3)
comp_descrip = 'AI to BII vs. AII to BI shift';
writePairedWilcoxAnimal_non_cor_pval(ActXWord,WordHandle,comp_descrip,test_stat,p_val,p_val_non_cor, dof, sample_n)

%close parenthesis
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Figure 3g - Non-significant cross zone shifts and intrazonal shifts stats

%description of statistics
txt_input = 'Fig 3g - Global remapping neurons intra-zone shifts and nonsig cross zonal';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%AI-BI vs. AII-BII
dof = table2array(fig3_data.table_list.t_frac_zone_remap(1,7));
test_stat = table2array(fig3_data.table_list.t_frac_zone_remap(1,8));
p_val = table2array(fig3_data.table_list.t_frac_zone_remap(1,10));
sample_n = table2array(fig3_data.table_list.t_frac_zone_remap(1,5));

%description of comparison
comp_descrip = 'AI-BI vs. AII-BII shift';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%AI-BI vs. AIII-BIII
dof = table2array(fig3_data.table_list.t_frac_zone_remap(2,7));
test_stat = table2array(fig3_data.table_list.t_frac_zone_remap(2,8));
p_val = table2array(fig3_data.table_list.t_frac_zone_remap(2,10));
sample_n = table2array(fig3_data.table_list.t_frac_zone_remap(2,5));

%description of comparison
comp_descrip = 'AI-BI vs. AIII-BIII shift';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%AII-BII vs. AIII-BIII
dof = table2array(fig3_data.table_list.t_frac_zone_remap(3,7));
test_stat = table2array(fig3_data.table_list.t_frac_zone_remap(3,8));
p_val = table2array(fig3_data.table_list.t_frac_zone_remap(3,10));
sample_n = table2array(fig3_data.table_list.t_frac_zone_remap(3,5));

%description of comparison
comp_descrip = 'AII-BII vs. AIII-BIII shift';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%inset semicolon break
txt_input = '; ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%AI-BIII vs. AIII-BI
dof = table2array(fig3_data.table_list.t_frac_zone_remap(6,7));
test_stat = table2array(fig3_data.table_list.t_frac_zone_remap(6,8));
p_val = table2array(fig3_data.table_list.t_frac_zone_remap(6,10));
sample_n = table2array(fig3_data.table_list.t_frac_zone_remap(6,5));

%description of comparison
comp_descrip = 'AI-BIII vs. AIII-BI shift';
writePairedWilcoxAnimal(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%close parenthesis
txt_input = ')';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,2);

%% Fig 3h - Distribution of partially remapping place field on A vs B trials

%description of statistics
txt_input = 'Fig 3h - distribution of partially remapping place fields on A vs. B trials';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
writeWordEnter(ActXWord,WordHandle,1);

%open parenthesis
txt_input = '(';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

dof = table2array(fig3_data.table_list.t_2ks_partial_pf_dist(1,7));
test_stat = table2array(fig3_data.table_list.t_2ks_partial_pf_dist(1,8));
p_val = table2array(fig3_data.table_list.t_2ks_partial_pf_dist(1,9));
sample_n = table2array(fig3_data.table_list.t_2ks_partial_pf_dist(1,5));

%description of comparison
comp_descrip = 'A vs. B partial remapping neurons place field centroid difference';
%modify this to run 2KS test
write2ksTest(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)

%close parenthesis
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
