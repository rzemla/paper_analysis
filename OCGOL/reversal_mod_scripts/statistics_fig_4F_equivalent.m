%% Load data

load('D:\OCGOL_reversal\figure4_matrices.mat')

%datapoints binned by behavior phase
binned_data = load('D:\OCGOL_reversal\figure4_matrices_bin.mat');

%% PV neighbor 

%days_included in 1 RM ANOVA analysis
day_range = [1:5,9];

%number of time points
nb_time_points = numel(day_range);
%number of animals
nb_animals = 3;

% input to function as animals x days
%A
[stats.neighbor.PV.RMlme.A] = one_way_RM_lme(PV_neighbor_A(:,day_range),nb_animals, nb_time_points);

%B
[stats.neighbor.PV.RMlme.B] = one_way_RM_lme(PV_neighbor_B(:,day_range),nb_animals, nb_time_points);

%A Day 1vs2 Vs. Day 3vs4
stats.neighbor.PV.ttest.A.d12Vs34 = paired_ttest(PV_neighbor_A(:,1), PV_neighbor_A(:,3));

%A Day 1vs2 Vs. Day 5vs6
stats.neighbor.PV.ttest.A.d12Vs56 = paired_ttest(PV_neighbor_A(:,1), PV_neighbor_A(:,5));

%B Day 1vs2 Vs. Day 3vs4
stats.neighbor.PV.ttest.B.d12Vs34 = paired_ttest(PV_neighbor_B(:,1), PV_neighbor_B(:,3));

%B Day 1vs2 Vs. Day 5vs6
stats.neighbor.PV.ttest.B.d12Vs56 = paired_ttest(PV_neighbor_B(:,1), PV_neighbor_B(:,5));

%% TC neighbor SI

% input to function as animal x days
%A
[stats.neighbor.TC.si.RMlme.A] = one_way_RM_lme(TC_neighbor_A_si(:,day_range),nb_animals, nb_time_points);

%B
[stats.neighbor.TC.si.RMlme.B] = one_way_RM_lme(TC_neighbor_B_si(:,day_range),nb_animals, nb_time_points);

%A Day 1vs2 Vs. Day 3vs4
stats.neighbor.TC.si.ttest.A.d12Vs34 = paired_ttest(TC_neighbor_A_si(:,1), TC_neighbor_A_si(:,3));

%A Day 1vs2 Vs. Day 5vs6
stats.neighbor.TC.si.ttest.A.d12Vs56 = paired_ttest(TC_neighbor_A_si(:,1), TC_neighbor_A_si(:,5));

%B Day 1vs2 Vs. Day 3vs4
stats.neighbor.TC.si.ttest.B.d12Vs34 = paired_ttest(TC_neighbor_B_si(:,1), TC_neighbor_B_si(:,3));

%B Day 1vs2 Vs. Day 5vs6
stats.neighbor.TC.si.ttest.B.d12Vs56 = paired_ttest(TC_neighbor_B_si(:,1), TC_neighbor_B_si(:,5));

%% TC neighbor TS

%A
[stats.neighbor.TC.ts.RMlme.A] = one_way_RM_lme(TC_corr_d1_A(:,day_range),nb_animals, nb_time_points);

%B
[stats.neighbor.TC.ts.RMlme.B] = one_way_RM_lme(TC_corr_d1_B(:,day_range),nb_animals, nb_time_points);

%A Day 1vs2 Vs. Day 3vs4
stats.neighbor.TC.ts.ttest.A.d12Vs34 = paired_ttest(TC_corr_d1_A(:,1), TC_corr_d1_A(:,3));

%A Day 1vs2 Vs. Day 5vs6
stats.neighbor.TC.ts.ttest.A.d12Vs56 = paired_ttest(TC_corr_d1_A(:,1), TC_corr_d1_A(:,5));

%A Day 1vs2 Vs. Day 3vs4
stats.neighbor.TC.ts.ttest.B.d12Vs34 = paired_ttest(TC_corr_d1_B(:,1), TC_corr_d1_B(:,3));

%B Day 1vs2 Vs. Day 5vs6
stats.neighbor.TC.ts.ttest.B.d12Vs56 = paired_ttest(TC_corr_d1_B(:,1), TC_corr_d1_B(:,5));

%% t tests for binned data with adjustment for mult. comparisons

%PV A  - D1
%Paired t-test Recall vs Rev
stats.binned.PV.ttest.A.rec_rev.d1 = paired_ttest(binned_data.PV_corr_d1_A(:,2), binned_data.PV_A_wrRev(:,2));

%Recall vs Rev Recall
stats.binned.PV.ttest.A.rec_revRec.d1 = paired_ttest(binned_data.PV_corr_d1_A(:,2), binned_data.rev_recall_pv_A(:,2));

%PV B  - D1
%Paired t-test Recall vs Rev
stats.binned.PV.ttest.B.rec_rev.d1 = paired_ttest(binned_data.PV_corr_d1_B(:,2), binned_data.PV_B_wrRev(:,2));

%Recall vs Rev Recall
stats.binned.PV.ttest.B.rec_revRec.d1 = paired_ttest(binned_data.PV_corr_d1_B(:,2), binned_data.rev_recall_pv_B(:,2));


%PV A - D2 - recall vs reversal
stats.binned.PV.ttest.A.rec_rev.d2 = paired_ttest(binned_data.PV_corr_d1_A(:,3), binned_data.PV_A_wrRev(:,3));

%PV B - D2 - recall vs reversal
stats.binned.PV.ttest.B.rec_rev.d2 = paired_ttest(binned_data.PV_corr_d1_B(:,3), binned_data.PV_B_wrRev(:,3));


%TC A - D1 
%Paired t-test Recall vs Rev
stats.binned.TC.ttest.A.rec_rev.d1 = paired_ttest(binned_data.TC_corr_d1_A(:,2), binned_data.TC_A_wrRev(:,2));

%Recall vs Rev Recall
stats.binned.TC.ttest.A.rec_revRec.d1 = paired_ttest(binned_data.TC_corr_d1_A(:,2), binned_data.rev_recall_tc_A(:,2));

%TC B - D1 
%Paired t-test Recall vs Rev
stats.binned.TC.ttest.B.rec_rev.d1 = paired_ttest(binned_data.TC_corr_d1_B(:,2), binned_data.TC_B_wrRev(:,2));

%Recall vs Rev Recall
stats.binned.TC.ttest.B.rec_revRec.d1 = paired_ttest(binned_data.TC_corr_d1_B(:,2), binned_data.rev_recall_tc_B(:,2));

%TC A - D2
%Paired t-test Recall vs Rev
stats.binned.TC.ttest.A.rec_rev.d2 = paired_ttest(binned_data.TC_corr_d1_A(:,3), binned_data.TC_A_wrRev(:,3));

%TC B - D2
%Paired t-test Recall vs Rev
stats.binned.TC.ttest.B.rec_rev.d2 = paired_ttest(binned_data.TC_corr_d1_B(:,3), binned_data.TC_B_wrRev(:,3));

%% Adjustments to multiple comparisons

%% Output tables for Excel export
% Neighbor table for PV/TC

%PV 
%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.neighbor.PV_A] = one_way_lme_table_entry(4,'sr','by animal',...
                'PV correlation to start of each phase - reversal - A (stable learn -> reversal --> rev recall)',...
                nb_animals,stats.neighbor.PV.RMlme.A);

[t_1_rm_lme.neighbor.PV_B] = one_way_lme_table_entry(4,'sr','by animal',...
                'PV correlation to start of each phase - reversal - B (stable learn -> reversal --> rev recall)',...
                nb_animals,stats.neighbor.PV.RMlme.B);

%TC TS
[t_1_rm_lme.neighbor.TC_ts_A] = one_way_lme_table_entry(4,'sr','by animal',...
                'TC TS correlation relative neighboring days A (stable learn -> reversal)',...
                nb_animals,stats.neighbor.TC.ts.RMlme.A);

[t_1_rm_lme.neighbor.TC_ts_B] = one_way_lme_table_entry(4,'sr','by animal',...
                'TC TS correlation relative neighboring days B (stable learn -> reversal)',...
                nb_animals,stats.neighbor.TC.ts.RMlme.B);

%% T-test output tables

%Figure, Subfigure, Data aggregation, Comparison, N, Test, Degrees of Freedom, 
%Test statistic, p-value, p-value adjusted, ad. method, Significance

%3-way Holm-Sidak correction for AUC RUN comparison
%significance level
alpha = 0.05;
%# of comparisons
c = 4;

fig_num = repmat(4,4,1);
fig_sub = repmat('sr',4,1);
data_agg = repmat('by animal',4,1);
comp_descrip = {'Neighboring day PV corr score 1vs2 Vs. 3vs4 - A';...
                'Neighboring day PV corr score 1vs2 Vs. 3vs4 - B';...
                'Neighboring day TC TS corr score 1vs2 Vs. 3vs4 - A';...
                'Neighboring day TC TS corr score 1vs2 Vs. 3vs4 - B'};
n_sample = [stats.neighbor.PV.ttest.A.d12Vs34(3), stats.neighbor.PV.ttest.B.d12Vs34(3),...
    stats.neighbor.TC.ts.ttest.A.d12Vs34(3), stats.neighbor.TC.ts.ttest.B.d12Vs34(3)]';
test_name = repmat('Paired t-test',4,1);
n_dof = [stats.neighbor.PV.ttest.A.d12Vs34(4), stats.neighbor.PV.ttest.B.d12Vs34(4),...
    stats.neighbor.TC.ts.ttest.A.d12Vs34(4), stats.neighbor.TC.ts.ttest.B.d12Vs34(4)]';
test_statistic = [stats.neighbor.PV.ttest.A.d12Vs34(2), stats.neighbor.PV.ttest.B.d12Vs34(2),...
    stats.neighbor.TC.ts.ttest.A.d12Vs34(2), stats.neighbor.TC.ts.ttest.B.d12Vs34(2)]';
adj_method = repmat('N/A', 4,1);
p = [stats.neighbor.PV.ttest.A.d12Vs34(1), stats.neighbor.PV.ttest.B.d12Vs34(1),...
    stats.neighbor.TC.ts.ttest.A.d12Vs34(1), stats.neighbor.TC.ts.ttest.B.d12Vs34(1)]';
p_adj = holm_sidak_p_adj(p',c,alpha);
%replace adj value with nan
p_adj = repmat('N/A', 4,1)';

sig_level = check_p_value_sig(p);
%
%sig_level = check_p_value_sig(p_adj);

%create RUN AUC/min table
t_paired_neighbor_all = table(fig_num, fig_sub, data_agg, comp_descrip, n_sample,...
            test_name, n_dof, test_statistic,p,p_adj', adj_method, sig_level,...
            'VariableNames',{'Figure','Subfigure','Data aggregation',...
            'Comparison','N', 'Test', 'Degrees of Freedom', 'Test statistic',...
            'p-value', 'p-value adjusted', 'Adjustment method','Significance'});

%% T-test comparison output for PV and TC relative to D1

%PV D1 A Recall vs Rev
stats.binned.PV.ttest.A.rec_rev.d1  
%PV D1 A Recall vs. Rev recall
stats.binned.PV.ttest.A.rec_revRec.d1  
%PV D2 A Recall vs Rev
stats.binned.PV.ttest.A.rec_rev.d2 

%PV D1 B Recall vs Rev
stats.binned.PV.ttest.B.rec_rev.d1 
%PV D1 B Recall vs. Rev recall
stats.binned.PV.ttest.B.rec_revRec.d1
%PV D2 B Recall vs Rev
stats.binned.PV.ttest.B.rec_rev.d2

%TC D1 A Recall vs Rev
stats.binned.TC.ttest.A.rec_rev.d1
%TC D1 A Recall vs. Rev recall
stats.binned.TC.ttest.A.rec_revRec.d1 
%TC D2 A Recall vs Rev
stats.binned.TC.ttest.A.rec_rev.d2 

%TC D1 B Recall vs Rev
stats.binned.TC.ttest.B.rec_rev.d1
%TC D1 B Recall vs. Rev recall
stats.binned.TC.ttest.B.rec_revRec.d1
%TC D2 B Recall vs Rev
stats.binned.TC.ttest.B.rec_rev.d2

%% combine into 1 matrix in the order above
rel_learn_rev_stats = [stats.binned.PV.ttest.A.rec_rev.d1;...
stats.binned.PV.ttest.A.rec_revRec.d1;...
stats.binned.PV.ttest.A.rec_rev.d2;...

stats.binned.PV.ttest.B.rec_rev.d1;... 
stats.binned.PV.ttest.B.rec_revRec.d1;...
stats.binned.PV.ttest.B.rec_rev.d2;...

stats.binned.TC.ttest.A.rec_rev.d1;...
stats.binned.TC.ttest.A.rec_revRec.d1;...
stats.binned.TC.ttest.A.rec_rev.d2;... 

stats.binned.TC.ttest.B.rec_rev.d1;...
stats.binned.TC.ttest.B.rec_revRec.d1;...
stats.binned.TC.ttest.B.rec_rev.d2];


%% Write to single table - correlation relative to day 1 of learn task phase

%3-way Holm-Sidak correction for AUC RUN comparison
%significance level
alpha = 0.05;
%# of comparisons
c = 4;

fig_num = repmat(4,12,1);
fig_sub = repmat('sr',12,1);
data_agg = repmat('by animal',12,1);
comp_descrip = {'PV correlation rel D1 task phase - A - Recall vs. Reversal';...
                'PV correlation rel D1 task phase - A - Recall vs. Reversal Recall';...
                'PV correlation rel D2 task phase - A - Recall vs. Reversal';...
                'PV correlation rel D1 task phase - B - Recall vs. Reversal';...
                'PV correlation rel D1 task phase - B - Recall vs. Reversal Recall';...
                'PV correlation rel D2 task phase - B - Recall vs. Reversal';...
                'TC (TS) correlation rel D1 task phase - A - Recall vs. Reversal';...
                'TC (TS)correlation rel D1 task phase - A - Recall vs. Reversal Recall';...
                'TC (TS)correlation rel D2 task phase - A - Recall vs. Reversal';...
                'TC (TS)correlation rel D1 task phase - B - Recall vs. Reversal';...
                'TC (TS)correlation rel D1 task phase - B - Recall vs. Reversal Recall';...
                'TC (TS)correlation rel D2 task phase - B - Recall vs. Reversal'};

n_sample = rel_learn_rev_stats(:,3);
test_name = repmat('Paired t-test',12,1);
n_dof = rel_learn_rev_stats(:,4);
test_statistic = rel_learn_rev_stats(:,2);
adj_method = repmat('N/A', 12,1);
p = rel_learn_rev_stats(:,1);
%p_adj = holm_sidak_p_adj(p',c,alpha);
%replace adj value with nan
p_adj = repmat('N/A', 12,1)';

sig_level = check_p_value_sig(p);
%
%sig_level = check_p_value_sig(p_adj);

%create RUN AUC/min table
t_paired_rel_phase_all = table(fig_num, fig_sub, data_agg, comp_descrip, n_sample,...
            test_name, n_dof, test_statistic,p,p_adj', adj_method, sig_level,...
            'VariableNames',{'Figure','Subfigure','Data aggregation',...
            'Comparison','N', 'Test', 'Degrees of Freedom', 'Test statistic',...
            'p-value', 'p-value adjusted', 'Adjustment method','Significance'});

%% Output to Excel spreadsheet
t1 = repmat({' '},1,12);

%spreadsheet name
spreadsheet_name = 'reversal_learning_stats.xlsx';

%sheet name
sheet_name = 'Figure 4 sup';

%write to Excel spreadsheet
%learn 4e
writetable(t_1_rm_lme.neighbor.PV_A  ,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','overwritesheet')
writetable(t_paired_neighbor_all(1,:),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(t_1_rm_lme.neighbor.PV_B  ,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_paired_neighbor_all(2,:),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(t_1_rm_lme.neighbor.TC_ts_A  ,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_paired_neighbor_all(3,:),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(t_1_rm_lme.neighbor.TC_ts_B  ,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_paired_neighbor_all(4,:),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(t_paired_rel_phase_all(1:3,:),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(t_paired_rel_phase_all(4:6,:),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(t_paired_rel_phase_all(7:9,:),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(t_paired_rel_phase_all(10:12,:),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%% ARCHIVE %%%%
%% 1way RM ANOVA PV - A
%{
%days_included in 1 RM ANOVA analysis
day_range = [2:6,9:11];

%number of time points
nb_time_points = numel(day_range);
%number of animals
nb_animals = 3;

%1-way RM LME analysis for learn
%input data arranged as animal x day

%A
[stats.PV.RMlme.A] = one_way_RM_lme(PV_corr_d1_A(:,day_range),nb_animals, nb_time_points);

%B
[stats.PV.RMlme.B] = one_way_RM_lme(PV_corr_d1_B(:,day_range),nb_animals, nb_time_points);

%A Day 2 vs. Day 6
stats.PV.ttest.A.d2v6 = paired_ttest(PV_corr_d1_A(:,2), PV_corr_d1_A(:,6));

%B Day 2 vs. Day 6
stats.PV.ttest.B.d2v6 = paired_ttest(PV_corr_d1_B(:,2), PV_corr_d1_B(:,6));

%% TC SI A and B

%A
[stats.TC.si.RMlme.A] = one_way_RM_lme(TC_corr_d1_A_si(:,day_range),nb_animals, nb_time_points);

%B
[stats.TC.si.RMlme.B] = one_way_RM_lme(TC_corr_d1_B_si(:,day_range),nb_animals, nb_time_points);

%A Day 2 vs. Day 6
stats.TC.si.ttest.A.d2v6 = paired_ttest(TC_corr_d1_A_si(:,2), TC_corr_d1_A_si(:,6));

%B Day 2 vs. Day 6
stats.TC.si.ttest.B.d2v6 = paired_ttest(TC_corr_d1_B_si(:,2), TC_corr_d1_B_si(:,6));


%% TC TS A and B

%A
[stats.TC.ts.RMlme.A] = one_way_RM_lme(TC_corr_d1_A(:,day_range),nb_animals, nb_time_points);

%B
[stats.TC.ts.RMlme.B] = one_way_RM_lme(TC_corr_d1_B(:,day_range),nb_animals, nb_time_points);

%A Day 2 vs. Day 6
stats.TC.ts.ttest.A.d2v6 = paired_ttest(TC_corr_d1_A(:,2), TC_corr_d1_A(:,6));

%B Day 2 vs. Day 6
stats.TC.ts.ttest.B.d2v6 = paired_ttest(TC_corr_d1_B(:,2), TC_corr_d1_B(:,6));


%% 
%TC SI
[t_1_rm_lme.neighbor.TC_si_A] = one_way_lme_table_entry(4,'sr','by animal',...
                'TC SI correlation relative neighboring days A (stable learn -> reversal)',...
                nb_animals,stats.neighbor.TC.si.RMlme.A);

[t_1_rm_lme.neighbor.TC_si_B] = one_way_lme_table_entry(4,'sr','by animal',...
                'TC SI correlation relative neighboring days B (stable learn -> reversal)',...
                nb_animals,stats.neighbor.TC.si.RMlme.B);

%% Code below does auto-adjustments of multi-way comparisons

data_input = [stats.PV.ttest.A.d2v6; stats.PV.ttest.B.d2v6];

%what comparisons are being made
comp_descrip_in = {'Pv correlation rel. day 1 A D1 vs. D6';...
                'A tuned fraction across time - SI - Recall D1 vs. D7'};

[t_ttest.PV.A] = paired_ttest_table_entry(data_input,...
        9, 'a', 'by animal', comp_descrip_in);

%learn B D6/D7
data_input = [stats.neighbor.PV.ttest.A.d12Vs34; stats.neighbor.PV.ttest.B.d12Vs34];

%what comparisons are being made
comp_descrip_in = {'B tuned fraction across time - SI - Recall D1 vs. D6';...
                'B tuned fraction across time - SI - Recall D1 vs. D7'};

[x.t_ttest.recallB_6_7] = paired_ttest_table_entry(data_input,...
        9, 'a', 'by animal', comp_descrip_in);

%recall AB D6/D7
data_input = [t_stats.recallAB6; t_stats.recallAB7];

%what comparisons are being made
comp_descrip_in = {'A&B tuned fraction across time - SI - Recall D1 vs. D6';...
                'A&B tuned fraction across time - SI - Recall D1 vs. D7'};

[t_ttest.recallAB_6_7] = paired_ttest_table_entry(data_input,...
        9, 'a', 'by animal', comp_descrip_in);

%paired t-test 2 vs. 6
%}
