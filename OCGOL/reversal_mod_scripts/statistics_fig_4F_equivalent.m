%% Load data

load('D:\OCGOL_reversal\figure4_matrices.mat')

%% 1way RM ANOVA PV - A

%days_included in 1 RM ANOVA analysis
day_range = [2:6,9:11];

%number of time points
nb_time_points = numel(day_range);
%number of animals
nb_animals = 3;

% input to function as days x animals
%A
[stats.PV.RMlme.A] = one_way_RM_lme(PV_corr_d1_A(:,day_range)',nb_animals, nb_time_points);

%B
[stats.PV.RMlme.B] = one_way_RM_lme(PV_corr_d1_B(:,day_range)',nb_animals, nb_time_points);

%A Day 2 vs. Day 6
stats.PV.ttest.A.d2v6 = paired_ttest(PV_corr_d1_A(:,2), PV_corr_d1_A(:,6));

%B Day 2 vs. Day 6
stats.PV.ttest.B.d2v6 = paired_ttest(PV_corr_d1_B(:,2), PV_corr_d1_B(:,6));

%% TC SI A and B

%A
[stats.TC.si.RMlme.A] = one_way_RM_lme(TC_corr_d1_A_si(:,day_range)',nb_animals, nb_time_points);

%B
[stats.TC.si.RMlme.B] = one_way_RM_lme(TC_corr_d1_B_si(:,day_range)',nb_animals, nb_time_points);

%A Day 2 vs. Day 6
stats.TC.si.ttest.A.d2v6 = paired_ttest(TC_corr_d1_A_si(:,2), TC_corr_d1_A_si(:,6));

%B Day 2 vs. Day 6
stats.TC.si.ttest.B.d2v6 = paired_ttest(TC_corr_d1_B_si(:,2), TC_corr_d1_B_si(:,6));


%% TC TS A and B

%A
[stats.TC.ts.RMlme.A] = one_way_RM_lme(TC_corr_d1_A(:,day_range)',nb_animals, nb_time_points);

%B
[stats.TC.ts.RMlme.B] = one_way_RM_lme(TC_corr_d1_B(:,day_range)',nb_animals, nb_time_points);

%A Day 2 vs. Day 6
stats.TC.ts.ttest.A.d2v6 = paired_ttest(TC_corr_d1_A(:,2), TC_corr_d1_A(:,6));

%B Day 2 vs. Day 6
stats.TC.ts.ttest.B.d2v6 = paired_ttest(TC_corr_d1_B(:,2), TC_corr_d1_B(:,6));

%% PV neighbor 

%days_included in 1 RM ANOVA analysis
day_range = [1:5,9];

%number of time points
nb_time_points = numel(day_range);
%number of animals
nb_animals = 3;

% input to function as days x animals
%A
[stats.neighbor.PV.RMlme.A] = one_way_RM_lme(PV_neighbor_A(:,day_range)',nb_animals, nb_time_points);

%B
[stats.neighbor.PV.RMlme.B] = one_way_RM_lme(PV_neighbor_B(:,day_range)',nb_animals, nb_time_points);

%A Day 1vs2 Vs. Day 5vs6
stats.neighbor.PV.ttest.A.d12Vs56 = paired_ttest(PV_neighbor_A(:,1), PV_neighbor_A(:,5));

%B Day 1vs2 Vs. Day 5vs6
stats.neighbor.PV.ttest.B.d12Vs56 = paired_ttest(PV_neighbor_B(:,1), PV_neighbor_B(:,5));

%% TC neighbor SI

% input to function as days x animals
%A
[stats.neighbor.TC.si.RMlme.A] = one_way_RM_lme(TC_corr_d1_A_si(:,day_range)',nb_animals, nb_time_points);

%B
[stats.neighbor.TC.si.RMlme.B] = one_way_RM_lme(TC_corr_d1_B_si(:,day_range)',nb_animals, nb_time_points);

%A Day 1vs2 Vs. Day 5vs6
stats.neighbor.TC.si.ttest.A.d12Vs56 = paired_ttest(TC_corr_d1_A_si(:,1), TC_corr_d1_A_si(:,5));

%B Day 1vs2 Vs. Day 5vs6
stats.neighbor.TC.si.ttest.B.d12Vs56 = paired_ttest(TC_corr_d1_B_si(:,1), TC_corr_d1_B_si(:,5));

%% TC neighbor TS

% input to function as days x animals
%A
[stats.neighbor.TC.ts.RMlme.A] = one_way_RM_lme(TC_corr_d1_A(:,day_range)',nb_animals, nb_time_points);

%B
[stats.neighbor.TC.ts.RMlme.B] = one_way_RM_lme(TC_corr_d1_B(:,day_range)',nb_animals, nb_time_points);

%A Day 1vs2 Vs. Day 5vs6
stats.neighbor.TC.ts.ttest.A.d12Vs56 = paired_ttest(TC_corr_d1_A(:,1), TC_corr_d1_A(:,5));

%B Day 1vs2 Vs. Day 5vs6
stats.neighbor.TC.ts.ttest.B.d12Vs56 = paired_ttest(TC_corr_d1_B(:,1), TC_corr_d1_B(:,5));



%% Output tables for Excel export

%PV 
%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.PV_A] = one_way_lme_table_entry(4,'sr','by animal',...
                'PV correlation relative for day 1 A (stable learn -> reversal)',...
                nb_animals,stats.PV.RMlme.A);

[t_1_rm_lme.PV_B] = one_way_lme_table_entry(4,'sr','by animal',...
                'PV correlation relative for day 1 B (stable learn -> reversal)',...
                nb_animals,stats.PV.RMlme.B);

%TC SI
[t_1_rm_lme.TC_si_A] = one_way_lme_table_entry(4,'sr','by animal',...
                'TC SI correlation relative for day 1 A (stable learn -> reversal)',...
                nb_animals,stats.TC.si.RMlme.A);

[t_1_rm_lme.TC_si_B] = one_way_lme_table_entry(4,'sr','by animal',...
                'TC SI correlation relative for day 1 B (stable learn -> reversal)',...
                nb_animals,stats.TC.si.RMlme.B);

%TC TS
[t_1_rm_lme.TC_ts_A] = one_way_lme_table_entry(4,'sr','by animal',...
                'TC TS correlation relative for day 1 A (stable learn -> reversal)',...
                nb_animals,stats.TC.ts.RMlme.A);

[t_1_rm_lme.TC_ts_B] = one_way_lme_table_entry(4,'sr','by animal',...
                'TC TS correlation relative for day 1 B (stable learn -> reversal)',...
                nb_animals,stats.TC.ts.RMlme.B);

%% CONTINUE WITH NEIGHBOR TABLE OUTPUT HERE


%%

%recall A D6/D7
data_input = [t_stats.recallA6; t_stats.recallA7];

%what comparisons are being made
comp_descrip_in = {'A tuned fraction across time - SI - Recall D1 vs. D6';...
                'A tuned fraction across time - SI - Recall D1 vs. D7'};

[t_ttest.recallA_6_7] = paired_ttest_table_entry(data_input,...
        9, 'a', 'by animal', comp_descrip_in);

%learn B D6/D7
data_input = [t_stats.recallB6; t_stats.recallB7];

%what comparisons are being made
comp_descrip_in = {'B tuned fraction across time - SI - Recall D1 vs. D6';...
                'B tuned fraction across time - SI - Recall D1 vs. D7'};

[t_ttest.recallB_6_7] = paired_ttest_table_entry(data_input,...
        9, 'a', 'by animal', comp_descrip_in);

%recall AB D6/D7
data_input = [t_stats.recallAB6; t_stats.recallAB7];

%what comparisons are being made
comp_descrip_in = {'A&B tuned fraction across time - SI - Recall D1 vs. D6';...
                'A&B tuned fraction across time - SI - Recall D1 vs. D7'};

[t_ttest.recallAB_6_7] = paired_ttest_table_entry(data_input,...
        9, 'a', 'by animal', comp_descrip_in);

%paired t-test 2 vs. 6
