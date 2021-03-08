%% Export stats for all supplementary figures in paper

%import supplementary data

%figure 2/3 sup data
%load('G:\Google_drive\task_selective_place_paper\matlab_data\source_data_sup_2_3.mat');

%figure 4/5 sup data
%load('G:\Google_drive\task_selective_place_paper\matlab_data\source_data_fig4_5_and_sup.mat');

%ex figure 10 sup data
%load('G:\Google_drive\task_selective_place_paper\matlab_data\source_data_ex10_sup.mat');

%laptop file directory
load('C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data\source_data_sup_2_3.mat');
load('C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data\source_data_fig4_5_and_sup.mat');
load('C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data\source_data_ex10_sup.mat');

%% Figure 4 A-B speed difference for A- and B- selective neurons

%speed difference for A and B- selective neurons
speed_diff_Asel = source_data_sup_2_3.event_speed_plot.Asel_speed_diff;
speed_diff_Bsel = source_data_sup_2_3.event_speed_plot.Bsel_speed_diff;

%remove any nan values from each vector
%A
nan_ct = find(isnan(speed_diff_Asel)== 1);
 if ~isempty(nan_ct)
     speed_diff_Asel(nan_ct) = [];
 end
%B
nan_ct = find(isnan(speed_diff_Bsel)== 1);
 if ~isempty(nan_ct)
     speed_diff_Bsel(nan_ct) = [];
 end
 
%for 1 sample Wilcoxon tests
paired_wilcoxon_stats_Asel = paired_wilcoxon_signrank(speed_diff_Asel,0);

%test comparison description for each table entry \
comp_descrip_in = {'A-B speed difference for A selective neurons'};

%table entries for paired wilcoxon test learning
%A selective
t_paired_wilcox_Asel = paired_wilcoxon_table_entry_no_adj(4,'b','pooled',...
                comp_descrip_in,paired_wilcoxon_stats_Asel);

%for 1 sample Wilcoxon tests
paired_wilcoxon_stats_Bsel = paired_wilcoxon_signrank(speed_diff_Bsel,0);

%test comparison description for each table entry \
comp_descrip_in = {'A-B speed difference for B selective neurons'};            
            
%B selective
t_paired_wilcox_Bsel = paired_wilcoxon_table_entry_no_adj(4,'b','pooled',...
                comp_descrip_in,paired_wilcoxon_stats_Bsel);
            
%% Figure 4 Place field count differences between task selective place cells

%normalized counts
Asel_norm = source_data_sup_2_3.pf_count.Asel_norm;
Bsel_norm = source_data_sup_2_3.pf_count.Bsel_norm;

AB.A_norm = source_data_sup_2_3.pf_count.AB.A_norm;
AB.B_norm = source_data_sup_2_3.pf_count.AB.B_norm;

%assemble into single field matrices (animal x A,B,AB-A, AB-B
single_field = [Asel_norm(:,1), Bsel_norm(:,1), AB.A_norm(:,1), AB.B_norm(:,1)];
double_field = [Asel_norm(:,2), Bsel_norm(:,2), AB.A_norm(:,2), AB.B_norm(:,2)];
triple_field = [Asel_norm(:,3), Bsel_norm(:,3), AB.A_norm(:,3), AB.B_norm(:,3)];

%single paired tests
%A vs B
paired_wilcoxon_stats.single.AvsB = paired_wilcoxon_signrank(single_field(:,1),single_field(:,2));
%A vs AB-A
paired_wilcoxon_stats.single.AvsAB_A = paired_wilcoxon_signrank(single_field(:,1),single_field(:,3));
%B vs AB-B
paired_wilcoxon_stats.single.BvsAB_B = paired_wilcoxon_signrank(single_field(:,2),single_field(:,4));

%double paired tests
paired_wilcoxon_stats.double.AvsB = paired_wilcoxon_signrank(double_field(:,1),double_field(:,2));
%A vs AB-A
paired_wilcoxon_stats.double.AvsAB_A = paired_wilcoxon_signrank(double_field(:,1),double_field(:,3));
%B vs AB-B
paired_wilcoxon_stats.double.BvsAB_B = paired_wilcoxon_signrank(double_field(:,2),double_field(:,4));

%triple paired tests
paired_wilcoxon_stats.triple.AvsB = paired_wilcoxon_signrank(triple_field(:,1),triple_field(:,2));
%A vs AB-A
paired_wilcoxon_stats.triple.AvsAB_A = paired_wilcoxon_signrank(triple_field(:,1),triple_field(:,3));
%B vs AB-B
paired_wilcoxon_stats.triple.BvsAB_B = paired_wilcoxon_signrank(triple_field(:,2),triple_field(:,4));

%single field stats combined
stats_in_single = [paired_wilcoxon_stats.single.AvsB; paired_wilcoxon_stats.single.AvsAB_A; paired_wilcoxon_stats.single.BvsAB_B];

comp_descrip = {'A vs B (1 vs 2) single field';'A vs A&B A (1 vs 3) single field';'B vs A&B (2 vs 4) single field'};

%do the adjustment in the function  for multo comp (3-way Holm-Sidak
%correction)
[t_out.wilcox_singlePF ] = paired_wilcoxon_table_entry(4,'c','by animal',...
                comp_descrip, stats_in_single);

%double field stats combined
stats_in_double = [paired_wilcoxon_stats.double.AvsB; paired_wilcoxon_stats.double.AvsAB_A; paired_wilcoxon_stats.double.BvsAB_B];

comp_descrip = {'A vs B (1 vs 2) double field';'A vs A&B A (1 vs 3) double field';'B vs A&B (2 vs 4) double field'};

%do the adjustment in the function  for multo comp (3-way Holm-Sidak
%correction)
[t_out.wilcox_doublePF ] = paired_wilcoxon_table_entry(4,'c','by animal',...
                comp_descrip, stats_in_double);

%triple field stats combined
stats_in_triple = [paired_wilcoxon_stats.triple.AvsB; paired_wilcoxon_stats.triple.AvsAB_A; paired_wilcoxon_stats.triple.BvsAB_B];

comp_descrip = {'A vs B (1 vs 2) triple field';'A vs A&B A (1 vs 3) triple field';'B vs A&B (2 vs 4) triple field'};

%do the adjustment in the function  for multo comp (3-way Holm-Sidak
%correction)
[t_out.wilcox_triplePF ] = paired_wilcoxon_table_entry(4,'c','by animal',...
                comp_descrip, stats_in_triple);
%% Figure 4 Place field width differences between task selective place cells
%2-sample KS test

%unload the data
width_Asel = source_data_sup_2_3.pf_width_pool.Asel;
width_Bsel = source_data_sup_2_3.pf_width_pool.Bsel;
width_ABA = source_data_sup_2_3.pf_width_pool.AB.A;
width_ABB = source_data_sup_2_3.pf_width_pool.AB.B;

%run all the KS tests

%Asel vs Bsel
kstest2_AB = kstest2_stats(width_Asel,width_Bsel);

%Asel vs ABAsel
kstest2_A_ABA = kstest2_stats(width_Asel,width_ABA);

%Asel vs ABBsel
kstest2_A_ABB = kstest2_stats(width_Asel,width_ABB);

%Bsel vs ABAsel
kstest2_B_ABA = kstest2_stats(width_Bsel,width_ABA);

%Bsel vs ABBsel
kstest2_B_ABB = kstest2_stats(width_Bsel,width_ABB);

%ABAsel vs ABBsel
kstest2_ABA_ABB = kstest2_stats(width_ABA,width_ABB);

%KS test mult compare correction and table entry generation
%setup data input for KS2 table function
input_data = [kstest2_AB; kstest2_A_ABA; kstest2_A_ABB;...
            kstest2_B_ABA; kstest2_B_ABB; kstest2_ABA_ABB];

%description of each KS test comparison
comp_descrip = {'A sel PF width vs. Bsel PF width (A vs. B)';...
    'A sel PF width vs AB.A sel PF width (A vs. C)';...
    'A sel PF width vs. AB.B sel PF width (A vs. D)';...
    'Bsel PF width vs. AB.A sel PF width (B vs. C)';...
    'Bsel PF width vs. AB.B sel PF width (B vs. D)';...
    'AB.A sel PF width vs. AB.B sel PF width (C vs. D)'};        

%generate table output for 2KS test with HS correction 
t_ks2_pf_width = kstest2_mult_compare(4,'d', 'pooled', comp_descrip, input_data);

%% Figure 7 Activity index score difference and activity rate diff. for activity remapping neurons

%input data 
com_idx_score = source_data_sup_2_3.remap_sup_plot_data.common_idx_values;
activity_idx_score = source_data_sup_2_3.remap_sup_plot_data.remap_idx_values;

%run KS test
kstest_com_vs_activity = kstest2_stats(com_idx_score,activity_idx_score);

comp_descrip = {'Common vs. global remapping activity index score'};
%single KS entry into stats table (not single entry, but no p val adj for
%mult comp. 
t_kstest_com_vs_act = kstest2_single_entry(7,'a','pooled',...
                comp_descrip, kstest_com_vs_activity);
            
%activity rate between activity remappers on A vs. B laps
AUC_min_act_remap = source_data_sup_2_3.remap_sup_plot_data.auc_min_activity_remap_rate;
%run the test
paired_wilcoxon_stats_act_remap_AUC = paired_wilcoxon_signrank(AUC_min_act_remap(1,:),AUC_min_act_remap(2,:));

%test comparison description for each table entry \
comp_descrip_in = {'Activity rate in run epochs for A vs B laps in all actvity remapping neurons'};

%table entries for paired wilcoxon test learning
t_paired_wilcox_activity_rate_com_activity = paired_wilcoxon_table_entry_no_adj(7,'c','pooled',...
                comp_descrip_in,paired_wilcoxon_stats_act_remap_AUC);


%% Figure 9 eqivalent long term correlation data (Fig 4) fraction tuned across time using SI criterion

%unload data (animal x day)
learn.A = source_data_short_learn_recall.si_frac_export.learning.A;
learn.B = source_data_short_learn_recall.si_frac_export.learning.B;
learn.AB = source_data_short_learn_recall.si_frac_export.learning.AB;

%1-way RM LME analysis for learn
%input data arranged as animal x day

%A
%number of time points
nb_time_points = 7;
%number of animals
nb_animals = 6;
[lme_stats.learnA] = one_way_RM_lme(learn.A,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.learnA] = one_way_lme_table_entry(9,'a','by animal',...
                'A tuned fraction across time - Learning D1-D7 - SI',...
                size(learn.A,1),lme_stats.learnA);
%B
[lme_stats.learnB] = one_way_RM_lme(learn.B,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.learnB] = one_way_lme_table_entry(9,'a','by animal',...
                'B tuned fraction across time - Learning D1-D7 - SI',...
                size(learn.B,1),lme_stats.learnB);

%AB
[lme_stats.learnAB] = one_way_RM_lme(learn.AB,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.learnAB] = one_way_lme_table_entry(9,'a','by animal',...
                'A&B tuned fraction across time - Learning D1-D7 - SI',...
                size(learn.AB,1),lme_stats.learnAB);

%run paired t-tests here for each class of neurons
%Learn TS A
%d1 vs d6 test
t_stats.learnA6 = paired_ttest(learn.A(:,1), learn.A(:,6));
%d1 vs d7 test
t_stats.learnA7 = paired_ttest(learn.A(:,1), learn.A(:,7));

%Learn TS B
%d1 vs d6 test
t_stats.learnB6 = paired_ttest(learn.B(:,1), learn.B(:,6));
%d1 vs d7 test
t_stats.learnB7 = paired_ttest(learn.B(:,1), learn.B(:,7));

%Learn TS AB
%d1 vs d6 test
t_stats.learnAB6 = paired_ttest(learn.AB(:,1), learn.AB(:,6));
%d1 vs d7 test
t_stats.learnAB7 = paired_ttest(learn.AB(:,1), learn.AB(:,7));

%learn A D6/D7
data_input = [t_stats.learnA6; t_stats.learnA7];

%what comparisons are being made
comp_descrip_in = {'A tuned fraction across time - SI - Learning D1 vs. D6';...
                'A tuned fraction across time - SI - Learning D1 vs. D7'};

[t_ttest.learnA_6_7] = paired_ttest_table_entry(data_input,...
        9, 'a', 'by animal', comp_descrip_in);

%learn B D6/D7
data_input = [t_stats.learnB6; t_stats.learnB7];

%what comparisons are being made
comp_descrip_in = {'B tuned fraction across time - SI - Learning D1 vs. D6';...
                'B tuned fraction across time - SI - Learning D1 vs. D7'};

[t_ttest.learnB_6_7] = paired_ttest_table_entry(data_input,...
        9, 'a', 'by animal', comp_descrip_in);

%learn AB D6/D7
data_input = [t_stats.learnAB6; t_stats.learnAB7];

%what comparisons are being made
comp_descrip_in = {'A&B tuned fraction across time - SI - Learning D1 vs. D6';...
                'A&B tuned fraction across time - I - Learning D1 vs. D7'};

[t_ttest.learnAB_6_7] = paired_ttest_table_entry(data_input,...
        9, 'a', 'by animal', comp_descrip_in);

%% Figure 9a SI Recall analysis - fraction tuned across time

%short recall data (animal x day)
recall.A = source_data_short_learn_recall.si_frac_export.recall.A;
recall.B = source_data_short_learn_recall.si_frac_export.recall.B;
recall.AB = source_data_short_learn_recall.si_frac_export.recall.AB;

%A
%number of time points
nb_time_points = 7;
%number of animals
nb_animals = 5;
[lme_stats.recallA] = one_way_RM_lme(recall.A,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.recallA] = one_way_lme_table_entry(9,'a','by animal',...
                'A tuned fraction across time - Recall D1-D9 - SI',...
                size(recall.A,1),lme_stats.recallA);
%B
[lme_stats.recallB] = one_way_RM_lme(recall.B,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.recallB] = one_way_lme_table_entry(9,'a','by animal',...
                'B tuned fraction across time - Recall D1-D9 - SI',...
                size(recall.B,1),lme_stats.recallB);

%AB
[lme_stats.recallAB] = one_way_RM_lme(recall.AB,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.recallAB] = one_way_lme_table_entry(9,'a','by animal',...
                'A&B tuned fraction across time - Recall D1-D9 - SI',...
                size(recall.AB,1),lme_stats.recallAB);

%run paired t-tests here for each class of neurons
%Recall TS A 
%d1 vs d6 test
t_stats.recallA6 = paired_ttest(recall.A(:,1), recall.A(:,4));
%d1 vs d7 test
t_stats.recallA7 = paired_ttest(recall.A(:,1), recall.A(:,5));

%Recall TS B
%d1 vs d6 test
t_stats.recallB6 = paired_ttest(recall.B(:,1), recall.B(:,4));
%d1 vs d7 test
t_stats.recallB7 = paired_ttest(recall.B(:,1), recall.B(:,5));

%Recall TS AB
%d1 vs d6 test
t_stats.recallAB6 = paired_ttest(recall.AB(:,1), recall.AB(:,4));
%d1 vs d7 test
t_stats.recallAB7 = paired_ttest(recall.AB(:,1), recall.AB(:,5));

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

%% TC S.I. Correlation A

TC_si_learn.A = source_data_short_learn_recall.TC.si.st_learn.A';
TC_si_recall.A = source_data_short_learn_recall.TC.si.st_recall.d4_d5_sub.A';

%format data for 2-way analysis
TC_si_recall.A(:,[4,5]) = nan;

%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = TC_si_learn.A(:,[2,3,6,7]);
data_in_2 = TC_si_recall.A(:,[2,3,6,7]);

%output stats for 2way RM LME analysis
lme_stats.TC_si_learn_recall.A = two_way_rm_lme(data_in_1,data_in_2);

%create formatting table for 2-way analysis
[t_2way_rm_lme.TC_si_learn_recallA] = two_way_lme_table_entry(9,'b','by animal',...
                'Tuning curve SI correlation relative to D1  - A laps',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.TC_si_learn_recall.A);       
            
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%PV day 6 learn vs recall
[ttest_stats.TC_si_learn_recall_6A] = unpaired_ttest(TC_si_learn.A(:,6), TC_si_recall.A(:,6));

%PV day 7 learn vs recall
[ttest_stats.TC_si_learn_recall_7A] = unpaired_ttest(TC_si_learn.A(:,7), TC_si_recall.A(:,7));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.TC_si_learn_recall_6A; ttest_stats.TC_si_learn_recall_7A];

comp_descrip_in = {'Learning vs recall Tuning curve SI correlation D6 - A laps';...
                'Learning vs recall Tuning curve SI correlation D7 - A laps'};

[t_ttest.TC_si_learn_recallA_6_7] = unpaired_ttest_table_entry(data_input,...
        9, 'b', 'by animal', comp_descrip_in);

%paired t-test comparisons

%PV d2 vs. d6 time learn A 
ttest_stats.TC_si_learnA2v6 = paired_ttest(TC_si_learn.A(:,2), TC_si_learn.A(:,6));

%PV d2 vs d7 time learn A
ttest_stats.TC_si_learnA2v7 = paired_ttest(TC_si_learn.A(:,2), TC_si_learn.A(:,7));

%PV d2 vs. d6 time recall A 
ttest_stats.TC_si_recallA2v6 = paired_ttest(TC_si_recall.A(:,2), TC_si_recall.A(:,6));

%PV d2 vs d7 time recall A
ttest_stats.TC_si_recallA2v7 = paired_ttest(TC_si_recall.A(:,2), TC_si_recall.A(:,7));

%generate paired t-test table entries
%learn A time D2 vs. D7
data_input = [ttest_stats.TC_si_learnA2v6; ttest_stats.TC_si_learnA2v7];

%what comparisons are being made
comp_descrip_in = {'Learning Tuning curve SI correlation - time -D2 vs. D6 - A laps';...
                'Learning Tuning curve SI correlation - time -D2 vs. D7 - A laps'};

[t_ttest.TC_si_learnA2v6_7] = paired_ttest_table_entry(data_input,...
        9, 'b', 'by animal', comp_descrip_in);
    
%recall A time D2 vs. D7
data_input = [ttest_stats.TC_si_recallA2v6; ttest_stats.TC_si_recallA2v7];

%what comparisons are being made
comp_descrip_in = {'Recall Tuning curve SI correlation - time -D2 vs. D6 - A laps';...
                'Recall Tuning curve SI correlation - time -D2 vs. D7 - A laps'};

[t_ttest.TC_si_recallA2v6_7] = paired_ttest_table_entry(data_input,...
        9, 'b', 'by animal', comp_descrip_in);
    
%% TC S.I. Correlation B

TC_si_learn.B = source_data_short_learn_recall.TC.si.st_learn.B';
TC_si_recall.B = source_data_short_learn_recall.TC.si.st_recall.d4_d5_sub.B';

%format data for 2-way analysis
TC_si_recall.B(:,[4,5]) = nan;

%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = TC_si_learn.B(:,[2,3,6,7]);
data_in_2 = TC_si_recall.B(:,[2,3,6,7]);

%output stats for 2way RM LME analysis
lme_stats.TC_si_learn_recall.B = two_way_rm_lme(data_in_1,data_in_2);

%create formatting table for 2-way analysis
[t_2way_rm_lme.TC_si_learn_recallB] = two_way_lme_table_entry(9,'b','by animal',...
                'Tuning curve SI correlation relative to D1  - B laps',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.TC_si_learn_recall.B);       
            
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%PV day 6 learn vs recall
[ttest_stats.TC_si_learn_recall_6B] = unpaired_ttest(TC_si_learn.B(:,6), TC_si_recall.B(:,6));

%PV day 7 learn vs recall
[ttest_stats.TC_si_learn_recall_7B] = unpaired_ttest(TC_si_learn.B(:,7), TC_si_recall.B(:,7));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.TC_si_learn_recall_6B; ttest_stats.TC_si_learn_recall_7B];

comp_descrip_in = {'Learning vs recall Tuning curve SI correlation D6 - B laps';...
                'Learning vs recall Tuning curve SI correlation D7 - B laps'};

[t_ttest.TC_si_learn_recallB_6_7] = unpaired_ttest_table_entry(data_input,...
        9, 'b', 'by animal', comp_descrip_in);

%paired t-test comparisons

%PV d2 vs. d6 time learn A 
ttest_stats.TC_si_learnB2v6 = paired_ttest(TC_si_learn.B(:,2), TC_si_learn.B(:,6));

%PV d2 vs d7 time learn A
ttest_stats.TC_si_learnB2v7 = paired_ttest(TC_si_learn.B(:,2), TC_si_learn.B(:,7));

%PV d2 vs. d6 time recall A 
ttest_stats.TC_si_recallB2v6 = paired_ttest(TC_si_recall.B(:,2), TC_si_recall.B(:,6));

%PV d2 vs d7 time recall A
ttest_stats.TC_si_recallB2v7 = paired_ttest(TC_si_recall.B(:,2), TC_si_recall.B(:,7));

%generate paired t-test table entries
%learn A time D2 vs. D7
data_input = [ttest_stats.TC_si_learnB2v6; ttest_stats.TC_si_learnB2v7];

%what comparisons are being made
comp_descrip_in = {'Learning Tuning curve SI correlation - time -D2 vs. D6 - B laps';...
                'Learning Tuning curve SI correlation - time -D2 vs. D7 - B laps'};

[t_ttest.TC_si_learnB2v6_7] = paired_ttest_table_entry(data_input,...
        9, 'b', 'by animal', comp_descrip_in);
    
%recall A time D2 vs. D7
data_input = [ttest_stats.TC_si_recallB2v6; ttest_stats.TC_si_recallB2v7];

%what comparisons are being made
comp_descrip_in = {'Recall Tuning curve SI correlation - time -D2 vs. D6 - B laps';...
                'Recall Tuning curve SI correlation - time -D2 vs. D7 - B laps'};

[t_ttest.TC_si_recallB2v6_7] = paired_ttest_table_entry(data_input,...
        9, 'b', 'by animal', comp_descrip_in);    
    

%% TC SI neighboring days A

%TC SI neighbor A

TC_si_n_learn.A = source_data_short_learn_recall.TC.neighbor.si.st_learn.A;
TC_si_n_recall.A = source_data_short_learn_recall.TC.neighbor.si.st_recall.A;

%format data for 2-way analysis
%PVn_recall.A(:,[4,5]) = nan;

%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = TC_si_n_learn.A(:,[1,2,6]);
data_in_2 = TC_si_n_recall.A(:,[1,2,6]);

%output stats for 2way RM LME analysis
lme_stats.TC_si_n_learn_recall.A = two_way_rm_lme(data_in_1,data_in_2);

%%%

%create formatting table for 2-way analysis
[t_2way_rm_lme.TC_si_n_learn_recallA] = two_way_lme_table_entry(9,'c','by animal',...
                'Tuning curve (SI) correlation neighboring days - learning vs. recall - A, 1 vs. 2, 2 vs. 3, 6 vs. 7',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.TC_si_n_learn_recall.A);       

%paired t-test comparisons
%PV neighbor corr 1 vs 2 vs. 6 vs 7 time learn A 
ttest_stats.TC_si_n_learnA12v67 = paired_ttest(TC_si_n_learn.A(:,1), TC_si_n_learn.A(:,6));

%PV neighbor corr 1vs 2 vs. 6 vs 7 time recall A 
ttest_stats.TC_si_n_recallA12v67 = paired_ttest(TC_si_n_recall.A(:,1), TC_si_n_recall.A(:,6));

%paired t-test single-entry table - PV neighbor Learn 1vs2 vs 6 vs 7
data_input = [ttest_stats.TC_si_n_learnA12v67];
comp_descrip_in = {'TC correlation learning A - time - 1 vs. 2 vs. 6 vs. 7'};
[t_ttest.TC_si_n_learnA12v67] = paired_ttest_single_table_entry(data_input,...
        9, 'c', 'by animal', comp_descrip_in);

%paired t-test single-entry table - PV neighbor recall 1vs2 vs 6 vs 7
data_input = [ttest_stats.TC_si_n_recallA12v67];
comp_descrip_in = {'TC correlation recall A - time - 1 vs. 2 vs. 6 vs. 7'};
[t_ttest.TC_si_n_recallA12v67] = paired_ttest_single_table_entry(data_input,...
        9, 'c', 'by animal', comp_descrip_in);
    
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%PV neighbor correlation day 1 vs. 2 learn vs recall
[ttest_stats.TC_si_n_learn_recall_1v2A] = unpaired_ttest(TC_si_n_learn.A(:,1), TC_si_n_recall.A(:,1));

%PV  neighbor day 6 vs. 7 learn vs recall
[ttest_stats.TC_si_n_learn_recall_6v7A] = unpaired_ttest(TC_si_n_learn.A(:,6), TC_si_n_recall.A(:,6));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.TC_si_n_learn_recall_1v2A; ttest_stats.TC_si_n_learn_recall_6v7A];

comp_descrip_in = {'TC correlation learning vs. recall A - behavior - 1 vs. 2 ';...
                'TC correlation learning vs. recall A - behavior -6 vs. 7 '};

[t_ttest.TC_si_n_learn_recallA_12_67] = unpaired_ttest_table_entry(data_input,...
        9, 'c', 'by animal', comp_descrip_in);


%% TC SI neighboring days B

%TC SI neighbor B

TC_si_n_learn.B = source_data_short_learn_recall.TC.neighbor.si.st_learn.B;
TC_si_n_recall.B = source_data_short_learn_recall.TC.neighbor.si.st_recall.B;

%format data for 2-way analysis
%PVn_recall.A(:,[4,5]) = nan;

%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = TC_si_n_learn.B(:,[1,2,6]);
data_in_2 = TC_si_n_recall.B(:,[1,2,6]);

%output stats for 2way RM LME analysis
lme_stats.TC_si_n_learn_recall.B = two_way_rm_lme(data_in_1,data_in_2);

%%%

%create formatting table for 2-way analysis
[t_2way_rm_lme.TC_si_n_learn_recallB] = two_way_lme_table_entry(9,'c','by animal',...
                'Tuning curve (SI) correlation neighboring days - learning vs. recall - B, 1 vs. 2, 2 vs. 3, 6 vs. 7',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.TC_si_n_learn_recall.B);       

%paired t-test comparisons
%PV neighbor corr 1 vs 2 vs. 6 vs 7 time learn A 
ttest_stats.TC_si_n_learnB12v67 = paired_ttest(TC_si_n_learn.B(:,1), TC_si_n_learn.B(:,6));

%PV neighbor corr 1vs 2 vs. 6 vs 7 time recall A 
ttest_stats.TC_si_n_recallB12v67 = paired_ttest(TC_si_n_recall.B(:,1), TC_si_n_recall.B(:,6));

%paired t-test single-entry table - PV neighbor Learn 1vs2 vs 6 vs 7
data_input = [ttest_stats.TC_si_n_learnB12v67];
comp_descrip_in = {'TC correlation learning B - time - 1 vs. 2 vs. 6 vs. 7'};
[t_ttest.TC_si_n_learnB12v67] = paired_ttest_single_table_entry(data_input,...
        9, 'c', 'by animal', comp_descrip_in);

%paired t-test single-entry table - PV neighbor recall 1vs2 vs 6 vs 7
data_input = [ttest_stats.TC_si_n_recallB12v67];
comp_descrip_in = {'TC correlation recall B - time - 1 vs. 2 vs. 6 vs. 7'};
[t_ttest.TC_si_n_recallB12v67] = paired_ttest_single_table_entry(data_input,...
        9, 'c', 'by animal', comp_descrip_in);
    
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%PV neighbor correlation day 1 vs. 2 learn vs recall
[ttest_stats.TC_si_n_learn_recall_1v2B] = unpaired_ttest(TC_si_n_learn.B(:,1), TC_si_n_recall.B(:,1));

%PV  neighbor day 6 vs. 7 learn vs recall
[ttest_stats.TC_si_n_learn_recall_6v7B] = unpaired_ttest(TC_si_n_learn.B(:,6), TC_si_n_recall.B(:,6));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.TC_si_n_learn_recall_1v2B; ttest_stats.TC_si_n_learn_recall_6v7B];

comp_descrip_in = {'TC correlation learning vs. recall B - behavior - 1 vs. 2 ';...
                'TC correlation learning vs. recall B - behavior -6 vs. 7 '};

[t_ttest.TC_si_n_learn_recallB_12_67] = unpaired_ttest_table_entry(data_input,...
        9, 'c', 'by animal', comp_descrip_in);
    
%% Centroid SI distance difference A
%Centroid SI neighbor A

cent_si_learn.A = source_data_short_learn_recall.angle_diff.si.st_learn.animal.mean.A;
%same data as in raw sub-struct
cent_si_recall.A = source_data_short_learn_recall.angle_diff.si.st_recall.animal.mean.A;


%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = cent_si_learn.A(:,[2,3,6,7]);
data_in_2 = cent_si_recall.A(:,[2,3,6,7]);

%output stats for 2way RM LME analysis
lme_stats.cent_si_learn_recall.A = two_way_rm_lme(data_in_1,data_in_2);

%create formatting table for 2-way analysis
[t_2way_rm_lme.cent_si_learn_recallA] = two_way_lme_table_entry(9,'d','by animal',...
                'Centroid difference relative to d1 - learning vs. recall - A - rel d1, rel d2, rel d5, rel d6',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.cent_si_learn_recall.A);       

%paired t-test comparisons (learn)
%rel d1 vs d5 learning
ttest_stats.cent_si_learnA1v5 = paired_ttest(cent_si_learn.A(:,2), cent_si_learn.A(:,6));

%rel d1 vs d6 learning
ttest_stats.cent_si_learnA1v6 = paired_ttest(cent_si_learn.A(:,2), cent_si_learn.A(:,7));

%paired t-test single-entry table -  Learn rel d1,d5 and rel d1,d6
data_input = [ttest_stats.cent_si_learnA1v5; ttest_stats.cent_si_learnA1v6];
comp_descrip_in = {'Centroid difference relative to d1 - learning  - A - rel d1 vs. rel d5';...
                    'Centroid difference relative to d1 - learning  - A - rel d1 vs. rel d6'};
[t_ttest.cent_si_learnA15v16] = paired_ttest_table_entry(data_input,...
        9, 'd', 'by animal', comp_descrip_in);

%paired t-test comparisons (recall)
%rel d1 vs d5 recall
ttest_stats.cent_si_recallA1v5 = paired_ttest(cent_si_recall.A(:,2), cent_si_recall.A(:,6));

%rel d1 vs d6 recall
ttest_stats.cent_si_recallA1v6 = paired_ttest(cent_si_recall.A(:,2), cent_si_recall.A(:,7));

%paired t-test single-entry table -  Learn rel d1,d5 and rel d1,d6
data_input = [ttest_stats.cent_si_recallA1v5; ttest_stats.cent_si_recallA1v6];
comp_descrip_in = {'Centroid difference relative to d1 - recall  - A - rel d1 vs. rel d5';...
                    'Centroid difference relative to d1 - recall  - A - rel d1 vs. rel d6'};
[t_ttest.cent_si_recallA15v16] = paired_ttest_table_entry(data_input,...
        9, 'd', 'by animal', comp_descrip_in);
    
    
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%cent diff learning vs recall A - d5
[ttest_stats.cent_si_learn_recall_15A] = unpaired_ttest(cent_si_learn.A(:,6), cent_si_recall.A(:,6));

%cent diff learning vs recall A - d6
[ttest_stats.cent_si_learn_recall_16A] = unpaired_ttest(cent_si_learn.A(:,7), cent_si_recall.A(:,7));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.cent_si_learn_recall_15A; ttest_stats.cent_si_learn_recall_16A];

comp_descrip_in = {'Centroid difference relative to d1 - learning vs. recall  - A - rel d5';...
                'Centroid difference relative to d1 - learning vs. recall  - A - rel d6'};

[t_ttest.cent_si_learn_recallA_15_16] = unpaired_ttest_table_entry(data_input,...
        9, 'd', 'by animal', comp_descrip_in);


%% Centroid SI distance difference B
%Centroid SI neighbor B

cent_si_learn.B = source_data_short_learn_recall.angle_diff.si.st_learn.animal.mean.B;
%same data as in raw sub-struct
cent_si_recall.B = source_data_short_learn_recall.angle_diff.si.st_recall.animal.mean.B;


%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = cent_si_learn.B(:,[2,3,6,7]);
data_in_2 = cent_si_recall.B(:,[2,3,6,7]);

%output stats for 2way RM LME analysis
lme_stats.cent_si_learn_recall.B = two_way_rm_lme(data_in_1,data_in_2);

%create formatting table for 2-way analysis
[t_2way_rm_lme.cent_si_learn_recallB] = two_way_lme_table_entry(9,'d','by animal',...
                'Centroid difference relative to d1 - learning vs. recall - B - rel d1, rel d2, rel d5, rel d6',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.cent_si_learn_recall.B);       

%paired t-test comparisons (learn)
%rel d1 vs d5 learning
ttest_stats.cent_si_learnB1v5 = paired_ttest(cent_si_learn.B(:,2), cent_si_learn.B(:,6));

%rel d1 vs d6 learning
ttest_stats.cent_si_learnB1v6 = paired_ttest(cent_si_learn.B(:,2), cent_si_learn.B(:,7));

%paired t-test single-entry table -  Learn rel d1,d5 and rel d1,d6
data_input = [ttest_stats.cent_si_learnB1v5; ttest_stats.cent_si_learnB1v6];
comp_descrip_in = {'Centroid difference relative to d1 - learning  - B - rel d1 vs. rel d5';...
                    'Centroid difference relative to d1 - learning  - B - rel d1 vs. rel d6'};
[t_ttest.cent_si_learnB15v16] = paired_ttest_table_entry(data_input,...
        9, 'd', 'by animal', comp_descrip_in);

%paired t-test comparisons (recall)
%rel d1 vs d5 recall
ttest_stats.cent_si_recallB1v5 = paired_ttest(cent_si_recall.B(:,2), cent_si_recall.B(:,6));

%rel d1 vs d6 recall
ttest_stats.cent_si_recallB1v6 = paired_ttest(cent_si_recall.B(:,2), cent_si_recall.B(:,7));

%paired t-test single-entry table -  Learn rel d1,d5 and rel d1,d6
data_input = [ttest_stats.cent_si_recallB1v5; ttest_stats.cent_si_recallB1v6];
comp_descrip_in = {'Centroid difference relative to d1 - recall  - B - rel d1 vs. rel d5';...
                    'Centroid difference relative to d1 - recall  - B - rel d1 vs. rel d6'};
[t_ttest.cent_si_recallB15v16] = paired_ttest_table_entry(data_input,...
        9, 'd', 'by animal', comp_descrip_in);
    
    
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%cent diff learning vs recall A - d5
[ttest_stats.cent_si_learn_recall_15B] = unpaired_ttest(cent_si_learn.B(:,6), cent_si_recall.B(:,6));

%cent diff learning vs recall A - d6
[ttest_stats.cent_si_learn_recall_16B] = unpaired_ttest(cent_si_learn.B(:,7), cent_si_recall.B(:,7));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.cent_si_learn_recall_15B; ttest_stats.cent_si_learn_recall_16B];

comp_descrip_in = {'Centroid difference relative to d1 - learning vs. recall  - B - rel d5';...
                'Centroid difference relative to d1 - learning vs. recall  - B - rel d6'};

[t_ttest.cent_si_learn_recallB_15_16] = unpaired_ttest_table_entry(data_input,...
        9, 'd', 'by animal', comp_descrip_in);

%% Ex. Figure 9 Export statistics - short term learn SI - 
% A&B correlation between A and B laps across time

%import TS data for each day
AB_si_corr_learn = source_data_short_learn_recall.corr_analysis.st_learn.AB_corr.si.raw.pooled.AB_corr_ratio;

%remove nan values from each column of input data
for ii=1:numel(AB_si_corr_learn)
    temp_idx = find(isnan(AB_si_corr_learn{ii}) == 1);
    if ~isempty(temp_idx)
        AB_si_corr_learn{ii}(temp_idx) = [];
    end
end

%input data for stats analysis - choose which days to use to KW test
data_input = AB_si_corr_learn(2:7)';

%run KW test
[kruskal_stats_learn] = kruskalwallis_test(data_input);

%generate single row entry for KW test
comp_descrip_in = {'SI - pooled neuron - group difference - Learn - D2,D3,D4,D5,D6,D7'};
[t_krusall_learn] = kruskal_wallis_single_table_entry(9,'e','pooled',...
                comp_descrip_in,kruskal_stats_learn);

%for 1 sample Wilcoxon tests
for ii=2:7
    paired_wilcoxon_stats_learn{ii} = paired_wilcoxon_signrank(AB_si_corr_learn{ii},1);
end

%test comparison description for each table entry \
comp_descrip_in = {'AB Corr vs. 1, Day 2 - Learn';...
                    'AB Corr vs. 1, Day 3 - Learn';...
                    'AB Corr vs. 1, Day 4 - Learn';...
                    'AB Corr vs. 1, Day 5 - Learn';...
                    'AB Corr vs. 1, Day 6 - Learn';...
                    'AB Corr vs. 1, Day 7 - Learn'};

%table entries for paired wilcoxon test learning
paired_wilcox_AB_learn = paired_wilcoxon_table_entry_no_adj(9,'e','pooled',...
                comp_descrip_in,paired_wilcoxon_stats_learn(2:7));

%% Ex. Figure 9 Export statistics - short term recall SI 
% A&B correlation between A and B laps across time

%recall
AB_si_corr_recall = source_data_short_learn_recall.corr_analysis.st_recall.AB_corr.si.raw.pooled.AB_corr_ratio;

%remove nan values from each column of input data
for ii=1:numel(AB_si_corr_recall)
    temp_idx = find(isnan(AB_si_corr_recall{ii}) == 1);
    if ~isempty(temp_idx)
        AB_si_corr_recall{ii}(temp_idx) = [];
    end
end

%input data for stats analysis - choose which days to use to KW test
data_input = AB_si_corr_recall([2,3,6,7,8,9])';

%run KW test
[kruskal_stats_recall] = kruskalwallis_test(data_input);

%generate single row entry for KW test
comp_descrip_in = {'SI - pooled neuron - group difference - Recall - D2,D3,D6,D7,D8,D9'};
[t_krusall_recall] = kruskal_wallis_single_table_entry(9,'e','pooled',...
                comp_descrip_in,kruskal_stats_recall);

%for 1 sample Wilcoxon tests
for ii=[2,3,6,7,8,9]
    paired_wilcoxon_stats_recall{ii} = paired_wilcoxon_signrank(AB_si_corr_recall{ii},1);
end

%test comparison description for each table entry \
comp_descrip_in = {'AB Corr vs. 1, Day 2 - Recall';...
                    'AB Corr vs. 1, Day 3 - Recall';...
                    'AB Corr vs. 1, Day 6 - Recall';...
                    'AB Corr vs. 1, Day 7 - Recall';...
                    'AB Corr vs. 1, Day 8 - Recall';...
                    'AB Corr vs. 1, Day 9 - Recall'};

%table entries for paired wilcoxon test learning
%table entries for paired wilcoxon test learning
paired_wilcox_AB_recall = paired_wilcoxon_table_entry_no_adj(9,'e','pooled',...
                comp_descrip_in,paired_wilcoxon_stats_recall([2,3,6,7,8,9]));

%% Figure 10 Long term recall data SI and TS criteria

%10a -fraction tuned across time - SI
%unload data (animal x day)
long_recall_si.A = source_data_short_learn_recall.si_frac_export.long_recall.A;
long_recall_si.B = source_data_short_learn_recall.si_frac_export.long_recall.B;
long_recall_si.AB = source_data_short_learn_recall.si_frac_export.long_recall.AB;

long_recall_ts.A = source_data_short_learn_recall.ts_frac_export.long_recall.A;
long_recall_ts.B = source_data_short_learn_recall.ts_frac_export.long_recall.B;
long_recall_ts.AB = source_data_short_learn_recall.ts_frac_export.long_recall.AB;

%1-way RM LME analysis for learn
%input data arranged as animal x day

%A
%number of time points
nb_time_points = 6;
%number of animals
nb_animals = 3;
[lme_stats.long_si_A] = one_way_RM_lme(long_recall_si.A,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_si_A] = one_way_lme_table_entry(10,'a','by animal',...
                'A tuned fraction across time - SI - Recall (long)',...
                nb_animals,lme_stats.long_si_A);
          
%paired t-test D1 vs. D25
t_stats.long_si_A_1v25 = paired_ttest(long_recall_si.A(:,1), long_recall_si.A(:,5));

%paired t-test D1 vs. D30
t_stats.long_si_A_1v30 = paired_ttest(long_recall_si.A(:,1), long_recall_si.A(:,6));


%paired t-test D1 vs. D30
data_input = [t_stats.long_si_A_1v25; t_stats.long_si_A_1v30];

%what comparisons are being made
comp_descrip_in = {'A tuned fraction across time - SI - Recall D1 vs. D25 (long)';...
                'A tuned fraction across time - SI - Recall D1 vs. D30 (long)'};

%generate table entry
[t_ttest.long_si_A_25_30] = paired_ttest_table_entry(data_input,...
        10, 'a', 'by animal', comp_descrip_in);    
    
%%%%B    
[lme_stats.long_si_B] = one_way_RM_lme(long_recall_si.B,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_si_B] = one_way_lme_table_entry(10,'a','by animal',...
                'B tuned fraction across time - SI - Recall (long)',...
                nb_animals,lme_stats.long_si_B);
          
%paired t-test D1 vs. D25
t_stats.long_si_B_1v25 = paired_ttest(long_recall_si.B(:,1), long_recall_si.B(:,5));

%paired t-test D1 vs. D30
t_stats.long_si_B_1v30 = paired_ttest(long_recall_si.B(:,1), long_recall_si.B(:,6));


%paired t-test D1 vs. D30
data_input = [t_stats.long_si_B_1v25; t_stats.long_si_B_1v30];

%what comparisons are being made
comp_descrip_in = {'B tuned fraction across time - SI - Recall D1 vs. D25 (long)';...
                'B tuned fraction across time - SI - Recall D1 vs. D30 (long)'};

[t_ttest.long_si_B_25_30] = paired_ttest_table_entry(data_input,...
        10, 'a', 'by animal', comp_descrip_in);        
    
%A&B    
[lme_stats.long_si_AB] = one_way_RM_lme(long_recall_si.AB,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_si_AB] = one_way_lme_table_entry(10,'a','by animal',...
                'AB tuned fraction across time - SI - Recall (long)',...
                nb_animals,lme_stats.long_si_AB);
          
%paired t-test D1 vs. D25
t_stats.long_si_AB_1v25 = paired_ttest(long_recall_si.AB(:,1), long_recall_si.AB(:,5));

%paired t-test D1 vs. D30
t_stats.long_si_AB_1v30 = paired_ttest(long_recall_si.AB(:,1), long_recall_si.AB(:,6));


%paired t-test D1 vs. D30
data_input = [t_stats.long_si_AB_1v25; t_stats.long_si_AB_1v30];

%what comparisons are being made
comp_descrip_in = {'AB tuned fraction across time - SI - Recall D1 vs. D25 (long)';...
                'AB tuned fraction across time - SI - Recall D1 vs. D30 (long)'};

[t_ttest.long_si_AB_25_30] = paired_ttest_table_entry(data_input,...
        10, 'a', 'by animal', comp_descrip_in);    
            
%% 10a - fraction tuned across time - TS
[lme_stats.long_ts_A] = one_way_RM_lme(long_recall_ts.A,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_ts_A] = one_way_lme_table_entry(10,'a','by animal',...
                'A tuned fraction across time - TS - Recall (long)',...
                nb_animals,lme_stats.long_ts_A);
          
%paired t-test D1 vs. D25
t_stats.long_ts_A_1v25 = paired_ttest(long_recall_ts.A(:,1), long_recall_ts.A(:,5));

%paired t-test D1 vs. D30
t_stats.long_ts_A_1v30 = paired_ttest(long_recall_ts.A(:,1), long_recall_ts.A(:,6));


%paired t-test D1 vs. D30
data_input = [t_stats.long_ts_A_1v25; t_stats.long_ts_A_1v30];

%what comparisons are being made
comp_descrip_in = {'A tuned fraction across time - TS - Recall D1 vs. D25 (long)';...
                'A tuned fraction across time - TS - Recall D1 vs. D30 (long)'};

%generate table entry
[t_ttest.long_ts_A_25_30] = paired_ttest_table_entry(data_input,...
        10, 'a', 'by animal', comp_descrip_in);    
    
%%%%B    
[lme_stats.long_ts_B] = one_way_RM_lme(long_recall_ts.B,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_ts_B] = one_way_lme_table_entry(10,'a','by animal',...
                'B tuned fraction across time - TS - Recall (long)',...
                nb_animals,lme_stats.long_ts_B);
          
%paired t-test D1 vs. D25
t_stats.long_ts_B_1v25 = paired_ttest(long_recall_ts.B(:,1), long_recall_ts.B(:,5));

%paired t-test D1 vs. D30
t_stats.long_ts_B_1v30 = paired_ttest(long_recall_ts.B(:,1), long_recall_ts.B(:,6));


%paired t-test D1 vs. D30
data_input = [t_stats.long_ts_B_1v25; t_stats.long_ts_B_1v30];

%what comparisons are being made
comp_descrip_in = {'B tuned fraction across time - TS - Recall D1 vs. D25 (long)';...
                'B tuned fraction across time - TS - Recall D1 vs. D30 (long)'};

[t_ttest.long_ts_B_25_30] = paired_ttest_table_entry(data_input,...
        10, 'a', 'by animal', comp_descrip_in);        
    
%A&B    
[lme_stats.long_ts_AB] = one_way_RM_lme(long_recall_ts.AB,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_ts_AB] = one_way_lme_table_entry(10,'a','by animal',...
                'A&B tuned fraction across time - TS - Recall (long)',...
                nb_animals,lme_stats.long_ts_AB);
          
%paired t-test D1 vs. D25
t_stats.long_ts_AB_1v25 = paired_ttest(long_recall_ts.AB(:,1), long_recall_ts.AB(:,5));

%paired t-test D1 vs. D30
t_stats.long_ts_AB_1v30 = paired_ttest(long_recall_ts.AB(:,1), long_recall_ts.AB(:,6));


%paired t-test D1 vs. D30
data_input = [t_stats.long_ts_AB_1v25; t_stats.long_ts_AB_1v30];

%what comparisons are being made
comp_descrip_in = {'A&B tuned fraction across time - TS - Recall D1 vs. D25 (long)';...
                'A&B tuned fraction across time - TS - Recall D1 vs. D30 (long)'};

[t_ttest.long_ts_AB_25_30] = paired_ttest_table_entry(data_input,...
        10, 'a', 'by animal', comp_descrip_in); 
    
    
%% 10c - PV A and B rel d1

%excluded 1 vs 1 correlation (d6, d16,d20,d25,d30)
long_recall_PV.A = source_data_short_learn_recall.corr_analysis.lt_recall.raw.day_PV_diag_mean.A(:,[6,16,20,25,30]);
long_recall_PV.B = source_data_short_learn_recall.corr_analysis.lt_recall.raw.day_PV_diag_mean.B(:,[6,16,20,25,30]);

%A
%number of time points (excluded self-day correction which = 1)
nb_time_points = 5;
%number of animals
nb_animals = 3;
%run 1-way LME
[lme_stats.long_PV_A] = one_way_RM_lme(long_recall_PV.A,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_PV_A] = one_way_lme_table_entry(10,'c','by animal',...
                'Population vector correlation relative to D1  - A laps - d6, d16,d20,d26,d30',...
                nb_animals,lme_stats.long_PV_A);
          
%paired t-test D6 vs. D25
t_stats.long_PV_A_6v25 = paired_ttest(long_recall_PV.A(:,1), long_recall_PV.A(:,4));

%paired t-test D6 vs. D30
t_stats.long_PV_A_6v30 = paired_ttest(long_recall_PV.A(:,1), long_recall_PV.A(:,5));


%paired t-test D1 vs. D30
data_input = [t_stats.long_PV_A_6v25; t_stats.long_PV_A_6v30];

%what comparisons are being made
comp_descrip_in = {'Recall PV correlation - time -D6 vs. D25 - A laps';...
                'Recall PV correlation - time -D6 vs. D30 - A laps)'};

%generate table entry
[t_ttest.long_PV_A_25_30] = paired_ttest_table_entry(data_input,...
        10, 'c', 'by animal', comp_descrip_in); 

%PV B
    %run 1-way LME
[lme_stats.long_PV_B] = one_way_RM_lme(long_recall_PV.B,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_PV_B] = one_way_lme_table_entry(10,'c','by animal',...
                'Population vector correlation relative to D1  - B laps - d6, d16,d20,d26,d30',...
                nb_animals,lme_stats.long_PV_B);
          
%paired t-test D6 vs. D25
t_stats.long_PV_B_6v25 = paired_ttest(long_recall_PV.B(:,1), long_recall_PV.B(:,4));

%paired t-test D6 vs. D30
t_stats.long_PV_B_6v30 = paired_ttest(long_recall_PV.B(:,1), long_recall_PV.B(:,5));


%paired t-test D1 vs. D30
data_input = [t_stats.long_PV_B_6v25; t_stats.long_PV_B_6v30];

%what comparisons are being made
comp_descrip_in = {'Recall PV correlation - time -D6 vs. D25 - B laps';...
                'Recall PV correlation - time -D6 vs. D30 - B laps)'};

%generate table entry
[t_ttest.long_PV_B_25_30] = paired_ttest_table_entry(data_input,...
        10, 'c', 'by animal', comp_descrip_in); 
    
    
%% 10c - TC SI - A and B rel d1

%excluded 1 vs 1 correlation (d6, d16,d20,d25,d30)
long_recall_TC_si.A = source_data_short_learn_recall.corr_analysis.lt_recall.si.raw.animal.day_TC_diag_mean.exp.A(:,[6,16,20,25,30]);
long_recall_TC_si.B = source_data_short_learn_recall.corr_analysis.lt_recall.si.raw.animal.day_TC_diag_mean.exp.B(:,[6,16,20,25,30]);

%A
%number of time points (excluded self-day correction which = 1)
nb_time_points = 5;
%number of animals
nb_animals = 3;
%run 1-way LME
[lme_stats.long_TC_si_A] = one_way_RM_lme(long_recall_TC_si.A,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_TC_si_A] = one_way_lme_table_entry(10,'c','by animal',...
                'Tuning curve SI correlation relative to D1  - A laps',...
                nb_animals,lme_stats.long_TC_si_A);
          
%paired t-test D6 vs. D25
t_stats.long_TC_si_A_6v25 = paired_ttest(long_recall_TC_si.A(:,1), long_recall_TC_si.A(:,4));

%paired t-test D6 vs. D30
t_stats.long_TC_si_A_6v30 = paired_ttest(long_recall_TC_si.A(:,1), long_recall_TC_si.A(:,5));


%paired t-test D1 vs. D30
data_input = [t_stats.long_TC_si_A_6v25; t_stats.long_TC_si_A_6v30];

%what comparisons are being made
comp_descrip_in = {'Recall Tuning curve SI correlation - time -D6 vs. D25 - A laps';...
                'Recall Tuning curve SI correlation - time -D6 vs. D30 - A laps'};

%generate table entry
[t_ttest.long_TC_si_A_25_30] = paired_ttest_table_entry(data_input,...
        10, 'c', 'by animal', comp_descrip_in); 

%TC_si B
    %run 1-way LME
[lme_stats.long_TC_si_B] = one_way_RM_lme(long_recall_TC_si.B,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_TC_si_B] = one_way_lme_table_entry(10,'c','by animal',...
                'Tuning curve SI correlation relative to D1  - B laps',...
                nb_animals,lme_stats.long_TC_si_B);
          
%paired t-test D6 vs. D25
t_stats.long_TC_si_B_6v25 = paired_ttest(long_recall_TC_si.B(:,1), long_recall_TC_si.B(:,4));

%paired t-test D6 vs. D30
t_stats.long_TC_si_B_6v30 = paired_ttest(long_recall_TC_si.B(:,1), long_recall_TC_si.B(:,5));


%paired t-test D1 vs. D30
data_input = [t_stats.long_TC_si_B_6v25; t_stats.long_TC_si_B_6v30];

%what comparisons are being made
comp_descrip_in = {'Recall Tuning curve SI correlation - time -D6 vs. D25 - B laps';...
                'Recall Tuning curve SI correlation - time -D6 vs. D30 - B laps'};

%generate table entry
[t_ttest.long_TC_si_B_25_30] = paired_ttest_table_entry(data_input,...
        10, 'c', 'by animal', comp_descrip_in); 
    
%% 10c - TC TS - A nad B rel d1
%excluded 1 vs 1 correlation (d6, d16,d20,d25,d30)
long_recall_TC_ts.A = source_data_short_learn_recall.corr_analysis.lt_recall.ts.raw.animal.day_TC_diag_mean.exp.A(:,[6,16,20,25,30]);
long_recall_TC_ts.B = source_data_short_learn_recall.corr_analysis.lt_recall.ts.raw.animal.day_TC_diag_mean.exp.B(:,[6,16,20,25,30]);

%A
%number of time points (excluded self-day correction which = 1)
nb_time_points = 5;
%number of animals
nb_animals = 3;
%run 1-way LME
[lme_stats.long_TC_ts_A] = one_way_RM_lme(long_recall_TC_ts.A,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_TC_ts_A] = one_way_lme_table_entry(10,'c','by animal',...
                'Tuning curve TS correlation relative to D1  - A laps',...
                nb_animals,lme_stats.long_TC_ts_A);
          
%paired t-test D6 vs. D25
t_stats.long_TC_ts_A_6v25 = paired_ttest(long_recall_TC_ts.A(:,1), long_recall_TC_ts.A(:,4));

%paired t-test D6 vs. D30
t_stats.long_TC_ts_A_6v30 = paired_ttest(long_recall_TC_ts.A(:,1), long_recall_TC_ts.A(:,5));


%paired t-test D1 vs. D30
data_input = [t_stats.long_TC_ts_A_6v25; t_stats.long_TC_ts_A_6v30];

%what comparisons are being made
comp_descrip_in = {'Recall Tuning curve TS correlation - time -D6 vs. D25 - A laps';...
                'Recall Tuning curve TS correlation - time -D6 vs. D30 - A laps'};

%generate table entry
[t_ttest.long_TC_ts_A_25_30] = paired_ttest_table_entry(data_input,...
        10, 'c', 'by animal', comp_descrip_in); 

%TC_ts B
    %run 1-way LME
[lme_stats.long_TC_ts_B] = one_way_RM_lme(long_recall_TC_ts.B,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_TC_ts_B] = one_way_lme_table_entry(10,'c','by animal',...
                'Tuning curve TS correlation relative to D1  - B laps',...
                nb_animals,lme_stats.long_TC_ts_B);
          
%paired t-test D6 vs. D25
t_stats.long_TC_ts_B_6v25 = paired_ttest(long_recall_TC_ts.B(:,1), long_recall_TC_ts.B(:,4));

%paired t-test D6 vs. D30
t_stats.long_TC_ts_B_6v30 = paired_ttest(long_recall_TC_ts.B(:,1), long_recall_TC_ts.B(:,5));


%paired t-test D1 vs. D30
data_input = [t_stats.long_TC_ts_B_6v25; t_stats.long_TC_ts_B_6v30];

%what comparisons are being made
comp_descrip_in = {'Recall Tuning curve TS correlation - time -D6 vs. D25 - B laps';...
                'Recall Tuning curve TS correlation - time -D6 vs. D30 - B laps'};

%generate table entry
[t_ttest.long_TC_ts_B_25_30] = paired_ttest_table_entry(data_input,...
        10, 'c', 'by animal', comp_descrip_in); 

%% 10d - PV neighbor

%excluded 1 vs 1 correlation (d6, d16,d20,d25,d30)
long_recall_PVn.A = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.raw.day_PV_diag_mean.A;
long_recall_PVn.B = source_data_short_learn_recall.corr_analysis.lt_recall.neighbor.raw.day_PV_diag_mean.B;

%A
%number of time points (excluded self-day correction which = 1)
nb_time_points = 5;
%number of animals
nb_animals = 3;
%run 1-way LME
[lme_stats.long_PVn_A] = one_way_RM_lme(long_recall_PVn.A,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_PVn_A] = one_way_lme_table_entry(10,'d','by animal',...
                'Population vector correlation neighboring days -  recall - A, 1 vs. 6, 6 vs 16, 16 vs 20, 20 vs 25, 25 vs 30',...
                nb_animals,lme_stats.long_PVn_A);
          
%paired t-test D1-6 vs. D20-25
t_stats.long_PVn_A_6v25 = paired_ttest(long_recall_PVn.A(:,1), long_recall_PVn.A(:,4));

%paired t-test D1-6 vs. D25-30
t_stats.long_PVn_A_6v30 = paired_ttest(long_recall_PVn.A(:,1), long_recall_PVn.A(:,5));


%paired t-test D1 vs. D25/D30
data_input = [t_stats.long_PVn_A_6v25; t_stats.long_PVn_A_6v30];

%what comparisons are being made
comp_descrip_in = {'PV correlation recall A - time - 1 vs. 6 vs. 20 vs. 25';...
                'PV correlation recall A - time - 1 vs. 6 vs. 25 vs. 30'};

%generate table entry
[t_ttest.long_PVn_A_25_30] = paired_ttest_table_entry(data_input,...
        10, 'd', 'by animal', comp_descrip_in); 

%PV B
    %run 1-way LME
[lme_stats.long_PVn_B] = one_way_RM_lme(long_recall_PVn.B,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_PVn_B] = one_way_lme_table_entry(10,'c','by animal',...
                'Population vector correlation neighboring days -  recall - B, 1 vs. 6, 6 vs 16, 16 vs 20, 20 vs 25, 25 vs 30',...
                nb_animals,lme_stats.long_PVn_B);
          
%paired t-test D6 vs. D25
t_stats.long_PVn_B_6v25 = paired_ttest(long_recall_PVn.B(:,1), long_recall_PVn.B(:,4));

%paired t-test D6 vs. D30
t_stats.long_PVn_B_6v30 = paired_ttest(long_recall_PVn.B(:,1), long_recall_PVn.B(:,5));


%paired t-test D1 vs. D30
data_input = [t_stats.long_PVn_B_6v25; t_stats.long_PVn_B_6v30];

%what comparisons are being made
comp_descrip_in = {'PV correlation recall B - time - 1 vs. 6 vs. 20 vs. 25';...
                'PV correlation recall B - time - 1 vs. 6 vs. 25 vs. 30'};

%generate table entry
[t_ttest.long_PVn_B_25_30] = paired_ttest_table_entry(data_input,...
        10, 'd', 'by animal', comp_descrip_in); 
    
%% 10d - TC SI neighbor
%excluded 1 vs 1 correlation (d6, d16,d20,d25,d30)
long_recall_TCn_si.A = source_data_short_learn_recall.TC.neighbor.si.lt_recall.A;
long_recall_TCn_si.B = source_data_short_learn_recall.TC.neighbor.si.lt_recall.B;

%A
%number of time points (excluded self-day correction which = 1)
nb_time_points = 5;
%number of animals
nb_animals = 3;
%run 1-way LME
[lme_stats.long_TCn_si_A] = one_way_RM_lme(long_recall_TCn_si.A,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_TCn_si_A] = one_way_lme_table_entry(10,'d','by animal',...
                'Tuning curve (SI) correlation neighboring days -  recall - A, 1 vs. 6, 6 vs 16, 16 vs 20, 20 vs 25, 25 vs 30',...
                nb_animals,lme_stats.long_TCn_si_A);
          
%paired t-test D6 vs. D25
t_stats.long_TCn_si_A_6v25 = paired_ttest(long_recall_TCn_si.A(:,1), long_recall_TCn_si.A(:,4));

%paired t-test D6 vs. D30
t_stats.long_TCn_si_A_6v30 = paired_ttest(long_recall_TCn_si.A(:,1), long_recall_TCn_si.A(:,5));


%paired t-test D1 vs. D30
data_input = [t_stats.long_TCn_si_A_6v25; t_stats.long_TCn_si_A_6v30];

%what comparisons are being made
comp_descrip_in = {'TC correlation recall A - time - 1 vs. 6 vs. 20 vs. 25';...
                'TC correlation recall A - time - 1 vs. 6 vs. 25 vs. 30'};

%generate table entry
[t_ttest.long_TCn_si_A_25_30] = paired_ttest_table_entry(data_input,...
        10, 'd', 'by animal', comp_descrip_in); 

%TCn_si B
    %run 1-way LME
[lme_stats.long_TCn_si_B] = one_way_RM_lme(long_recall_TCn_si.B,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_TCn_si_B] = one_way_lme_table_entry(10,'d','by animal',...
                'Tuning curve (SI) correlation neighboring days -  recall - B, 1 vs. 6, 6 vs 16, 16 vs 20, 20 vs 25, 25 vs 30',...
                nb_animals,lme_stats.long_TCn_si_B);
          
%paired t-test D6 vs. D25
t_stats.long_TCn_si_B_6v25 = paired_ttest(long_recall_TCn_si.B(:,1), long_recall_TCn_si.B(:,4));

%paired t-test D6 vs. D30
t_stats.long_TCn_si_B_6v30 = paired_ttest(long_recall_TCn_si.B(:,1), long_recall_TCn_si.B(:,5));


%paired t-test D1 vs. D30
data_input = [t_stats.long_TCn_si_B_6v25; t_stats.long_TCn_si_B_6v30];

%what comparisons are being made
comp_descrip_in = {'TC correlation recall B - time - 1 vs. 6 vs. 20 vs. 25';...
                'TC correlation recall B - time - 1 vs. 6 vs. 25 vs. 30'};

%generate table entry
[t_ttest.long_TCn_si_B_25_30] = paired_ttest_table_entry(data_input,...
        10, 'd', 'by animal', comp_descrip_in); 
%% 10d - TC TS neighbor
%excluded 1 vs 1 correlation (d6, d16,d20,d25,d30)
long_recall_TCn_ts.A = source_data_short_learn_recall.TC.neighbor.ts.lt_recall.A;
long_recall_TCn_ts.B = source_data_short_learn_recall.TC.neighbor.ts.lt_recall.B;

%A
%number of time points (excluded self-day correction which = 1)
nb_time_points = 5;
%number of animals
nb_animals = 3;
%run 1-way LME
[lme_stats.long_TCn_ts_A] = one_way_RM_lme(long_recall_TCn_ts.A,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_TCn_ts_A] = one_way_lme_table_entry(10,'d','by animal',...
                'Tuning curve (TS) correlation neighboring days -  recall - A, 1 vs. 6, 6 vs 16, 16 vs 20, 20 vs 25, 25 vs 30',...
                nb_animals,lme_stats.long_TCn_ts_A);
          
%paired t-test D6 vs. D25
t_stats.long_TCn_ts_A_6v25 = paired_ttest(long_recall_TCn_ts.A(:,1), long_recall_TCn_ts.A(:,4));

%paired t-test D6 vs. D30
t_stats.long_TCn_ts_A_6v30 = paired_ttest(long_recall_TCn_ts.A(:,1), long_recall_TCn_ts.A(:,5));


%paired t-test D1 vs. D30
data_input = [t_stats.long_TCn_ts_A_6v25; t_stats.long_TCn_ts_A_6v30];

%what comparisons are being made
comp_descrip_in = {'TC correlation recall A - time - 1 vs. 6 vs. 20 vs. 25';...
                'TC correlation recall A - time - 1 vs. 6 vs. 25 vs. 30'};

%generate table entry
[t_ttest.long_TCn_ts_A_25_30] = paired_ttest_table_entry(data_input,...
        10, 'd', 'by animal', comp_descrip_in); 

%TCn_si B
    %run 1-way LME
[lme_stats.long_TCn_ts_B] = one_way_RM_lme(long_recall_TCn_ts.B,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.long_TCn_ts_B] = one_way_lme_table_entry(10,'d','by animal',...
                'Tuning curve (TS) correlation neighboring days -  recall - B, 1 vs. 6, 6 vs 16, 16 vs 20, 20 vs 25, 25 vs 30',...
                nb_animals,lme_stats.long_TCn_ts_B);
          
%paired t-test D6 vs. D25
t_stats.long_TCn_ts_B_6v25 = paired_ttest(long_recall_TCn_ts.B(:,1), long_recall_TCn_ts.B(:,4));

%paired t-test D6 vs. D30
t_stats.long_TCn_ts_B_6v30 = paired_ttest(long_recall_TCn_ts.B(:,1), long_recall_TCn_ts.B(:,5));


%paired t-test D1 vs. D30
data_input = [t_stats.long_TCn_ts_B_6v25; t_stats.long_TCn_ts_B_6v30];

%what comparisons are being made
comp_descrip_in = {'TC correlation recall B - time - 1 vs. 6 vs. 20 vs. 25';...
                'TC correlation recall B - time - 1 vs. 6 vs. 25 vs. 30'};

%generate table entry
[t_ttest.long_TCn_ts_B_25_30] = paired_ttest_table_entry(data_input,...
        10, 'd', 'by animal', comp_descrip_in);
    
%% 10e A&B tuned A vs B lap corr recall SI
%recall
AB_si_corr_recall = source_data_short_learn_recall.corr_analysis.lt_recall.AB_corr.si.raw.pooled.AB_corr_ratio;

%remove nan values from each column of input data
for ii=1:numel(AB_si_corr_recall)
    temp_idx = find(isnan(AB_si_corr_recall{ii}) == 1);
    if ~isempty(temp_idx)
        AB_si_corr_recall{ii}(temp_idx) = [];
    end
end

%input data for stats analysis - choose which days to use to KW test
data_input = AB_si_corr_recall([2,3,4,5,6])';

%run KW test
[kruskal_stats_recall] = kruskalwallis_test(data_input);

%generate single row entry for KW test
comp_descrip_in = {'SI - pooled neuron - group difference - Recall - D6,D16,D20,D25,D30'};
[t_krusall_recall_lt_si] = kruskal_wallis_single_table_entry(10,'e','pooled',...
                comp_descrip_in,kruskal_stats_recall);

clear paired_wilcoxon_stats_recall
%for 1 sample Wilcoxon tests
for ii=[2,3,4,5,6]
    paired_wilcoxon_stats_recall{ii} = paired_wilcoxon_signrank(AB_si_corr_recall{ii},1);
end

%test comparison description for each table entry \
comp_descrip_in = {'AB Corr vs. 1, Day 6 - Recall';...
                    'AB Corr vs. 1, Day 16 - Recall';...
                    'AB Corr vs. 1, Day 20 - Recall';...
                    'AB Corr vs. 1, Day 25 - Recall';...
                    'AB Corr vs. 1, Day 30 - Recall'};

%table entries for paired wilcoxon test learning
%table entries for paired wilcoxon test learning
t_paired_wilcox_AB_recall_lt_si = paired_wilcoxon_table_entry_no_adj(10,'e','pooled',...
                comp_descrip_in,paired_wilcoxon_stats_recall([2,3,4,5,6]));
            
%% 10e A&B tuned A vs B lap corr recall TS
%recall
AB_ts_corr_recall = source_data_short_learn_recall.corr_analysis.lt_recall.AB_corr.ts.raw.pooled.AB_corr_ratio;

%remove nan values from each column of input data
for ii=1:numel(AB_ts_corr_recall)
    temp_idx = find(isnan(AB_ts_corr_recall{ii}) == 1);
    if ~isempty(temp_idx)
        AB_ts_corr_recall{ii}(temp_idx) = [];
    end
end

%input data for stats analysis - choose which days to use to KW test
data_input = AB_ts_corr_recall([2,3,4,5,6])';

%run KW test
[kruskal_stats_recall] = kruskalwallis_test(data_input);

%generate single row entry for KW test
comp_descrip_in = {'TS - pooled neuron - group difference - Recall - D6,D16,D20,D25,D30'};
[t_krusall_recall_lt_ts] = kruskal_wallis_single_table_entry(10,'e','pooled',...
                comp_descrip_in,kruskal_stats_recall);

clear paired_wilcoxon_stats_recall
%for 1 sample Wilcoxon tests
for ii=[2,3,4,5,6]
    paired_wilcoxon_stats_recall{ii} = paired_wilcoxon_signrank(AB_ts_corr_recall{ii},1);
end

%test comparison description for each table entry \
comp_descrip_in = {'AB Corr vs. 1, Day 6 - Recall';...
                    'AB Corr vs. 1, Day 16 - Recall';...
                    'AB Corr vs. 1, Day 20 - Recall';...
                    'AB Corr vs. 1, Day 25 - Recall';...
                    'AB Corr vs. 1, Day 30 - Recall'};

%table entries for paired wilcoxon test learning
t_paired_wilcox_AB_recall_lt_ts = paired_wilcoxon_table_entry_no_adj(10,'e','pooled',...
                comp_descrip_in,paired_wilcoxon_stats_recall([2,3,4,5,6]));
            
%% Create Ex. Figure 4 stats export spreadsheet
%spreadsheet name
spreadsheet_name = 'statistics_summary.xlsx';

%sheet name
sheet_name = 'Extended Figure 4';

%empty row
t1 = repmat({' '},1,12);

%exported Excel spreadsheet
%write to Excel spreadsheet
%Ex Fig 4b
writetable(t_paired_wilcox_Asel,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','overwritesheet')
writetable(t_paired_wilcox_Bsel,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 4c - place field count between selective neurons
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_out.wilcox_singlePF,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_out.wilcox_doublePF,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_out.wilcox_triplePF,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 4c - place field width between task selective neurons
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ks2_pf_width,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%sheet name
sheet_name = 'Extended Figure 7';

%Ex Fig 7a and c
%writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_kstest_com_vs_act,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','overwritesheet')
writetable(t_paired_wilcox_activity_rate_com_activity,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%sheet name
sheet_name = 'Extended Figure 9';
%Ex Fig 9a - learn
writetable(t_1_rm_lme.learnA,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','overwritesheet')
writetable(t_ttest.learnA_6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.learnB,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.learnB_6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.learnAB,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.learnAB_6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%9a - recall 
writetable(t_1_rm_lme.recallA,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.recallA_6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_1_rm_lme.recallB,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.recallB_6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_1_rm_lme.recallAB,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.recallAB_6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%9b TC SI correlation across days - A trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.TC_si_learn_recallA,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_si_learn_recallA_6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_si_learnA2v6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_si_recallA2v6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%9b TC SI correlation across days -B trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.TC_si_learn_recallB,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_si_learn_recallB_6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_si_learnB2v6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_si_recallB2v6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%9c TC SI neighbor correlation across days - A trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.TC_si_n_learn_recallA,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_si_n_learnA12v67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_si_n_recallA12v67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_si_n_learn_recallA_12_67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%9c TC SI neighbor correlation across days - B trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.TC_si_n_learn_recallB,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_si_n_learnB12v67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_si_n_recallB12v67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_si_n_learn_recallB_12_67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%9d Centroid difference SI neighbor correlation across days - A trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.cent_si_learn_recallA,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.cent_si_learnA15v16,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.cent_si_recallA15v16,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.cent_si_learn_recallA_15_16,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%9d Centroid difference SI neighbor correlation across days - B trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.cent_si_learn_recallB,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.cent_si_learnB15v16,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.cent_si_recallB15v16,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.cent_si_learn_recallB_15_16,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%9e - A&B tuned A vs B lap tuning curve correlations 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_krusall_learn,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(paired_wilcox_AB_learn,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
%recall 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_krusall_recall,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(paired_wilcox_AB_recall,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
 

%sheet name
sheet_name = 'Extended Figure 10';
%Ex Fig 10a - SI
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_si_A,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','overwritesheet')
writetable(t_ttest.long_si_A_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_si_B,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_si_B_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_si_AB,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_si_AB_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10a - TS
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_ts_A,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_ts_A_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_ts_B,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_ts_B_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_ts_AB,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_ts_AB_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10c PV - A
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_PV_A,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_PV_A_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10c PV - B
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_PV_B,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_PV_B_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10c TC SI - A
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_TC_si_A,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_TC_si_A_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10c TC SI - B
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_TC_si_B,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_TC_si_B_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10c TC TS - A
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_TC_ts_A,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_TC_ts_A_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10c TC TS - B
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_TC_ts_B,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_TC_ts_B_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10d PV neighbor - A
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_PVn_A,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_PVn_A_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10d PV neighbor - B
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_PVn_B,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_PVn_B_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10d TC neighbor SI - A
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_TCn_si_A,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_TCn_si_A_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10d TC neighbor SI - B
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_TCn_si_B,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_TCn_si_B_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10d TC neighbor TS - A
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_TCn_ts_A,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_TCn_ts_A_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10d TC neighbor TS - B
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_1_rm_lme.long_TCn_ts_B,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_ttest.long_TCn_ts_B_25_30,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10e A&B corr -SI
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_krusall_recall_lt_si,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_paired_wilcox_AB_recall_lt_si,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 10e A&B corr - TS
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_krusall_recall_lt_ts,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_paired_wilcox_AB_recall_lt_ts,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
