%% Master script for generating stats for paper
%source data, import and format Prism analysis, and export to Excel, and
%Word legend data

%import data for each figure

%figure 2 source data
load('G:\Google_drive\task_selective_place_paper\matlab_data\source_data_fig4_5_and_sup.mat')

%laptop directory
%load('C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data\source_data_fig4_5_and_sup.mat')


%% Figure 4e Learning Analysis 
%fraction of task selective place cells over days in ST learn and recall

%unload data (animal x day)
learn.A = source_data_short_learn_recall.ts_frac_export.learning.A;
learn.B = source_data_short_learn_recall.ts_frac_export.learning.B;
learn.AB = source_data_short_learn_recall.ts_frac_export.learning.AB;

%1-way RM LME analysis for learn
%input data arranged as animal x day

%A
%number of time points
nb_time_points = 7;
%number of animals
nb_animals = 6;
[lme_stats.learnA] = one_way_RM_lme(learn.A,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.learnA] = one_way_lme_table_entry(4,'e','by animal',...
                'A tuned fraction across time - Learning D1-D7 - TS',...
                size(learn.A,1),lme_stats.learnA);
%B
[lme_stats.learnB] = one_way_RM_lme(learn.B,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.learnB] = one_way_lme_table_entry(4,'e','by animal',...
                'B tuned fraction across time - Learning D1-D7 - TS',...
                size(learn.B,1),lme_stats.learnB);

%AB
[lme_stats.learnAB] = one_way_RM_lme(learn.AB,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.learnAB] = one_way_lme_table_entry(4,'e','by animal',...
                'A&B tuned fraction across time - Learning D1-D7 - TS',...
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
comp_descrip_in = {'A tuned fraction across time - TS - Learning D1 vs. D6';...
                'A tuned fraction across time - TS - Learning D1 vs. D7'};

[t_ttest.learnA_6_7] = paired_ttest_table_entry(data_input,...
        4, 'e', 'by animal', comp_descrip_in);

%learn B D6/D7
data_input = [t_stats.learnB6; t_stats.learnB7];

%what comparisons are being made
comp_descrip_in = {'B tuned fraction across time - TS - Learning D1 vs. D6';...
                'B tuned fraction across time - TS - Learning D1 vs. D7'};

[t_ttest.learnB_6_7] = paired_ttest_table_entry(data_input,...
        4, 'e', 'by animal', comp_descrip_in);

%learn AB D6/D7
data_input = [t_stats.learnAB6; t_stats.learnAB7];

%what comparisons are being made
comp_descrip_in = {'A&B tuned fraction across time - TS - Learning D1 vs. D6';...
                'A&B tuned fraction across time - TS - Learning D1 vs. D7'};

[t_ttest.learnAB_6_7] = paired_ttest_table_entry(data_input,...
        4, 'e', 'by animal', comp_descrip_in);

%% Figure 4e Recall analysis

%short recall data (animal x day)
recall.A = source_data_short_learn_recall.ts_frac_export.recall.A;
recall.B = source_data_short_learn_recall.ts_frac_export.recall.B;
recall.AB = source_data_short_learn_recall.ts_frac_export.recall.AB;

%A
%number of time points
nb_time_points = 7;
%number of animals
nb_animals = 5;
[lme_stats.recallA] = one_way_RM_lme(recall.A,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.recallA] = one_way_lme_table_entry(4,'e','by animal',...
                'A tuned fraction across time - Recall D1-D9 - TS',...
                size(recall.A,1),lme_stats.recallA);
%B
[lme_stats.recallB] = one_way_RM_lme(recall.B,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.recallB] = one_way_lme_table_entry(4,'e','by animal',...
                'B tuned fraction across time - Recall D1-D9 - TS',...
                size(recall.B,1),lme_stats.recallB);

%AB
[lme_stats.recallAB] = one_way_RM_lme(recall.AB,nb_animals, nb_time_points);

%generate 1-row table entry with 1 RM LME analysis
[t_1_rm_lme.recallAB] = one_way_lme_table_entry(4,'e','by animal',...
                'A&B tuned fraction across time - Recall D1-D9 - TS',...
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
comp_descrip_in = {'A tuned fraction across time - TS - Recall D1 vs. D6';...
                'A tuned fraction across time - TS - Recall D1 vs. D7'};

[t_ttest.recallA_6_7] = paired_ttest_table_entry(data_input,...
        4, 'e', 'by animal', comp_descrip_in);

%learn B D6/D7
data_input = [t_stats.recallB6; t_stats.recallB7];

%what comparisons are being made
comp_descrip_in = {'B tuned fraction across time - TS - Recall D1 vs. D6';...
                'B tuned fraction across time - TS - Recall D1 vs. D7'};

[t_ttest.recallB_6_7] = paired_ttest_table_entry(data_input,...
        4, 'e', 'by animal', comp_descrip_in);

%recall AB D6/D7
data_input = [t_stats.recallAB6; t_stats.recallAB7];

%what comparisons are being made
comp_descrip_in = {'A&B tuned fraction across time - TS - Recall D1 vs. D6';...
                'A&B tuned fraction across time - TS - Recall D1 vs. D7'};

[t_ttest.recallAB_6_7] = paired_ttest_table_entry(data_input,...
        4, 'e', 'by animal', comp_descrip_in);


%% 2-way repeated measures mixed effects ANOVA (equivalent to PRISM output)

%PV A

PV_learn.A = source_data_short_learn_recall.PV.st_learn.A';
PV_recall.A = source_data_short_learn_recall.PV.st_recall.d4_d5_sub.A';

%format data for 2-way analysis
PV_recall.A(:,[4,5]) = nan;

%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = PV_learn.A(:,[2,3,6,7]);
data_in_2 = PV_recall.A(:,[2,3,6,7]);

%output stats for 2way RM LME analysis
lme_stats.PV_learn_recall.A = two_way_rm_lme(data_in_1,data_in_2);

%create formatting table for 2-way analysis
[t_2way_rm_lme.PV_learn_recallA] = two_way_lme_table_entry(4,'f','by animal',...
                'Population vector correlation relative to D1  - A laps',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.PV_learn_recall.A);       
            
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%PV day 6 learn vs recall
[ttest_stats.PV_learn_recall_6A] = unpaired_ttest(PV_learn.A(:,6), PV_recall.A(:,6));

%PV day 7 learn vs recall
[ttest_stats.PV_learn_recall_7A] = unpaired_ttest(PV_learn.A(:,7), PV_recall.A(:,7));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.PV_learn_recall_6A; ttest_stats.PV_learn_recall_7A];

comp_descrip_in = {'Learning vs recall PV correlation D6 - A laps';...
                'Learning vs recall PV correlation D7 - A laps'};

[t_ttest.learn_recallA_6_7] = unpaired_ttest_table_entry(data_input,...
        4, 'f', 'by animal', comp_descrip_in);

%paired t-test comparisons

%PV d2 vs. d6 time learn A 
ttest_stats.PV_learnA2v6 = paired_ttest(PV_learn.A(:,2), PV_learn.A(:,6));

%PV d2 vs d7 time learn A
ttest_stats.PV_learnA2v7 = paired_ttest(PV_learn.A(:,2), PV_learn.A(:,7));

%PV d2 vs. d6 time recall A 
ttest_stats.PV_recallA2v6 = paired_ttest(PV_recall.A(:,2), PV_recall.A(:,6));

%PV d2 vs d7 time recall A
ttest_stats.PV_recallA2v7 = paired_ttest(PV_recall.A(:,2), PV_recall.A(:,7));

%generate paired t-test table entries
%learn A time D2 vs. D7
data_input = [ttest_stats.PV_learnA2v6; ttest_stats.PV_learnA2v7];

%what comparisons are being made
comp_descrip_in = {'Learning PV correlation - time -D2 vs. D6 - A laps';...
                'Learning PV correlation - time -D2 vs. D7 - A laps'};

[t_ttest.PV_learnA2v6_7] = paired_ttest_table_entry(data_input,...
        4, 'f', 'by animal', comp_descrip_in);
    
%recall A time D2 vs. D7
data_input = [ttest_stats.PV_recallA2v6; ttest_stats.PV_recallA2v7];

%what comparisons are being made
comp_descrip_in = {'Recall PV correlation - time -D2 vs. D6 - A laps';...
                'Recall PV correlation - time -D2 vs. D7 - A laps'};

[t_ttest.PV_recallA2v6_7] = paired_ttest_table_entry(data_input,...
        4, 'f', 'by animal', comp_descrip_in);

    
%% PV Correlation B

PV_learn.B = source_data_short_learn_recall.PV.st_learn.B';
PV_recall.B = source_data_short_learn_recall.PV.st_recall.d4_d5_sub.B';

%format data for 2-way analysis
PV_recall.B(:,[4,5]) = nan;

%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = PV_learn.B(:,[2,3,6,7]);
data_in_2 = PV_recall.B(:,[2,3,6,7]);

%output stats for 2way RM LME analysis
lme_stats.PV_learn_recall.B = two_way_rm_lme(data_in_1,data_in_2);

%create formatting table for 2-way analysis
[t_2way_rm_lme.PV_learn_recallB] = two_way_lme_table_entry(4,'f','by animal',...
                'Population vector correlation relative to D1  - B laps',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.PV_learn_recall.B);       
            
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%PV day 6 learn vs recall
[ttest_stats.PV_learn_recall_6B] = unpaired_ttest(PV_learn.B(:,6), PV_recall.B(:,6));

%PV day 7 learn vs recall
[ttest_stats.PV_learn_recall_7B] = unpaired_ttest(PV_learn.B(:,7), PV_recall.B(:,7));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.PV_learn_recall_6B; ttest_stats.PV_learn_recall_7B];

comp_descrip_in = {'Learning vs recall PV correlation D6 - B laps';...
                'Learning vs recall PV correlation D7 - B laps'};

[t_ttest.learn_recallB_6_7] = unpaired_ttest_table_entry(data_input,...
        4, 'f', 'by animal', comp_descrip_in);

%paired t-test comparisons

%PV d2 vs. d6 time learn A 
ttest_stats.PV_learnB2v6 = paired_ttest(PV_learn.B(:,2), PV_learn.B(:,6));

%PV d2 vs d7 time learn A
ttest_stats.PV_learnB2v7 = paired_ttest(PV_learn.B(:,2), PV_learn.B(:,7));

%PV d2 vs. d6 time recall A 
ttest_stats.PV_recallB2v6 = paired_ttest(PV_recall.B(:,2), PV_recall.B(:,6));

%PV d2 vs d7 time recall A
ttest_stats.PV_recallB2v7 = paired_ttest(PV_recall.B(:,2), PV_recall.B(:,7));

%generate paired t-test table entries
%learn A time D2 vs. D7
data_input = [ttest_stats.PV_learnB2v6; ttest_stats.PV_learnB2v7];

%what comparisons are being made
comp_descrip_in = {'Learning PV correlation - time -D2 vs. D6 - B laps';...
                'Learning PV correlation - time -D2 vs. D7 - B laps'};

[t_ttest.PV_learnB2v6_7] = paired_ttest_table_entry(data_input,...
        4, 'f', 'by animal', comp_descrip_in);
    
%recall A time D2 vs. D7
data_input = [ttest_stats.PV_recallB2v6; ttest_stats.PV_recallB2v7];

%what comparisons are being made
comp_descrip_in = {'Recall PV correlation - time -D2 vs. D6 - B laps';...
                'Recall PV correlation - time -D2 vs. D7 - B laps'};

[t_ttest.PV_recallB2v6_7] = paired_ttest_table_entry(data_input,...
        4, 'f', 'by animal', comp_descrip_in);

%% TC T.S. Correlation A

TC_ts_learn.A = source_data_short_learn_recall.TC.ts.st_learn.A';
TC_ts_recall.A = source_data_short_learn_recall.TC.ts.st_recall.d4_d5_sub.A';

%format data for 2-way analysis
TC_ts_recall.A(:,[4,5]) = nan;

%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = TC_ts_learn.A(:,[2,3,6,7]);
data_in_2 = TC_ts_recall.A(:,[2,3,6,7]);

%output stats for 2way RM LME analysis
lme_stats.TC_ts_learn_recall.A = two_way_rm_lme(data_in_1,data_in_2);

%create formatting table for 2-way analysis
[t_2way_rm_lme.TC_ts_learn_recallA] = two_way_lme_table_entry(4,'f','by animal',...
                'Tuning curve TS correlation relative to D1  - A laps',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.TC_ts_learn_recall.A);       
            
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%PV day 6 learn vs recall
[ttest_stats.TC_ts_learn_recall_6A] = unpaired_ttest(TC_ts_learn.A(:,6), TC_ts_recall.A(:,6));

%PV day 7 learn vs recall
[ttest_stats.TC_ts_learn_recall_7A] = unpaired_ttest(TC_ts_learn.A(:,7), TC_ts_recall.A(:,7));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.TC_ts_learn_recall_6A; ttest_stats.TC_ts_learn_recall_7A];

comp_descrip_in = {'Learning vs recall Tuning curve TS correlation D6 - A laps';...
                'Learning vs recall Tuning curve TS correlation D7 - A laps'};

[t_ttest.TC_ts_learn_recallA_6_7] = unpaired_ttest_table_entry(data_input,...
        4, 'f', 'by animal', comp_descrip_in);

%paired t-test comparisons

%PV d2 vs. d6 time learn A 
ttest_stats.TC_ts_learnA2v6 = paired_ttest(TC_ts_learn.A(:,2), TC_ts_learn.A(:,6));

%PV d2 vs d7 time learn A
ttest_stats.TC_ts_learnA2v7 = paired_ttest(TC_ts_learn.A(:,2), TC_ts_learn.A(:,7));

%PV d2 vs. d6 time recall A 
ttest_stats.TC_ts_recallA2v6 = paired_ttest(TC_ts_recall.A(:,2), TC_ts_recall.A(:,6));

%PV d2 vs d7 time recall A
ttest_stats.TC_ts_recallA2v7 = paired_ttest(TC_ts_recall.A(:,2), TC_ts_recall.A(:,7));

%generate paired t-test table entries
%learn A time D2 vs. D7
data_input = [ttest_stats.TC_ts_learnA2v6; ttest_stats.TC_ts_learnA2v7];

%what comparisons are being made
comp_descrip_in = {'Learning Tuning curve TS correlation - time -D2 vs. D6 - A laps';...
                'Learning Tuning curve TS correlation - time -D2 vs. D7 - A laps'};

[t_ttest.TC_ts_learnA2v6_7] = paired_ttest_table_entry(data_input,...
        4, 'f', 'by animal', comp_descrip_in);
    
%recall A time D2 vs. D7
data_input = [ttest_stats.TC_ts_recallA2v6; ttest_stats.TC_ts_recallA2v7];

%what comparisons are being made
comp_descrip_in = {'Recall Tuning curve TS correlation - time -D2 vs. D6 - A laps';...
                'Recall Tuning curve TS correlation - time -D2 vs. D7 - A laps'};

[t_ttest.TC_ts_recallA2v6_7] = paired_ttest_table_entry(data_input,...
        4, 'f', 'by animal', comp_descrip_in);
    
%% TC T.S. Correlation B

TC_ts_learn.B = source_data_short_learn_recall.TC.ts.st_learn.B';
TC_ts_recall.B = source_data_short_learn_recall.TC.ts.st_recall.d4_d5_sub.B';

%format data for 2-way analysis
TC_ts_recall.B(:,[4,5]) = nan;

%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = TC_ts_learn.B(:,[2,3,6,7]);
data_in_2 = TC_ts_recall.B(:,[2,3,6,7]);

%output stats for 2way RM LME analysis
lme_stats.TC_ts_learn_recall.B = two_way_rm_lme(data_in_1,data_in_2);

%create formatting table for 2-way analysis
[t_2way_rm_lme.TC_ts_learn_recallB] = two_way_lme_table_entry(4,'f','by animal',...
                'Tuning curve TS correlation relative to D1  - B laps',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.TC_ts_learn_recall.B);       
            
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%PV day 6 learn vs recall
[ttest_stats.TC_ts_learn_recall_6B] = unpaired_ttest(TC_ts_learn.B(:,6), TC_ts_recall.B(:,6));

%PV day 7 learn vs recall
[ttest_stats.TC_ts_learn_recall_7B] = unpaired_ttest(TC_ts_learn.B(:,7), TC_ts_recall.B(:,7));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.TC_ts_learn_recall_6B; ttest_stats.TC_ts_learn_recall_7B];

comp_descrip_in = {'Learning vs recall Tuning curve TS correlation D6 - B laps';...
                'Learning vs recall Tuning curve TS correlation D7 - B laps'};

[t_ttest.TC_ts_learn_recallB_6_7] = unpaired_ttest_table_entry(data_input,...
        4, 'f', 'by animal', comp_descrip_in);

%paired t-test comparisons

%PV d2 vs. d6 time learn A 
ttest_stats.TC_ts_learnB2v6 = paired_ttest(TC_ts_learn.B(:,2), TC_ts_learn.B(:,6));

%PV d2 vs d7 time learn A
ttest_stats.TC_ts_learnB2v7 = paired_ttest(TC_ts_learn.B(:,2), TC_ts_learn.B(:,7));

%PV d2 vs. d6 time recall A 
ttest_stats.TC_ts_recallB2v6 = paired_ttest(TC_ts_recall.B(:,2), TC_ts_recall.B(:,6));

%PV d2 vs d7 time recall A
ttest_stats.TC_ts_recallB2v7 = paired_ttest(TC_ts_recall.B(:,2), TC_ts_recall.B(:,7));

%generate paired t-test table entries
%learn A time D2 vs. D7
data_input = [ttest_stats.TC_ts_learnB2v6; ttest_stats.TC_ts_learnB2v7];

%what comparisons are being made
comp_descrip_in = {'Learning Tuning curve TS correlation - time -D2 vs. D6 - B laps';...
                'Learning Tuning curve TS correlation - time -D2 vs. D7 - B laps'};

[t_ttest.TC_ts_learnB2v6_7] = paired_ttest_table_entry(data_input,...
        4, 'f', 'by animal', comp_descrip_in);
    
%recall A time D2 vs. D7
data_input = [ttest_stats.TC_ts_recallB2v6; ttest_stats.TC_ts_recallB2v7];

%what comparisons are being made
comp_descrip_in = {'Recall Tuning curve TS correlation - time -D2 vs. D6 - B laps';...
                'Recall Tuning curve TS correlation - time -D2 vs. D7 - B laps'};

[t_ttest.TC_ts_recallB2v6_7] = paired_ttest_table_entry(data_input,...
        4, 'f', 'by animal', comp_descrip_in);
    
%% PV neighboring days A
%PV neighbor A

PVn_learn.A = source_data_short_learn_recall.PV.neighbor.st_learn.A;
PVn_recall.A = source_data_short_learn_recall.PV.neighbor.st_recall.A;

%format data for 2-way analysis
%PVn_recall.A(:,[4,5]) = nan;

%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = PVn_learn.A(:,[1,2,6]);
data_in_2 = PVn_recall.A(:,[1,2,6]);

%output stats for 2way RM LME analysis
lme_stats.PVn_learn_recall.A = two_way_rm_lme(data_in_1,data_in_2);

%create formatting table for 2-way analysis
[t_2way_rm_lme.PVn_learn_recallA] = two_way_lme_table_entry(4,'g','by animal',...
                'Population vector correlation neighboring days - learning vs. recall - A, 1 vs. 2, 2 vs. 3, 6 vs. 7',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.PVn_learn_recall.A);       

%paired t-test comparisons
%PV neighbor corr 1 vs 2 vs. 6 vs 7 time learn A 
ttest_stats.PVn_learnA12v67 = paired_ttest(PVn_learn.A(:,1), PVn_learn.A(:,6));

%PV neighbor corr 1vs 2 vs. 6 vs 7 time recall A 
ttest_stats.PVn_recallA12v67 = paired_ttest(PVn_recall.A(:,1), PVn_recall.A(:,6));

%paired t-test single-entry table - PV neighbor Learn 1vs2 vs 6 vs 7
data_input = [ttest_stats.PVn_learnA12v67];
comp_descrip_in = {'PV correlation learning A - time - 1 vs. 2 vs. 6 vs. 7'};
[t_ttest.PVn_learnA12v67] = paired_ttest_single_table_entry(data_input,...
        4, 'g', 'by animal', comp_descrip_in);

%paired t-test single-entry table - PV neighbor recall 1vs2 vs 6 vs 7
data_input = [ttest_stats.PVn_recallA12v67];
comp_descrip_in = {'PV correlation recall A - time - 1 vs. 2 vs. 6 vs. 7'};
[t_ttest.PVn_recallA12v67] = paired_ttest_single_table_entry(data_input,...
        4, 'g', 'by animal', comp_descrip_in);
    
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%PV neighbor correlation day 1 vs. 2 learn vs recall
[ttest_stats.PVn_learn_recall_1v2A] = unpaired_ttest(PVn_learn.A(:,1), PVn_recall.A(:,1));

%PV  neighbor day 6 vs. 7 learn vs recall
[ttest_stats.PVn_learn_recall_6v7A] = unpaired_ttest(PVn_learn.A(:,6), PVn_recall.A(:,6));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.PVn_learn_recall_1v2A; ttest_stats.PVn_learn_recall_6v7A];

comp_descrip_in = {'PV correlation learning vs. recall A - behavior - 1 vs. 2';...
                'PV correlation learning vs. recall A - behavior -6 vs. 7'};

[t_ttest.PVn_learn_recallA_12_67] = unpaired_ttest_table_entry(data_input,...
        4, 'g', 'by animal', comp_descrip_in);

%% PV neighboring days B
%PV neighbor B

PVn_learn.B = source_data_short_learn_recall.PV.neighbor.st_learn.B;
PVn_recall.B = source_data_short_learn_recall.PV.neighbor.st_recall.B;

%format data for 2-way analysis
%PVn_recall.A(:,[4,5]) = nan;

%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = PVn_learn.B(:,[1,2,6]);
data_in_2 = PVn_recall.B(:,[1,2,6]);

%output stats for 2way RM LME analysis
lme_stats.PVn_learn_recall.B = two_way_rm_lme(data_in_1,data_in_2);

%create formatting table for 2-way analysis
[t_2way_rm_lme.PVn_learn_recallB] = two_way_lme_table_entry(4,'g','by animal',...
                'Population vector correlation neighboring days - learning vs. recall - B, 1 vs. 2, 2 vs. 3, 6 vs. 7',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.PVn_learn_recall.B);       

%paired t-test comparisons
%PV neighbor corr 1 vs 2 vs. 6 vs 7 time learn A 
ttest_stats.PVn_learnB12v67 = paired_ttest(PVn_learn.B(:,1), PVn_learn.B(:,6));

%PV neighbor corr 1vs 2 vs. 6 vs 7 time recall A 
ttest_stats.PVn_recallB12v67 = paired_ttest(PVn_recall.B(:,1), PVn_recall.B(:,6));

%paired t-test single-entry table - PV neighbor Learn 1vs2 vs 6 vs 7
data_input = [ttest_stats.PVn_learnB12v67];
comp_descrip_in = {'PV correlation learning B - time - 1 vs. 2 vs. 6 vs. 7'};
[t_ttest.PVn_learnB12v67] = paired_ttest_single_table_entry(data_input,...
        4, 'g', 'by animal', comp_descrip_in);

%paired t-test single-entry table - PV neighbor recall 1vs2 vs 6 vs 7
data_input = [ttest_stats.PVn_recallB12v67];
comp_descrip_in = {'PV correlation recall B - time - 1 vs. 2 vs. 6 vs. 7'};
[t_ttest.PVn_recallB12v67] = paired_ttest_single_table_entry(data_input,...
        4, 'g', 'by animal', comp_descrip_in);
    
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%PV neighbor correlation day 1 vs. 2 learn vs recall
[ttest_stats.PVn_learn_recall_1v2B] = unpaired_ttest(PVn_learn.B(:,1), PVn_recall.B(:,1));

%PV  neighbor day 6 vs. 7 learn vs recall
[ttest_stats.PVn_learn_recall_6v7B] = unpaired_ttest(PVn_learn.B(:,6), PVn_recall.B(:,6));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.PVn_learn_recall_1v2B; ttest_stats.PVn_learn_recall_6v7B];

comp_descrip_in = {'PV correlation learning vs. recall B - behavior - 1 vs. 2';...
                'PV correlation learning vs. recall B - behavior -6 vs. 7'};

[t_ttest.PVn_learn_recallB_12_67] = unpaired_ttest_table_entry(data_input,...
        4, 'g', 'by animal', comp_descrip_in);
%% TC TS neighboring days A

%TC TS neighbor A

TC_ts_n_learn.A = source_data_short_learn_recall.TC.neighbor.ts.st_learn.A;
TC_ts_n_recall.A = source_data_short_learn_recall.TC.neighbor.ts.st_recall.A;

%format data for 2-way analysis
%PVn_recall.A(:,[4,5]) = nan;

%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = TC_ts_n_learn.A(:,[1,2,6]);
data_in_2 = TC_ts_n_recall.A(:,[1,2,6]);

%output stats for 2way RM LME analysis
lme_stats.TC_ts_n_learn_recall.A = two_way_rm_lme(data_in_1,data_in_2);

%%%

%create formatting table for 2-way analysis
[t_2way_rm_lme.TC_ts_n_learn_recallA] = two_way_lme_table_entry(4,'g','by animal',...
                'Tuning curve (TS) correlation neighboring days - learning vs. recall - A, 1 vs. 2, 2 vs. 3, 6 vs. 7',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.TC_ts_n_learn_recall.A);       

%paired t-test comparisons
%PV neighbor corr 1 vs 2 vs. 6 vs 7 time learn A 
ttest_stats.TC_ts_n_learnA12v67 = paired_ttest(TC_ts_n_learn.A(:,1), TC_ts_n_learn.A(:,6));

%PV neighbor corr 1vs 2 vs. 6 vs 7 time recall A 
ttest_stats.TC_ts_n_recallA12v67 = paired_ttest(TC_ts_n_recall.A(:,1), TC_ts_n_recall.A(:,6));

%paired t-test single-entry table - PV neighbor Learn 1vs2 vs 6 vs 7
data_input = [ttest_stats.TC_ts_n_learnA12v67];
comp_descrip_in = {'TC correlation learning A - time - 1 vs. 2 vs. 6 vs. 7'};
[t_ttest.TC_ts_n_learnA12v67] = paired_ttest_single_table_entry(data_input,...
        4, 'g', 'by animal', comp_descrip_in);

%paired t-test single-entry table - PV neighbor recall 1vs2 vs 6 vs 7
data_input = [ttest_stats.TC_ts_n_recallA12v67];
comp_descrip_in = {'TC correlation recall A - time - 1 vs. 2 vs. 6 vs. 7'};
[t_ttest.TC_ts_n_recallA12v67] = paired_ttest_single_table_entry(data_input,...
        4, 'g', 'by animal', comp_descrip_in);
    
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%PV neighbor correlation day 1 vs. 2 learn vs recall
[ttest_stats.TC_ts_n_learn_recall_1v2A] = unpaired_ttest(TC_ts_n_learn.A(:,1), TC_ts_n_recall.A(:,1));

%PV  neighbor day 6 vs. 7 learn vs recall
[ttest_stats.TC_ts_n_learn_recall_6v7A] = unpaired_ttest(TC_ts_n_learn.A(:,6), TC_ts_n_recall.A(:,6));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.TC_ts_n_learn_recall_1v2A; ttest_stats.TC_ts_n_learn_recall_6v7A];

comp_descrip_in = {'TC correlation learning vs. recall A - behavior - 1 vs. 2 ';...
                'TC correlation learning vs. recall A - behavior -6 vs. 7 '};

[t_ttest.TC_ts_n_learn_recallA_12_67] = unpaired_ttest_table_entry(data_input,...
        4, 'g', 'by animal', comp_descrip_in);


%% TC TS neighboring days B

%TC TS neighbor B

TC_ts_n_learn.B = source_data_short_learn_recall.TC.neighbor.ts.st_learn.B;
TC_ts_n_recall.B = source_data_short_learn_recall.TC.neighbor.ts.st_recall.B;

%format data for 2-way analysis
%PVn_recall.A(:,[4,5]) = nan;

%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = TC_ts_n_learn.B(:,[1,2,6]);
data_in_2 = TC_ts_n_recall.B(:,[1,2,6]);

%output stats for 2way RM LME analysis
lme_stats.TC_ts_n_learn_recall.B = two_way_rm_lme(data_in_1,data_in_2);

%%%

%create formatting table for 2-way analysis
[t_2way_rm_lme.TC_ts_n_learn_recallB] = two_way_lme_table_entry(4,'g','by animal',...
                'Tuning curve (TS) correlation neighboring days - learning vs. recall - B, 1 vs. 2, 2 vs. 3, 6 vs. 7',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.TC_ts_n_learn_recall.B);       

%paired t-test comparisons
%PV neighbor corr 1 vs 2 vs. 6 vs 7 time learn A 
ttest_stats.TC_ts_n_learnB12v67 = paired_ttest(TC_ts_n_learn.B(:,1), TC_ts_n_learn.B(:,6));

%PV neighbor corr 1vs 2 vs. 6 vs 7 time recall A 
ttest_stats.TC_ts_n_recallB12v67 = paired_ttest(TC_ts_n_recall.B(:,1), TC_ts_n_recall.B(:,6));

%paired t-test single-entry table - PV neighbor Learn 1vs2 vs 6 vs 7
data_input = [ttest_stats.TC_ts_n_learnB12v67];
comp_descrip_in = {'TC correlation learning B - time - 1 vs. 2 vs. 6 vs. 7'};
[t_ttest.TC_ts_n_learnB12v67] = paired_ttest_single_table_entry(data_input,...
        4, 'g', 'by animal', comp_descrip_in);

%paired t-test single-entry table - PV neighbor recall 1vs2 vs 6 vs 7
data_input = [ttest_stats.TC_ts_n_recallB12v67];
comp_descrip_in = {'TC correlation recall B - time - 1 vs. 2 vs. 6 vs. 7'};
[t_ttest.TC_ts_n_recallB12v67] = paired_ttest_single_table_entry(data_input,...
        4, 'g', 'by animal', comp_descrip_in);
    
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%PV neighbor correlation day 1 vs. 2 learn vs recall
[ttest_stats.TC_ts_n_learn_recall_1v2B] = unpaired_ttest(TC_ts_n_learn.B(:,1), TC_ts_n_recall.B(:,1));

%PV  neighbor day 6 vs. 7 learn vs recall
[ttest_stats.TC_ts_n_learn_recall_6v7B] = unpaired_ttest(TC_ts_n_learn.B(:,6), TC_ts_n_recall.B(:,6));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.TC_ts_n_learn_recall_1v2B; ttest_stats.TC_ts_n_learn_recall_6v7B];

comp_descrip_in = {'TC correlation learning vs. recall B - behavior - 1 vs. 2 ';...
                'TC correlation learning vs. recall B - behavior -6 vs. 7 '};

[t_ttest.TC_ts_n_learn_recallB_12_67] = unpaired_ttest_table_entry(data_input,...
        4, 'g', 'by animal', comp_descrip_in);
    
%% Centroid TS distance difference A
%Centroid TS neighbor A

cent_ts_learn.A = source_data_short_learn_recall.angle_diff.ts.st_learn.animal.mean.A;
%same data as in raw sub-struct
cent_ts_recall.A = source_data_short_learn_recall.angle_diff.ts.st_recall.animal.mean.A;


%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = cent_ts_learn.A(:,[2,3,6,7]);
data_in_2 = cent_ts_recall.A(:,[2,3,6,7]);

%output stats for 2way RM LME analysis
lme_stats.cent_ts_learn_recall.A = two_way_rm_lme(data_in_1,data_in_2);

%create formatting table for 2-way analysis
[t_2way_rm_lme.cent_ts_learn_recallA] = two_way_lme_table_entry(4,'h','by animal',...
                'Centroid difference relative to d1 - learning vs. recall - A - rel d1, rel d2, rel d5, rel d6',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.cent_ts_learn_recall.A);       

%paired t-test comparisons (learn)
%rel d1 vs d5 learning
ttest_stats.cent_ts_learnA1v5 = paired_ttest(cent_ts_learn.A(:,2), cent_ts_learn.A(:,6));

%rel d1 vs d6 learning
ttest_stats.cent_ts_learnA1v6 = paired_ttest(cent_ts_learn.A(:,2), cent_ts_learn.A(:,7));

%paired t-test single-entry table -  Learn rel d1,d5 and rel d1,d6
data_input = [ttest_stats.cent_ts_learnA1v5; ttest_stats.cent_ts_learnA1v6];
comp_descrip_in = {'Centroid difference relative to d1 - learning  - A - rel d1 vs. rel d5';...
                    'Centroid difference relative to d1 - learning  - A - rel d1 vs. rel d6'};
[t_ttest.cent_ts_learnA15v16] = paired_ttest_table_entry(data_input,...
        4, 'h', 'by animal', comp_descrip_in);

%paired t-test comparisons (recall)
%rel d1 vs d5 recall
ttest_stats.cent_ts_recallA1v5 = paired_ttest(cent_ts_recall.A(:,2), cent_ts_recall.A(:,6));

%rel d1 vs d6 recall
ttest_stats.cent_ts_recallA1v6 = paired_ttest(cent_ts_recall.A(:,2), cent_ts_recall.A(:,7));

%paired t-test single-entry table -  Learn rel d1,d5 and rel d1,d6
data_input = [ttest_stats.cent_ts_recallA1v5; ttest_stats.cent_ts_recallA1v6];
comp_descrip_in = {'Centroid difference relative to d1 - recall  - A - rel d1 vs. rel d5';...
                    'Centroid difference relative to d1 - recall  - A - rel d1 vs. rel d6'};
[t_ttest.cent_ts_recallA15v16] = paired_ttest_table_entry(data_input,...
        4, 'h', 'by animal', comp_descrip_in);
    
    
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%cent diff learning vs recall A - d5
[ttest_stats.cent_ts_learn_recall_15A] = unpaired_ttest(cent_ts_learn.A(:,6), cent_ts_recall.A(:,6));

%cent diff learning vs recall A - d6
[ttest_stats.cent_ts_learn_recall_16A] = unpaired_ttest(cent_ts_learn.A(:,7), cent_ts_recall.A(:,7));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.cent_ts_learn_recall_15A; ttest_stats.cent_ts_learn_recall_16A];

comp_descrip_in = {'Centroid difference relative to d1 - learning vs. recall  - A - rel d5';...
                'Centroid difference relative to d1 - learning vs. recall  - A - rel d6'};

[t_ttest.cent_ts_learn_recallA_15_16] = unpaired_ttest_table_entry(data_input,...
        4, 'h', 'by animal', comp_descrip_in);


%% Centroid TS distance difference B
%Centroid TS neighbor B

cent_ts_learn.B = source_data_short_learn_recall.angle_diff.ts.st_learn.animal.mean.B;
%same data as in raw sub-struct
cent_ts_recall.B = source_data_short_learn_recall.angle_diff.ts.st_recall.animal.mean.B;


%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = cent_ts_learn.B(:,[2,3,6,7]);
data_in_2 = cent_ts_recall.B(:,[2,3,6,7]);

%output stats for 2way RM LME analysis
lme_stats.cent_ts_learn_recall.B = two_way_rm_lme(data_in_1,data_in_2);

%create formatting table for 2-way analysis
[t_2way_rm_lme.cent_ts_learn_recallB] = two_way_lme_table_entry(4,'h','by animal',...
                'Centroid difference relative to d1 - learning vs. recall - B - rel d1, rel d2, rel d5, rel d6',...
                [size(data_in_1,1), size(data_in_2,1)],lme_stats.cent_ts_learn_recall.B);       

%paired t-test comparisons (learn)
%rel d1 vs d5 learning
ttest_stats.cent_ts_learnB1v5 = paired_ttest(cent_ts_learn.B(:,2), cent_ts_learn.B(:,6));

%rel d1 vs d6 learning
ttest_stats.cent_ts_learnB1v6 = paired_ttest(cent_ts_learn.B(:,2), cent_ts_learn.B(:,7));

%paired t-test single-entry table -  Learn rel d1,d5 and rel d1,d6
data_input = [ttest_stats.cent_ts_learnB1v5; ttest_stats.cent_ts_learnB1v6];
comp_descrip_in = {'Centroid difference relative to d1 - learning  - B - rel d1 vs. rel d5';...
                    'Centroid difference relative to d1 - learning  - B - rel d1 vs. rel d6'};
[t_ttest.cent_ts_learnB15v16] = paired_ttest_table_entry(data_input,...
        4, 'h', 'by animal', comp_descrip_in);

%paired t-test comparisons (recall)
%rel d1 vs d5 recall
ttest_stats.cent_ts_recallB1v5 = paired_ttest(cent_ts_recall.B(:,2), cent_ts_recall.B(:,6));

%rel d1 vs d6 recall
ttest_stats.cent_ts_recallB1v6 = paired_ttest(cent_ts_recall.B(:,2), cent_ts_recall.B(:,7));

%paired t-test single-entry table -  Learn rel d1,d5 and rel d1,d6
data_input = [ttest_stats.cent_ts_recallB1v5; ttest_stats.cent_ts_recallB1v6];
comp_descrip_in = {'Centroid difference relative to d1 - recall  - B - rel d1 vs. rel d5';...
                    'Centroid difference relative to d1 - recall  - B - rel d1 vs. rel d6'};
[t_ttest.cent_ts_recallB15v16] = paired_ttest_table_entry(data_input,...
        4, 'h', 'by animal', comp_descrip_in);
    
    
%unpaired t-tests
%output: p-val, t-statistic, n1, n2, dof
%cent diff learning vs recall A - d5
[ttest_stats.cent_ts_learn_recall_15B] = unpaired_ttest(cent_ts_learn.B(:,6), cent_ts_recall.B(:,6));

%cent diff learning vs recall A - d6
[ttest_stats.cent_ts_learn_recall_16B] = unpaired_ttest(cent_ts_learn.B(:,7), cent_ts_recall.B(:,7));

%create formatting table for unpaired t-test comparison
data_input = [ttest_stats.cent_ts_learn_recall_15B; ttest_stats.cent_ts_learn_recall_16B];

comp_descrip_in = {'Centroid difference relative to d1 - learning vs. recall  - B - rel d5';...
                'Centroid difference relative to d1 - learning vs. recall  - B - rel d6'};

[t_ttest.cent_ts_learn_recallB_15_16] = unpaired_ttest_table_entry(data_input,...
        4, 'h', 'by animal', comp_descrip_in); 
    
%% Assemble Figure 4 stats export table
t1 = repmat({' '},1,12);
blank_row = cell2table(t1);

%spreadsheet name
spreadsheet_name = 'statistics_summary.xlsx';

%sheet name
sheet_name = 'Figure 4';

%exported Excel spreadsheet
%write to Excel spreadsheet
%learn 4e
insert_table_rows(t_1_rm_lme.learnA,spreadsheet_name,sheet_name,'overwritesheet')
insert_table_rows(t_ttest.learnA_6_7,spreadsheet_name,sheet_name,'append')

insert_table_rows(blank_row,spreadsheet_name,sheet_name,'append')
insert_table_rows(t_1_rm_lme.learnB,spreadsheet_name,sheet_name,'append')
insert_table_rows(t_ttest.learnB_6_7,spreadsheet_name,sheet_name,'append')

insert_table_rows(blank_row,spreadsheet_name,sheet_name,'append')
insert_table_rows(t_1_rm_lme.learnAB,spreadsheet_name,sheet_name,'append')
insert_table_rows(t_ttest.learnAB_6_7,spreadsheet_name,sheet_name,'append')

insert_table_rows(blank_row,spreadsheet_name,sheet_name,'append')

%recall 4e
insert_table_rows(t_1_rm_lme.recallA,spreadsheet_name,sheet_name,'append')
insert_table_rows(t_ttest.recallA_6_7,spreadsheet_name,sheet_name,'append')

insert_table_rows(blank_row,spreadsheet_name,sheet_name,'append')
insert_table_rows(t_1_rm_lme.recallB,spreadsheet_name,sheet_name,'append')
insert_table_rows(t_ttest.recallB_6_7,spreadsheet_name,sheet_name,'append')

insert_table_rows(blank_row,spreadsheet_name,sheet_name,'append')
insert_table_rows(t_1_rm_lme.recallAB,spreadsheet_name,sheet_name,'append');
insert_table_rows(t_ttest.recallAB_6_7,spreadsheet_name,sheet_name,'append');

%4f PV correlation across days -A trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.PV_learn_recallA,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.learn_recallA_6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.PV_learnA2v6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.PV_recallA2v6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%4f PV correlation across days -B trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.PV_learn_recallB,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.learn_recallB_6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.PV_learnB2v6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.PV_recallB2v6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%4f TC TS correlation across days - A trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.TC_ts_learn_recallA,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_ts_learn_recallA_6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_ts_learnA2v6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_ts_recallA2v6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%4f TC TS correlation across days -B trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.TC_ts_learn_recallB,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_ts_learn_recallB_6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_ts_learnB2v6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_ts_recallB2v6_7,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%4g PV neighbor correlation across days - A trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.PVn_learn_recallA,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.PVn_learnA12v67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.PVn_recallA12v67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.PVn_learn_recallA_12_67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%4g PV neighbor correlation across days - B trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.PVn_learn_recallB,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.PVn_learnB12v67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.PVn_recallB12v67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.PVn_learn_recallB_12_67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%4g TC TS neighbor correlation across days - A trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.TC_ts_n_learn_recallA,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_ts_n_learnA12v67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_ts_n_recallA12v67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_ts_n_learn_recallA_12_67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%4g TC TS neighbor correlation across days - B trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.TC_ts_n_learn_recallB,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_ts_n_learnB12v67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_ts_n_recallB12v67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.TC_ts_n_learn_recallB_12_67,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%4h Centroid difference TS neighbor correlation across days - A trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.cent_ts_learn_recallA,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.cent_ts_learnA15v16,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.cent_ts_recallA15v16,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.cent_ts_learn_recallA_15_16,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')

%4h Centroid difference TS neighbor correlation across days - B trials
%tables: 
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_2way_rm_lme.cent_ts_learn_recallB,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.cent_ts_learnB15v16,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.cent_ts_recallB15v16,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')
writetable(t_ttest.cent_ts_learn_recallB_15_16,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true','WriteMode','append')


%% Figure 5 Export statistics - short term learn TS 


% A&B correlation between A and B laps across time

%import TS data for each day
AB_ts_corr_learn = source_data_short_learn_recall.corr_analysis.st_learn.AB_corr.ts.raw.pooled.AB_corr_ratio;

%remove nan values from each column of input data
for ii=1:numel(AB_ts_corr_learn)
    temp_idx = find(isnan(AB_ts_corr_learn{ii}) == 1);
    if ~isempty(temp_idx)
        AB_ts_corr_learn{ii}(temp_idx) = [];
    end
end

%input data for stats analysis - choose which days to use to KW test
data_input = AB_ts_corr_learn(2:7)';

%run KW test
[kruskal_stats_learn] = kruskalwallis_test(data_input);

%generate single row entry for KW test
comp_descrip_in = {'TS - pooled neuron - group difference - Learn - D2,D3,D4,D5,D6,D7'};
[t_krusall_learn] = kruskal_wallis_single_table_entry(5,'c','pooled',...
                comp_descrip_in,kruskal_stats_learn);

%for 1 sample Wilcoxon tests
for ii=2:7
    paired_wilcoxon_stats_learn{ii} = paired_wilcoxon_signrank(AB_ts_corr_learn{ii},1);
end

%test comparison description for each table entry \
comp_descrip_in = {'AB Corr vs. 1, Day 2 - Learn';...
                    'AB Corr vs. 1, Day 3 - Learn';...
                    'AB Corr vs. 1, Day 4 - Learn';...
                    'AB Corr vs. 1, Day 5 - Learn';...
                    'AB Corr vs. 1, Day 6 - Learn';...
                    'AB Corr vs. 1, Day 7 - Learn'};

%table entries for paired wilcoxon test learning
paired_wilcox_AB_learn = paired_wilcoxon_table_entry_no_adj(5,'c','pooled',...
                comp_descrip_in,paired_wilcoxon_stats_learn(2:7));

%% Figure 5 Export statistics - short term recall TS 
% A&B correlation between A and B laps across time

%recall
AB_ts_corr_recall = source_data_short_learn_recall.corr_analysis.st_recall.AB_corr.ts.raw.pooled.AB_corr_ratio;

%remove nan values from each column of input data
for ii=1:numel(AB_ts_corr_recall)
    temp_idx = find(isnan(AB_ts_corr_recall{ii}) == 1);
    if ~isempty(temp_idx)
        AB_ts_corr_recall{ii}(temp_idx) = [];
    end
end

%input data for stats analysis - choose which days to use to KW test
data_input = AB_ts_corr_recall([2,3,6,7,8,9])';

%run KW test
[kruskal_stats_recall] = kruskalwallis_test(data_input);

%generate single row entry for KW test
comp_descrip_in = {'TS - pooled neuron - group difference - Recall - D2,D3,D6,D7,D8,D9'};
[t_krusall_recall] = kruskal_wallis_single_table_entry(5,'d','pooled',...
                comp_descrip_in,kruskal_stats_recall);

%for 1 sample Wilcoxon tests
for ii=[2,3,6,7,8,9]
    paired_wilcoxon_stats_recall{ii} = paired_wilcoxon_signrank(AB_ts_corr_recall{ii},1);
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
paired_wilcox_AB_recall = paired_wilcoxon_table_entry_no_adj(5,'d','pooled',...
                comp_descrip_in,paired_wilcoxon_stats_recall([2,3,6,7,8,9]));

%% Create figure 5 stats export spreadsheet
%spreadsheet name
spreadsheet_name = 'statistics_summary.xlsx';

%sheet name
sheet_name = 'Figure 5';

%empty row
t1 = repmat({' '},1,12);

%exported Excel spreadsheet
%write to Excel spreadsheet
%learn 5c
writetable(t_krusall_learn,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','overwritesheet')
writetable(paired_wilcox_AB_learn,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
%recall 5c
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_krusall_recall,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(paired_wilcox_AB_recall,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
            
%% 1-way Repeated Measures ANOVA - deal with this later (first figure - do with GG correction - works with MATLAB exchange function


%rm_sample_data = [];

% training_groups = repmat([1:4],4,1);
% training_groups = training_groups(:);
% 
% rm_table = table(rm_sample_data(:),training_groups,'VariableNames', {'score','t_grps'});
% 
% 
% rm_model = fitrm(rm_table,'score~t_grps','WithinModel','separatemeans')
% 
% rm_model = fitrm(rm_table,'score~t_grps','WithinDesign',table(rm_sample_data'))
% 
% rm_anova = ranova(rm_model)
% 
% anova_rm(rm_sample_data)

%repeated measures 1-way ANOVA
%paired t-test

%unpaired t-test
%kruskall wallis test
%write up functions for each test
%import all matlab data here

