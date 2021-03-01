%% Master script for generating stats for paper
%source data, import and format Prism analysis, and export to Excel, and
%Word legend data

%import data for each figure

%figure 2 source data
%load('G:\Google_drive\task_selective_place_paper\matlab_data\source_data_fig4_5_and_sup.mat')

%laptop directory
load('C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data\source_data_fig4_5_and_sup.mat')


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
%test_data  = [];

PV_learn.A = source_data_short_learn_recall.PV.st_learn.A';
PV_recall.A = source_data_short_learn_recall.PV.st_recall.d4_d5_sub.A';

%format data for 2-way analysis
PV_recall.A(:,[4,5]) = nan;

%data input (D2, D3, D6, D7 matching for learning and recall)
%assemble learning and recall matrices seperate and then concatenate

data_in_1 = PV_learn.A(:,[2,3,6,7]);
data_in_2= PV_recall.A(:,[2,3,6,7]);



%corr score vector
temp = test_data';
temp = temp(:);
corr_score = temp;

%subject data
temp = categorical(repmat(1:11,6,1))';
subject = temp(:);

%trial_type vector
trial_type = 3*ones(6,11);
trial_type(:,1:6) = 2;
trial_type = trial_type';
trial_type = categorical(trial_type(:)); 

%time vector
%define ordinal time rankining
time_value = {'1','2','3','4','5','6'};


time_vec = repmat([1:6]',1,11)';
time_vec = categorical(categorical(time_vec(:)),time_value,'Ordinal',true);

% %organize into columns for input to lme
% Corr_score
% Subject 
% Time
% Trial_type

%organize test data into table
lme_tbl = table(corr_score, subject, time_vec, trial_type, 'VariableNames',{'corr', 'subject','time','trial'});

%fit LME with subjects being random factor
lme = fitlme(lme_tbl,'corr ~ 1+ time*trial + (1|subject)','FitMethod','REML','DummyVarCoding','effects');

%fit LME with symmetrical covariance matrix (sphericity adjustment) - no
%sig difference
%compound symmetry structure/isotropic symmetry structure 'CompSymm'
%lme = fitlme(lme_tbl,'corr ~ 1+ time*trial + (1|subject)','FitMethod','REML','DummyVarCoding','effects');

%anova on lme
lme_stats = anova(lme,'DFMethod','satterthwaite')%% Figure 4f PV correlation



%% Assemble Figure 4 stats export table
t1 = repmat({' '},1,12);

%spreadsheet name
spreadsheet_name = 'statistics_summary.xlsx';

%sheet name
sheet_name = 'Figure 4';

%exported Excel spreadsheet
%write to Excel spreadsheet
%learn 4e
writetable(t_1_rm_lme.learnA,spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','overwrite')
writetable(t_ttest.learnA_6_7,spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')
writetable(t_1_rm_lme.learnB,spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')
writetable(t_ttest.learnB_6_7,spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')
writetable(t_1_rm_lme.learnAB,spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')
writetable(t_ttest.learnAB_6_7,spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')

%recall 4e
writetable(t_1_rm_lme.recallA,spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')
writetable(t_ttest.recallA_6_7,spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')
writetable(t_1_rm_lme.recallB,spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')
writetable(t_ttest.recallB_6_7,spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')

writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')
writetable(t_1_rm_lme.recallAB,spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')
writetable(t_ttest.recallAB_6_7,spreadsheet_name,'Sheet',sheet_name,"AutoFitWidth",true,'WriteMode','append')




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

