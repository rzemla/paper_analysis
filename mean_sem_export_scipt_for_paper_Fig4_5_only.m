%% import prism fraction of licks data on A and B trials

%% import data for figure 4 stats

dirpath = 'C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data';

load(fullfile(dirpath,'source_data_fig4_5_and_sup.mat'))
load(fullfile(dirpath,'fig4_5_table_data.mat'))

%% Mean and sem data for PV 
%%% All PV data A
%PV A learn 2 vs 7
[mean_2_7.PV.learn.A,sem_2_7.PV.learn.A,nb_comp_2_7.PV.learn.A] = learn_2_7_mean_sem_PV(source_data_short_learn_recall.PV.st_learn.A);

%PV B learn 2 vs 7
[mean_2_7.PV.learn.B,sem_2_7.PV.learn.B,nb_comp_2_7.PV.learn.B] = learn_2_7_mean_sem_PV(source_data_short_learn_recall.PV.st_learn.B);

%PV A recall 2 vs 7
[mean_2_7.PV.recall.A,sem_2_7.PV.recall.A,nb_comp_2_7.PV.recall.A] = recall_2_7_mean_sem_PV(source_data_short_learn_recall.PV.st_recall.d4_d5_sub.A);

%PV B recall 2 vs 7
[mean_2_7.PV.recall.B,sem_2_7.PV.recall.B,nb_comp_2_7.PV.recall.B] = recall_2_7_mean_sem_PV(source_data_short_learn_recall.PV.st_recall.d4_d5_sub.B);

%PV learn vs. recall A
%learn data
s1 = source_data_short_learn_recall.PV.st_learn.A;
%recall data
s2 = source_data_short_learn_recall.PV.st_recall.d4_d5_sub.A;

%mean and sem comp - learn vs. recall A
[mean_7_lr.PV.A, sem_7_lr.PV.A,nb_comp_7.PV.A] = learn_recall_7_mean_sem_PV(s1,s2);

%PV learn vs. recall B
%learn data
s1 = source_data_short_learn_recall.PV.st_learn.B;
%recall data
s2 = source_data_short_learn_recall.PV.st_recall.d4_d5_sub.B;

%mean and sem comp - learn vs. recall B
[mean_7_lr.PV.B, sem_7_lr.PV.B,nb_comp_7.PV.B] = learn_recall_7_mean_sem_PV(s1,s2);

%% Mean and sem data for TC TS (PV function works the same for TC corr values) 
%%% All TS data A
%TS A learn 2 vs 7
[mean_2_7.TC.ts.learn.A,sem_2_7.TC.ts.learn.A,nb_comp_2_7.TC.ts.learn.A] = learn_2_7_mean_sem_PV(source_data_short_learn_recall.TC.ts.st_learn.A);

%TC.ts B learn 2 vs 7
[mean_2_7.TC.ts.learn.B,sem_2_7.TC.ts.learn.B,nb_comp_2_7.TC.ts.learn.B] = learn_2_7_mean_sem_PV(source_data_short_learn_recall.TC.ts.st_learn.B);

%TC.ts A recall 2 vs 7
[mean_2_7.TC.ts.recall.A,sem_2_7.TC.ts.recall.A,nb_comp_2_7.TC.ts.recall.A] = recall_2_7_mean_sem_PV(source_data_short_learn_recall.TC.ts.st_recall.d4_d5_sub.A);

%TC.ts B recall 2 vs 7
[mean_2_7.TC.ts.recall.B,sem_2_7.TC.ts.recall.B,nb_comp_2_7.TC.ts.recall.B] = recall_2_7_mean_sem_PV(source_data_short_learn_recall.TC.ts.st_recall.d4_d5_sub.B);

%TC.ts learn vs. recall A
%learn data
s1 = source_data_short_learn_recall.TC.ts.st_learn.A;
%recall data
s2 = source_data_short_learn_recall.TC.ts.st_recall.d4_d5_sub.A;

%mean and sem comp - learn vs. recall A
[mean_7_lr.TC.ts.A, sem_7_lr.TC.ts.A,nb_comp_7.TC.ts.A] = learn_recall_7_mean_sem_PV(s1,s2);

%TC.ts learn vs. recall B
%learn data
s1 = source_data_short_learn_recall.TC.ts.st_learn.B;
%recall data
s2 = source_data_short_learn_recall.TC.ts.st_recall.d4_d5_sub.B;

%mean and sem comp - learn vs. recall B
[mean_7_lr.TC.ts.B, sem_7_lr.TC.ts.B,nb_comp_7.TC.ts.B] = learn_recall_7_mean_sem_PV(s1,s2);

%% Mean and sem data for TC SI (PV function works the same for TC corr values) 

%%% All TC data A
%TC A learn 2 vs 7
[mean_2_7.TC.si.learn.A,sem_2_7.TC.si.learn.A,nb_comp_2_7.TC.si.learn.A] = learn_2_7_mean_sem_PV(source_data_short_learn_recall.TC.si.st_learn.A);

%TC.si B learn 2 vs 7
[mean_2_7.TC.si.learn.B,sem_2_7.TC.si.learn.B,nb_comp_2_7.TC.si.learn.B] = learn_2_7_mean_sem_PV(source_data_short_learn_recall.TC.si.st_learn.B);

%TC.si A recall 2 vs 7
[mean_2_7.TC.si.recall.A,sem_2_7.TC.si.recall.A,nb_comp_2_7.TC.si.recall.A] = recall_2_7_mean_sem_PV(source_data_short_learn_recall.TC.si.st_recall.d4_d5_sub.A);

%TC.si B recall 2 vs 7
[mean_2_7.TC.si.recall.B,sem_2_7.TC.si.recall.B,nb_comp_2_7.TC.si.recall.B] = recall_2_7_mean_sem_PV(source_data_short_learn_recall.TC.si.st_recall.d4_d5_sub.B);

%TC.si learn vs. recall A
%learn data
s1 = source_data_short_learn_recall.TC.si.st_learn.A;
%recall data
s2 = source_data_short_learn_recall.TC.si.st_recall.d4_d5_sub.A;

%mean and sem comp - learn vs. recall A
[mean_7_lr.TC.si.A, sem_7_lr.TC.si.A,nb_comp_7.TC.si.A] = learn_recall_7_mean_sem_PV(s1,s2);

%TC.si learn vs. recall B
%learn data
s1 = source_data_short_learn_recall.TC.si.st_learn.B;
%recall data
s2 = source_data_short_learn_recall.TC.si.st_recall.d4_d5_sub.B;

%mean and sem comp - learn vs. recall B
[mean_7_lr.TC.si.B, sem_7_lr.TC.si.B,nb_comp_7.TC.si.B] = learn_recall_7_mean_sem_PV(s1,s2);

%% Mean and sem for neighbor analysis - PV and TC

%nb_animal x day comp score
%1 column - 1 vs. 2 ; 6th column - 6th vs. 7th day

%PV neigh learn A 1,2 vs 6,7
nei_mat = source_data_short_learn_recall.PV.neighbor.st_learn.A;
[mean_out, sem_out, nb_animals] = neighbor_corr_mean_sem(nei_mat);




%% Export the data to word

%% Start Word document that will contain the formatted stats data

WordFileName='fig_4_5_mean_sem.doc';
CurDir=pwd;
FileSpec = fullfile(CurDir,WordFileName);
%active X handle for manipulating document (ActXWord)
[ActXWord,WordHandle]=StartWord(FileSpec);

fprintf('Document will be saved in %s\n',FileSpec);

%% PV learning

%A 2 7 learning
txt_input = 'PV Day 2 vs. 7 A trial learning';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_2_7.PV.learn.A(1);
sem_txt = sem_2_7.PV.learn.A(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_2_7.PV.learn.A(2);
sem_txt = sem_2_7.PV.learn.A(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%number of samples
txt_input = ['nb samp ', num2str(nb_comp_2_7.PV.learn.A)];
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);

%B 2 7 learning
txt_input = 'PV Day 2 vs. 7 B trial learning';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);

mean_txt = mean_2_7.PV.learn.B(1);
sem_txt = sem_2_7.PV.learn.B(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

mean_txt = mean_2_7.PV.learn.B(2);
sem_txt = sem_2_7.PV.learn.B(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%number of samples
txt_input = ['nb samp ', num2str(nb_comp_2_7.PV.learn.B)];
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%% PV recall
writeWordEnter(ActXWord,WordHandle,1);

%A 2 7 recalling
txt_input = 'PV Day 2 vs. 7 A recall recalling';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_2_7.PV.recall.A(1);
sem_txt = sem_2_7.PV.recall.A(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_2_7.PV.recall.A(2);
sem_txt = sem_2_7.PV.recall.A(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%number of samples
txt_input = ['nb samp ', num2str(nb_comp_2_7.PV.recall.A)];
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);

%B 2 7 recalling
txt_input = 'PV Day 2 vs. 7 B trial recalling';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);

mean_txt = mean_2_7.PV.recall.B(1);
sem_txt = sem_2_7.PV.recall.B(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

mean_txt = mean_2_7.PV.recall.B(2);
sem_txt = sem_2_7.PV.recall.B(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%number of samples
txt_input = ['nb samp ', num2str(nb_comp_2_7.PV.recall.B)];
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);


%% PV learn vs. recall d7
writeWordEnter(ActXWord,WordHandle,1);

%A 2 7 learn vs. recalling
txt_input = 'PV Day 2 vs. 7 A learn vs recall';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_7_lr.PV.A(1);
sem_txt = sem_7_lr.PV.A(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_7_lr.PV.A(2);
sem_txt = sem_7_lr.PV.A(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%number of samples
txt_input = ['nb samp ', num2str(nb_comp_7.PV.A)];
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);


%B 7 learn vs. recalling
txt_input = 'PV Day 2 vs. 7 B learn vs recall';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);

mean_txt = mean_7_lr.PV.B(1);
sem_txt = sem_7_lr.PV.B(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

mean_txt = mean_7_lr.PV.B(2);
sem_txt = sem_7_lr.PV.B(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%number of samples
txt_input = ['nb samp ', num2str(nb_comp_7.PV.B)];
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);


%% TC ts learning

%A 2 7 learning
txt_input = 'TC Day 2 vs. 7 A trial learning';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_2_7.TC.ts.learn.A(1);
sem_txt = sem_2_7.TC.ts.learn.A(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_2_7.TC.ts.learn.A(2);
sem_txt = sem_2_7.TC.ts.learn.A(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%number of samples
txt_input = ['nb samp ', num2str(nb_comp_2_7.TC.ts.learn.A)];
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);

%B 2 7 learning
txt_input = 'TC ts Day 2 vs. 7 B trial learning';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);

mean_txt = mean_2_7.TC.ts.learn.B(1);
sem_txt = sem_2_7.TC.ts.learn.B(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

mean_txt = mean_2_7.TC.ts.learn.B(2);
sem_txt = sem_2_7.TC.ts.learn.B(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%number of samples
txt_input = ['nb samp ', num2str(nb_comp_2_7.TC.ts.learn.B)];
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%% TC.ts recall
writeWordEnter(ActXWord,WordHandle,1);

%A 2 7 recalling
txt_input = 'TC ts Day 2 vs. 7 A recall recalling';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_2_7.TC.ts.recall.A(1);
sem_txt = sem_2_7.TC.ts.recall.A(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_2_7.TC.ts.recall.A(2);
sem_txt = sem_2_7.TC.ts.recall.A(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%number of samples
txt_input = ['nb samp ', num2str(nb_comp_2_7.TC.ts.recall.A)];
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);

%B 2 7 recalling
txt_input = 'TC ts Day 2 vs. 7 B trial recalling';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);

mean_txt = mean_2_7.TC.ts.recall.B(1);
sem_txt = sem_2_7.TC.ts.recall.B(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

mean_txt = mean_2_7.TC.ts.recall.B(2);
sem_txt = sem_2_7.TC.ts.recall.B(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%number of samples
txt_input = ['nb samp ', num2str(nb_comp_2_7.TC.ts.recall.B)];
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);


%% TC ts learn vs. recall d7
writeWordEnter(ActXWord,WordHandle,1);

%A 2 7 learn vs. recalling
txt_input = 'TC ts Day 2 vs. 7 A learn vs recall';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_7_lr.TC.ts.A(1);
sem_txt = sem_7_lr.TC.ts.A(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_7_lr.TC.ts.A(2);
sem_txt = sem_7_lr.TC.ts.A(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%number of samples
txt_input = ['nb samp ', num2str(nb_comp_7.TC.ts.A)];
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);


%B 7 learn vs. recalling
txt_input = 'TC ts Day 2 vs. 7 B learn vs recall';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);

mean_txt = mean_7_lr.TC.ts.B(1);
sem_txt = sem_7_lr.TC.ts.B(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

mean_txt = mean_7_lr.TC.ts.B(2);
sem_txt = sem_7_lr.TC.ts.B(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%number of samples
txt_input = ['nb samp ', num2str(nb_comp_7.TC.ts.B)];
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%% Close Word document
%CloseWord(ActXWord,WordHandle,FileSpec);

