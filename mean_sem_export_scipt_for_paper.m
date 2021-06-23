%% import prism fraction of licks data on A and B trials

dirpath = 'C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data\fig1_csv_data_prism_export';

%fraction of licks in reward zone
frac_licks_A = readmatrix(fullfile(dirpath,'A fraction of licks in reward zone.csv'),'Range','B2:E5');
frac_licks_B = readmatrix(fullfile(dirpath,'B fraction of licks in reward zone.csv'),'Range','B2:E5');

%fraction correct trials
frac_corr_A = readmatrix(fullfile(dirpath,'Fraction correct A trials.csv'),'Range','B2:E5');
frac_corr_B = readmatrix(fullfile(dirpath,'Fraction correct B trials.csv'),'Range','B2:E5');

%% Mean and sem for licking fraction of correct trials across training sessions

nb_animals = size(frac_licks_A,2);
%licking mean and sem
mean_lick.A = mean(frac_licks_A,1);
mean_lick.B = mean(frac_licks_B,1);

sem_lick.A = std(frac_licks_A,0,1)./sqrt(nb_animals);
sem_lick.B = std(frac_licks_B,0,1)./sqrt(nb_animals);

%fraction correct mean and sem
mean_corr.A = mean(frac_corr_A,1);
mean_corr.B = mean(frac_corr_B,1);

sem_corr.A = std(frac_corr_A,0,1)./sqrt(nb_animals);
sem_corr.B = std(frac_corr_B,0,1)./sqrt(nb_animals);

