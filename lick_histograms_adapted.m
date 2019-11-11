

%% load in licks, behavior (trial order from each day) into struct 
%list directories

%I45 LT (1,1,1,1)
dirList{1} = {'G:\Figure_1_OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_rand_d1_052218',...
            'G:\Figure_1_OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_5A5B_053018',...
            'G:\Figure_1_OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_3A3B_060518',...
            'G:\Figure_1_OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_AB_061418'};
%I46 (1,1,1,1)
dirList{2} = {'G:\Figure_1_OCGOL_learning_long_term\I46\behavior_only\I46_rand_d1_052918',...
                'G:\Figure_1_OCGOL_learning_long_term\I46\behavior_only\I46_5A5B_060118',...
                'G:\Figure_1_OCGOL_learning_long_term\I46\behavior_only\I46_3A3B_060718',...
                'G:\Figure_1_OCGOL_learning_long_term\I46\behavior_only\I46_AB_061518'};
            
%I47 RS (1,0,1,1)
dirList{3} = {'G:\Figure_1_OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_rand_d2_051518',...
                'G:\Figure_1_OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_5AB_d7_052218',...
                'G:\Figure_1_OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_3AB_d8_052418',...
                'G:\Figure_1_OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_AB_061418'};
            
%I47LP (1,1,1,1)
dirList{4} = {'G:\Figure_1_OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_rand_d2_051518',...
                'G:\Figure_1_OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_5AB_d1_051718',...
                'G:\Figure_1_OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_3AB_d8_052418',...
                'G:\Figure_1_OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_AB_061418'};


%% Load in Behavior struct with lick data

%for each animal
tic
for aa=1:4
    for ss = 1:4
        lick_data{aa}{ss} = load(fullfile(dirList{aa}{ss},'output','Behavior.mat'),'Behavior');
    end
end
toc;

%% Take all lick positions and bin 


% figure
% hold on
% ylim([0 0.5])
% b= bar(sum(lick_data{1, 1}{1, 1}.Behavior.lick.bin_licks_lap,1)./...
%     sum(sum(lick_data{1, 1}{1, 1}.Behavior.lick.bin_licks_lap,1)),1,'hist')
% b.FaceColor = 'g'


%% Normalize lap positions in lick struct to this matrix

%get max pos in cm for each lap
for aa=1:4
    for ss=1:4
        max_pos{aa}{ss} = lick_data{aa}{ss}.Behavior.position_lap(:,2);
    end
end
%get normalized position of each lick on each lap
for aa=1:4
    %for each session
    for ss=1:4
        %for each lap
        for ll=1:size(lick_data{aa}{ss}.Behavior.lick.lick_lap,2)
            lick_lap_norm{aa}{ss}{ll} = lick_data{aa}{ss}.Behavior.lick.lick_lap{ll}.position./max_pos{aa}{ss}(ll);
        end
    end
end

histogram(cell2mat(lick_lap_norm{aa}{ss}'),0:0.025:1,'Normalization','probability')


%% Split into A and B normalized lick positions

for aa=1:4
    %for each session (skip RF sessions)
    for ss=2:4
        %correct and incorrect A trials
        A_idx = find(lick_data{aa}{ss}.Behavior.performance.trialOrder == 2 | lick_data{aa}{ss}.Behavior.performance.trialOrder ==20);
        lick_lap_norm_A{aa}{ss} = lick_lap_norm{aa}{ss}(A_idx);

        %correct and incorrect B trials
        B_idx = find(lick_data{aa}{ss}.Behavior.performance.trialOrder == 3 | lick_data{aa}{ss}.Behavior.performance.trialOrder ==30);
        lick_lap_norm_B{aa}{ss} = lick_lap_norm{aa}{ss}(B_idx);
        
    end
end

%% Get mean/median normalized A/B reward position for each animal

for aa=1:4
    %for each session (skip RF sessions)
    for ss=2:4
        %reward A normalized position onset
        rew_A_mean(aa,ss) = mean(lick_data{aa}{ss}.Behavior.rewards{2}.position_norm);
        %reward B normalized position onset
        rew_B_mean(aa,ss) = mean(lick_data{aa}{ss}.Behavior.rewards{1}.position_norm);
    end
end
%mean from all animals on each session
rew_A_mean_all = mean(rew_A_mean,1);
rew_B_mean_all = mean(rew_B_mean,1);

%reward range for plotting (50 bins)
rew_A_mean_all(2)*50, (rew_A_mean_all(2)+0.051)*50;
rew_B_mean_all(2)*50, (rew_B_mean_all(2)+0.051)*50;

%% Calculate fraction of licks in reward zone and plot (Figure 1C)

%get mean of max positon for each animal, session
for aa=1:4
        mean_max_pos(:,aa) = cellfun(@mean,max_pos{aa});
end

%10 cm distance on normalized position
range_norm = (ones(4,4)*10)./mean_max_pos;
mean_range_norm = mean(mean(range_norm));

%get reward ranges for
rew_A_norm = mean(rew_A_mean_all(2:end));
rew_B_norm = mean(rew_B_mean_all(2:end));


%merge lick positions across all laps
for aa=1:4
    for ss=1:4
        if ss ~= 1
            lick_lap_norm_A_all_laps{aa}{ss} = cell2mat(lick_lap_norm_A{aa}{ss}');
            lick_lap_norm_B_all_laps{aa}{ss} = cell2mat(lick_lap_norm_B{aa}{ss}');
        elseif ss==1
            lick_lap_norm_all_laps{aa}{1} = cell2mat(lick_lap_norm{aa}{1}');
        end
    end
end

%find number of licks in reward zones A for A laps; B for B laps
for aa=1:4
    for ss=1:4
        if ss~=1
            lick_count_Azone(ss,aa) = size(find(lick_lap_norm_A_all_laps{aa}{ss} >= rew_A_norm & lick_lap_norm_A_all_laps{aa}{ss} <= (rew_A_norm +mean_range_norm)),1);
            lick_count_Bzone(ss,aa) = size(find(lick_lap_norm_B_all_laps{aa}{ss} >= rew_B_norm & lick_lap_norm_B_all_laps{aa}{ss} <= (rew_B_norm +mean_range_norm)),1);
        elseif ss==1
            lick_count_Azone(1,aa) = size(find(lick_lap_norm_all_laps{aa}{ss} >= rew_A_norm & lick_lap_norm_all_laps{aa}{ss} <= (rew_A_norm +mean_range_norm)),1);
            lick_count_Bzone(1,aa) = size(find(lick_lap_norm_all_laps{aa}{ss} >= rew_B_norm & lick_lap_norm_all_laps{aa}{ss} <= (rew_B_norm +mean_range_norm)),1);   
        end
    end
end

%get total A lap and B lap licks
for aa=1:4
    for ss=1:4
        if ss~=1
            lick_count_A(ss,aa) = size(lick_lap_norm_A_all_laps{aa}{ss},1);
            lick_count_B(ss,aa) = size(lick_lap_norm_B_all_laps{aa}{ss},1);
        elseif ss==1
            lick_count_A(1,aa) = size(lick_lap_norm_all_laps{aa}{ss},1);
            lick_count_B(1,aa) = size(lick_lap_norm_all_laps{aa}{ss},1);   
        end
    end
end

%fraction of licks in reward zone on A and B laps
frac_zone_A = lick_count_Azone./lick_count_A;
frac_zone_B = lick_count_Bzone./lick_count_B;

%get mean fraction of lick in A and B laps
mean_zone_A = mean(frac_zone_A,2);
mean_zone_B = mean(frac_zone_B,2);

%sem of fraction of licks by animal
sem_zone_A = std(frac_zone_A,0,2)./sqrt(4);
sem_zone_B = std(frac_zone_B,0,2)./sqrt(4);

%line plot
figure
hold on
xlim([0.5 4.5])
xticks(1:4)
xticklabels({'Random\newline foraging' ,'5A5B','3A3B','Random\newline AB'})
ylim([0 1.2])
yticks(0:0.2:1)
ylabel('Fraction of licks in reward zone')
ae = errorbar(mean_zone_A,sem_zone_A,'Color',[65,105,225]./255,'LineWidth',2);
be = errorbar(mean_zone_B,sem_zone_B,'Color',[220,20,60]./255,'LineWidth',2);

legend([ae be],{'A','B'},'Location','northwest')

set(gca,'FontSize',16)
set(gca,'LineWidth',2)

%% Statistics for reward zone lick - perform 1-way anova and t-test with Bonferroni correction for mc

%each column is diff session - A
[p_zoneA,tbl_zoneA,stats_zoneA] = anova1(frac_zone_A',{'RF','5A5B','3A3B','Rand AB'});
%each coumns is diff session - B 
[p_zoneB,tbl_zoneB,stats_zoneB] = anova1(frac_zone_B',{'RF','5A5B','3A3B','Rand AB'});

%get columns so that they are sessions
frac_zone_A_ses = frac_zone_A'; 
frac_zone_B_ses = frac_zone_B';

%do 3 paired t-tests for R- vs randAB; 3A3B vs. randAB, 5A5B vs. rand AB
%RF vs. rand AB; sig after Holm-Sidak correction
[~,p_Azone(1)] =ttest(frac_zone_A_ses(:,1),frac_zone_A_ses(:,4));
%5A5B vs. rand AB
[~,p_Azone(2)] =ttest(frac_zone_A_ses(:,2),frac_zone_A_ses(:,4));
%3A3V vs. rand AB
[~,p_Azone(3)] =ttest(frac_zone_A_ses(:,3),frac_zone_A_ses(:,4));

%RF vs. rand AB; sig after Holm-Sidak correction
[~,p_Bzone(1)] =ttest(frac_zone_B_ses(:,1),frac_zone_B_ses(:,4));
%5A5B vs. rand AB
[~,p_Bzone(2)] =ttest(frac_zone_B_ses(:,2),frac_zone_B_ses(:,4));
%3A3V vs. rand AB
[~,p_Bzone(3)] =ttest(frac_zone_B_ses(:,3),frac_zone_B_ses(:,4));

%post hoc Holm-Sidak multi comp correction for 3 tests (Prism)

%Dunn-Sidak correction (below),  Holm-Sidak (Prism analysis)
%p_sidak = 1-(1-0.05)^(1/3)
%p_sidak = 1-(1-0.01)^(1/3)
%p_sidak = 1-(1-0.001)^(1/3)

%% Calculate percent of trials correct and construct line plot (Figure 1D)

for aa=1:4
    for ss=2:4
        trialData = lick_data{aa}{ss}.Behavior.performance.trialOrder;
        %frac A corr, frac B corr, all corr
        corr_mat(ss,1,aa) = size(find(trialData == 2),1)./size(find(trialData == 2 | trialData == 20),1);
        corr_mat(ss,2,aa) = size(find(trialData == 3),1)./size(find(trialData == 3 | trialData == 30),1);
        corr_mat(ss,3,aa) = size(find(trialData == 2 | trialData == 3),1)./size(trialData,1);
    end
end

%get mean and sem for A and B
mean_corr = mean(corr_mat,3);
sem_corr = std(corr_mat,0,3)./sqrt(4);

%line plot
figure
hold on
xlim([0.5 4.5])
xticks(1:4)
xticklabels({'Random\newline foraging' ,'5A5B','3A3B','Random\newline AB'})
xtickangle(45)
ylim([0 1.2])
yticks([0:0.2:1])
ylabel('Fraction correct')
%mean A corr
ae = errorbar(mean_corr(:,1),sem_corr(:,1),'Color',[65,105,225]./255,'LineWidth',2);
%mean B corr
be = errorbar(mean_corr(:,2),sem_corr(:,2),'Color',[220,20,60]./255,'LineWidth',2);
%mean all corr
abe = errorbar(mean_corr(:,3),sem_corr(:,3),'Color',[139, 0, 139]./255,'LineWidth',2);

legend([ae be abe],{'A','B','All'},'Location','northwest')

set(gca,'FontSize',16)
set(gca,'LineWidth',2)

%% Generate histogram for figure animal

%which animal to generate plot from
aa=2;
%50 bins
edges = (0:0.02:1);

%magenta = [0.85, 0.42, 0.68];
green = [50,205,50]./255;
blue = [65,105,225]./255;
red =  [220,20,60]./255;

figure('Position', [100 100 600 800]);

subplot(4,2,[1 2])
hold on
title('Random')
plotHisto(cell2mat(lick_lap_norm{aa}{1}'),edges,green);
hold off

subplot(4,2,3)
hold on
title('A')
plotHisto(cell2mat(lick_lap_norm_A{aa}{2}'),edges,blue);
hold off

subplot(4,2,4)
hold on
title('B')
plotHisto(cell2mat(lick_lap_norm_B{aa}{2}'),edges,red);
hold off

subplot(4,2,5)
hold on
title('A')
plotHisto(cell2mat(lick_lap_norm_A{aa}{3}'),edges,blue);
hold off

subplot(4,2,6)
hold on
title('B')
plotHisto(cell2mat(lick_lap_norm_B{aa}{3}'),edges,red);
hold off

subplot(4,2,7)
hold on
title('A')
plotHisto(cell2mat(lick_lap_norm_A{aa}{4}'),edges,blue);
hold off

subplot(4,2,8)
hold on
title('B')
plotHisto(cell2mat(lick_lap_norm_B{aa}{4}'),edges,red);
hold off

