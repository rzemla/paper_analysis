

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



%% Get percent correct plot
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

