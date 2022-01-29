%% Load in all datasets

%I47LP
dirlist{1} = {'E:\I47_LP_rand_d2_051518_lick\output',...
    'E:\I47_LP_5AB_d1_051718_lick\output',...
    'E:\I47_LP_3AB_d8_052418_lick\output',...
    'E:\I47_LP_AB_061418_lick\output'};

%I47 RS
dirlist{2} = {'E:\I47_RS_rand_d2_051518_lick\output',...
    'E:\I47_RS_5AB_d1_051618_lick\output',...
    'E:\I47_RS_3AB_d8_052418_lick\output',...
    'E:\I47_RS_AB_061418_lick\output'};

%I46
dirlist{3} = {'E:\I46_rand_d1_052918_lick\output','E:\I46_5A5B_060118_lick\output',...
    'E:\I46_3A3B_060718_lick\output','E:\I46_AB_061518_lick\output'};
%I45 LT
dirlist{4} = {'E:\I45_RT_rand_d1_052218_lick\output','E:\I45_RT_5A5B_053018_lick\output',...
    'E:\I45_RT_3A3B_060518_lick\output','E:\I45_RT_AB_061418_lick\output'};

%load in Behaviors and licks variables
%for each animal
for jj = 1:size(dirlist,2)
    %for each session
    for ii = 1:size(dirlist{jj},2)
        stages{jj}{ii} = load([dirlist{jj}{ii} '\matlab.mat'],'licks','Behavior');
    end
end

%% Combine lick position into one vector

%for each animal
for jj = 1:size(dirlist,2)
    %for each session
    for ss = 1:size(dirlist,2)
        
        if ss ~= 1
            
            trialOrder = stages{jj}{ss}.Behavior{1}.performance.trialOrder;
            
            %take only correct A trials
            stages{jj}{ss}.lickCombined{1} = [];
            stages{jj}{ss}.lickCombined{2} = [];
            
            %for each lap
            for ii=1:size(stages{jj}{ss}.licks,2)
                if trialOrder(ii) == 2 || trialOrder(ii) == 20
                    %combine all lick positions into one vector
                    stages{jj}{ss}.lickCombined{1} = [stages{jj}{ss}.lickCombined{1}; stages{jj}{ss}.licks{ii}(:,2)];
                elseif trialOrder(ii) == 3 || trialOrder(ii) == 30
                    stages{jj}{ss}.lickCombined{2} = [stages{jj}{ss}.lickCombined{2}; stages{jj}{ss}.licks{ii}(:,2)];
                end
                
            end
            
        elseif ss == 1
            stages{jj}{ss}.lickCombined = [];
            for ii=1:size(stages{jj}{ss}.licks,2)
                stages{jj}{ss}.lickCombined = [stages{jj}{ss}.lickCombined; stages{jj}{ss}.licks{ii}(:,2)];
            end
            
        end
    end
end

%% Calculate fraction of licks in reward zones

%figure these out for correct reward range
% %reward range for trials A (by position)
 rewardRange{1} = [135 149];
 rewardRange{2} = [55 69];

 
%reward range for trials B (by position)
% rewardRange{1} = [135 154];
% rewardRange{2} = [55 74];

%for each animal
for jj = 1:size(dirlist,2)
    %for each session
    for ss = 1:size(dirlist,2)
        if ss ~= 1
            % A reward zone
            fracRew(ss,1,jj) = numel(find(stages{jj}{ss}.lickCombined{1} >= rewardRange{1}(1) & ...
                stages{jj}{ss}.lickCombined{1} <= rewardRange{1}(2)))/numel(stages{jj}{ss}.lickCombined{1});
            
            % B reward zone
            fracRew(ss,2,jj) = numel(find(stages{jj}{ss}.lickCombined{2} >= rewardRange{2}(1) & ...
                stages{jj}{ss}.lickCombined{2} <= rewardRange{2}(2)))/numel(stages{jj}{ss}.lickCombined{2});
        
        elseif ss == 1
            % A reward zone
            fracRew(ss,1,jj) = numel(find(stages{jj}{ss}.lickCombined >= rewardRange{1}(1) & ...
                stages{jj}{ss}.lickCombined <= rewardRange{1}(2)))/numel(stages{jj}{ss}.lickCombined);
            % B reward zone
            fracRew(ss,2,jj) = numel(find(stages{jj}{ss}.lickCombined >= rewardRange{2}(1) & ...
                stages{jj}{ss}.lickCombined <= rewardRange{2}(2)))/numel(stages{jj}{ss}.lickCombined);
        end
    end
end

%% Plot the line graph

%colors
magenta = [0.85, 0.42, 0.68];
blue = [0.25,0.41,0.88];
red =  [0.85,0.07,0.23];

%Mean and SEM for each data point across animals
%mean
fracMean = mean(fracRew,3);

%sem
fracSem = std(fracRew,0,3)./sqrt(size(fracRew,3));

%Wilcoxon signed rank test (paired Mann Whitney U test) on 5A5B tests vs rand AB
%too few samples
% [p,h] = signrank(squeeze(fracRew(2,1,:)), squeeze(fracRew(4,1,:)))
% [p,h] = signrank(squeeze(fracRew(2,2,:)), squeeze(fracRew(4,2,:)))

%Kruskall Wallis H test - did RM ANOVA in PRISM

%A trials
%extended range = 0.038
%narrow range  = 0.2857
% [p,tbl, stats] = kruskalwallis(squeeze(fracRew(:,1,:))')
% mc = multcompare(stats)
% 
% %B trials
% %extended range = 0.026
% %narrow range  = 0.285
% [p,tbl, stats] = anova1(squeeze(fracRew(:,1,:))')
% mc = multcompare(stats)
% 
% 
% %B trials
% [p,tbl, stats]  = kruskalwallis(squeeze(fracRew(:,2,:))')
% mc = multcompare(stats)

f = figure('Position', [100 100 900 600]);
hold on
ax = gca;
ax.FontSize = 10.5;
%A fraction
%plot([1:4], fracMean(:,1),'o','MarkerEdgeColor', blue, 'MarkerFaceColor', blue)
eA = errorbar([1:4],fracMean(:,1),fracSem(:,1),'-o','MarkerFaceColor', blue,'MarkerEdgeColor', blue, 'Color', blue);

%B fraction
%plot([1:4], fracMean(:,2),'s','MarkerEdgeColor', red, 'MarkerFaceColor', red)
eB = errorbar([1:4],fracMean(:,2),fracSem(:,2),'-s','MarkerFaceColor', red,'MarkerEdgeColor', red, 'Color', red);

xlabel('Learning stage')
xlim([0 5])
xticks([1:4]);
%xticklabels({'Random \newline reward','5A5B','3A3B','Random AB'})
xticklabels({'Stage 1','Stage 2','Stage 3','Final'})
xtickangle(45)

ylabel('Fraction of licks in reward zone')
ylim([0 1.1])
yticks([0.1 0.4 0.7 1])
%number of animals
text([4.5],[0.1],'n = 4','FontSize',12); 

%add significance comparison horizontal bars
%hA = sigstar({[2,4]},[0.05]);
%set(hA,'Color',blue)

%hB = sigstar({[2,4]},[0.05]);
%set(hB,'Color',red);

hold off

legend([eA,eB], {'A', 'B'},'Location','northwest', 'FontSize', 11)
legend('boxoff')

%% Calculate fraction of correct trials on each day

% for each OCGOL days get fraction correct A and B trials

for jj=1:4
    for ss=2:4
        fracCorr(ss-1,1,jj) = numel(find(stages{jj}{ss}.Behavior{1}.performance.trialOrder ==2))/...
            numel(find(stages{jj}{ss}.Behavior{1}.performance.trialOrder ==2 | stages{jj}{ss}.Behavior{1}.performance.trialOrder ==20));
        fracCorr(ss-1,2,jj) = numel(find(stages{jj}{ss}.Behavior{1}.performance.trialOrder ==3))/...
            numel(find(stages{jj}{ss}.Behavior{1}.performance.trialOrder ==3 | stages{jj}{ss}.Behavior{1}.performance.trialOrder ==30));
    end
end

fracCorrMean = mean(fracCorr,3);
fracCorrSem = std(fracCorr,0,3)/sqrt(4);


%% Plot fraction of correct trials

f = figure('Position', [100 100 900 600]);
hold on
ax = gca;
ax.FontSize = 10.5;
%A fraction
%plot([1:4], fracMean(:,1),'o','MarkerEdgeColor', blue, 'MarkerFaceColor', blue)
eA = errorbar([1:3],fracCorrMean(:,1),fracCorrSem(:,1),'-o','MarkerFaceColor', blue,'MarkerEdgeColor', blue, 'Color', blue);

%B fraction
%plot([1:4], fracMean(:,2),'s','MarkerEdgeColor', red, 'MarkerFaceColor', red)
eB = errorbar([1:3],fracCorrMean(:,2),fracCorrSem(:,2),'-s','MarkerFaceColor', red,'MarkerEdgeColor', red, 'Color', red);

xlabel('Learning stage')
xlim([0 4])
xticks([1:3]);
xticklabels({'Stage 2','Stage 3','Final'})
xtickangle(45)

ylabel('Fraction of correct trials')
ylim([0 1.1])
yticks([0.1 0.4 0.7 1])

%number of animals
text([3.5],[0.1],'n = 4','FontSize',12); 

hold off

legend([eA,eB], {'A', 'B'},'Location','northwest', 'FontSize', 11)
legend('boxoff')

%ANOVA stats
%5A5B vs rand AB comparison
squeeze(fracCorr(:,1,:))'



