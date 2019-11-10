

%% load in licks, behavior (trial order from each day) into struct 
%list directories

%I45 LT (1,1,1,1)
dirlist = {'E:\I45_RT_rand_d1_052218_lick\output',...
            'E:\I45_RT_5A5B_053018_lick\output',...
            'E:\I45_RT_3A3B_060518_lick\output',...
            'E:\I45_RT_AB_061418_lick\output'};

%I46 (1,1,1,1)
dirlist = {'E:\I46_rand_d1_052918_lick\output',...
            'E:\I46_5A5B_060118_lick\output',...
            'E:\I46_3A3B_060718_lick\output',...
            'E:\I46_AB_061518_lick\output'};

%I47 RS (1,0,1,1)
dirlist = {'E:\I47_RS_rand_d2_051518_lick\output',...
            'E:\I47_RS_5AB_d1_051618_lick\output',...
            'E:\I47_RS_3AB_d8_052418_lick\output',...
            'E:\I47_RS_AB_061418_lick\output'};

%I47LP (1,1,1,1)
dirlist = {'E:\I47_LP_rand_d2_051518_lick\output',...
            'E:\I47_LP_5AB_d1_051718_lick\output',...
            'E:\I47_LP_3AB_d8_052418_lick\output',...
            'E:\I47_LP_AB_061418_lick\output'};

%load in Behaviors and licks variables
for ii = 1:size(dirlist,2)
   stages{ii} = load([dirlist{ii} '\matlab.mat'],'licks','Behavior');    
end

%% iterate the sequence below

edges = (0:5:200);

for ss = 1:size(dirlist,2)

    if ss ~= 1
        
        trialOrder = stages{ss}.Behavior{1}.performance.trialOrder;
        
        %take only correct A trials
        stages{ss}.lickCombined{1} = [];
        stages{ss}.lickCombined{2} = [];
        
        for ii=1:size(stages{ss}.licks,2)
            if trialOrder(ii) == 2 || trialOrder(ii) == 20
                %combine all lick positions into one vector
                stages{ss}.lickCombined{1} = [stages{ss}.lickCombined{1}; stages{ss}.licks{ii}(:,2)];
            elseif trialOrder(ii) == 3 || trialOrder(ii) == 30
                stages{ss}.lickCombined{2} = [stages{ss}.lickCombined{2}; stages{ss}.licks{ii}(:,2)];
            end
            
        end
        
    elseif ss == 1
        stages{ss}.lickCombined = [];
        for ii=1:size(stages{ss}.licks,2)
            stages{ss}.lickCombined = [stages{ss}.lickCombined; stages{ss}.licks{ii}(:,2)];
        end
        
    end

end


%% Generate histograms
magenta = [0.85, 0.42, 0.68];
blue = [0.25,0.41,0.88];
red =  [0.85,0.07,0.23];

figure('Position', [100 100 600 800]);

subplot(4,2,[1 2])
hold on
title('Random')
plotHisto(stages{1}.lickCombined,edges,magenta);
hold off

subplot(4,2,3)
hold on
title('A')
plotHisto(stages{2}.lickCombined{1},edges,blue);
hold off

subplot(4,2,4)
hold on
title('B')
plotHisto(stages{2}.lickCombined{2},edges,red);
hold off

subplot(4,2,5)
hold on
title('A')
plotHisto(stages{3}.lickCombined{1},edges,blue);
hold off

subplot(4,2,6)
hold on
title('B')
plotHisto(stages{3}.lickCombined{2},edges,red);
hold off

subplot(4,2,7)
hold on
title('A')
plotHisto(stages{4}.lickCombined{1},edges,blue);
hold off

subplot(4,2,8)
hold on
title('B')
plotHisto(stages{4}.lickCombined{2},edges,red);
hold off

%GOAL PLOTS
% fraction lick histograms for each behavior stage for each animal on A and B trials


%fraction of correct A and B trials at each timepoint for each animal

