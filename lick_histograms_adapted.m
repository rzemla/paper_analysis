

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


figure
hold on
ylim([0 0.5])
b= bar(sum(lick_data{1, 1}{1, 1}.Behavior.lick.bin_licks_lap,1)./...
    sum(sum(lick_data{1, 1}{1, 1}.Behavior.lick.bin_licks_lap,1)),1,'hist')
b.FaceColor = 'g'


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


%% iterate the sequence below

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
    
        %if RF session
    elseif ss == 1
        stages{ss}.lickCombined = [];
        for ii=1:size(stages{ss}.licks,2)
            stages{ss}.lickCombined = [stages{ss}.lickCombined; stages{ss}.licks{ii}(:,2)];
        end
        
    end

end


%% Generate histograms

%which animal to generate plot from
aa=2;
%40 bins
edges = (0:0.025:1);

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

