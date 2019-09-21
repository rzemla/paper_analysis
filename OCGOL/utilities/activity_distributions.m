function [outputArg1,outputArg2] = activity_distributions(session_vars,task_selective_ROIs,options)
%plot AUC rate and event rate for animal across imaging sessions

%% Define variables
selectTrial = options.selectTrial;
sessionSelect = options.sessionSelect;

%% Plot scatter of task-selective ROIs vs AUC/min rate on opposing trials
%as function of AUC/min
figure
for ss=1:6
    subplot(1,6,ss)
    hold on
    xlim([0 3])
    ylim([0 25])
    %AUC/min in B for A-selective ROIs
    scatter(ones(1,length(task_selective_ROIs{ss}.A.idx)),session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.AUC_min(task_selective_ROIs{ss}.A.idx),...
        12,'filled')
    scatter(1,median(session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.AUC_min(task_selective_ROIs{ss}.A.idx)),20,'filled','MarkerFaceColor','k')
    %AUC/min in A for B-selective ROIs
    scatter(2*ones(1,length(task_selective_ROIs{ss}.B.idx)),session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.AUC_min(task_selective_ROIs{ss}.B.idx),...
        12,'filled')
    scatter(2,median(session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.AUC_min(task_selective_ROIs{ss}.B.idx)),20,'filled','MarkerFaceColor','k')
end

%% Plot scatter of task-selective ROIs vs event frequency on opposing trials
figure
for ss=1:6
    subplot(1,6,ss)
    hold on
    xlim([0 3])
    ylim([0 15])
    %AUC/min in B for A-selective ROIs
    scatter(ones(1,length(task_selective_ROIs{ss}.A.idx)),session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.event_freq(task_selective_ROIs{ss}.A.idx),...
        12,'filled')
    scatter(1,median(session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.event_freq(task_selective_ROIs{ss}.A.idx)),20,'filled','MarkerFaceColor','k')
    %AUC/min in A for B-selective ROIs
    scatter(2*ones(1,length(task_selective_ROIs{ss}.B.idx)),session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.event_freq(task_selective_ROIs{ss}.B.idx),...
        12,'filled')
    scatter(2,median(session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.event_freq(task_selective_ROIs{ss}.B.idx)),20,'filled','MarkerFaceColor','k')
end
%% Plot AUC/min of A vs on scatter for A-selective and B-selective neurons

%scatter
figure;
for ss=1:6
subplot(1,6,ss)
hold on
axis square
title('AUC/min')
xlim([0 40])
ylim([0 40])
xlabel('A')
ylabel('B')
scatter(session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.AUC_min(task_selective_ROIs{ss}.A.idx),session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.AUC_min(task_selective_ROIs{ss}.A.idx),14,'filled','b')
scatter(session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.AUC_min(task_selective_ROIs{ss}.B.idx),session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.AUC_min(task_selective_ROIs{ss}.B.idx),14,'filled','r')
plot([0,40],[0 40],'k--','LineWidth',1.5)
end

%% Plot AUC/min of A vs on scatter for all neurons - run

%scatter
figure;
for ss=1:6
subplot(1,6,ss)
hold on
axis square
title('AUC/min')
xlim([0 40])
ylim([0 40])
xlabel('A')
ylabel('B')
scatter(session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.AUC_min,session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.AUC_min,7,'filled','k')
plot([0,40],[0 40],'k--','LineWidth',1.5)
end

%% Plot AUC/min of A vs on scatter for all neurons - NO run
figure;
for ss=1:6
subplot(1,6,ss)
hold on
axis square
title('AUC/min')
xlim([0 40])
ylim([0 40])
xlabel('A')
ylabel('B')
scatter(session_vars{ss}.Events_split{selectTrial(1)}.NoRun.properties.AUC_min,session_vars{ss}.Events_split{selectTrial(2)}.NoRun.properties.AUC_min,7,'filled','k')
plot([0,40],[0 40],'k--','LineWidth',1.5)
end

%% Plot AUC/min of A vs on scatter for all neurons - Run vs no run
figure;
for ss=sessionSelect
subplot(1,6,ss)
hold on
axis square
title('Run vs. No run AUC/min')
xlim([0 40])
ylim([0 40])
xlabel('No Run')
ylabel('Run')
scatter(session_vars{ss}.Events_split{selectTrial(1)}.NoRun.properties.AUC_min,session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.AUC_min,7,'filled','b')
scatter(session_vars{ss}.Events_split{selectTrial(2)}.NoRun.properties.AUC_min,session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.AUC_min,7,'filled','r')
plot([0,40],[0 40],'k--','LineWidth',1.5)
end


%% 
%% Plot Event freq of A vs on scatter for all neurons - Run vs no run
figure;
for ss=sessionSelect
subplot(1,6,ss)
hold on
axis square
title('Run vs. No run Event freq.')
xlim([0 40])
ylim([0 40])
xlabel('No run')
ylabel('Run')
scatter(session_vars{ss}.Events_split{selectTrial(1)}.NoRun.properties.event_freq,session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.event_freq,7,'filled','b')
scatter(session_vars{ss}.Events_split{selectTrial(2)}.NoRun.properties.event_freq,session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.event_freq,7,'filled','r')
plot([0,40],[0 40],'k--','LineWidth',1.5)
end

end

