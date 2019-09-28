function [outputArg1,outputArg2] = SCE_rate_plots(session_vars_learn,session_vars_recall,SCE_learning,SCE_recall, perf_learning,perf_recall)

%% Get minutes total spent in no-run sessions

%get minutes spent in each no run epoch for each session for each animal
%(total)
%learn
for aa= 1:size(session_vars_learn,2)
    for ss=1:6
        frame_norun_learn(aa,ss) = length(find((session_vars_learn{aa}{ss}.Behavior_split{3}.run_ones  ==0) ==1));
    end
end
%recall
for aa= 1:size(session_vars_recall,2)
    for ss=1:7
        frame_norun_recall(aa,ss) = length(find((session_vars_recall{aa}{ss}.Behavior_split{3}.run_ones  ==0) ==1));
    end
end

%learn - get no-run epoch frame count for each session
for aa= 1:size(session_vars_learn,2)
    for ss=1:6
        norun_fr_A_learn(aa,ss) =  length(find((session_vars_learn{aa}{ss}.Behavior_split{4}.run_ones ==0) == 1));
        norun_fr_B_learn(aa,ss) =  length(find((session_vars_learn{aa}{ss}.Behavior_split{5}.run_ones ==0) == 1));
    end
end

%recall - get no-run epoch frame count for each session
for aa= 1:size(session_vars_recall,2)
    for ss=1:7
        norun_fr_A_recall(aa,ss) =  length(find((session_vars_recall{aa}{ss}.Behavior_split{1}.run_ones ==0) == 1));
        norun_fr_B_recall(aa,ss) =  length(find((session_vars_recall{aa}{ss}.Behavior_split{2}.run_ones ==0) == 1));
    end
end

%frames period in second (first animal, first session, recall)
dt = session_vars_recall{1}{1}.Imaging_split{1, 3}.dt;

%convert no run total frames to minutes (all recording session)
norun_learn_time = frame_norun_learn.*(dt/60);
norun_recall_time = frame_norun_recall.*(dt/60);

%convert no run frames to min (A vs. B(
norun_time_A_learn = norun_fr_A_learn.*(dt/60);
norun_time_B_learn = norun_fr_B_learn.*(dt/60);

norun_time_A_recall = norun_fr_A_recall.*(dt/60);
norun_time_B_recall = norun_fr_B_recall.*(dt/60);


%% Get SCE count and # of id'd neurons for each animal on each session
for aa= 1:size(session_vars_learn,2)
    for ss=1:6
        %all
        SCE_learn_count(aa,ss) = SCE_learning{aa}.SCE{ss}.nbSCE;
        %A
        SCE_learn_count_A(aa,ss) = size(SCE_learning{aa}.SCE{ss}.sce_activity.A,2);
        %B
        SCE_learn_count_B(aa,ss) = size(SCE_learning{aa}.SCE{ss}.sce_activity.B,2);
        %get total # of neurons in each session
        id_neurons_learn(aa,ss) =  size(session_vars_learn{aa}{ss}.Imaging_split{3}.trace_restricted,2);
    end
end

for aa= 1:size(session_vars_recall,2)
    for ss=1:7
        %all
        SCE_recall_count(aa,ss) = SCE_recall{aa}.SCE{ss}.nbSCE;
        %A
        SCE_recall_count_A(aa,ss) = size(SCE_recall{aa}.SCE{ss}.sce_activity.A,2);
        %B
        SCE_recall_count_B(aa,ss) = size(SCE_recall{aa}.SCE{ss}.sce_activity.B,2);
        %get total # of neurons in each session
        id_neurons_recall(aa,ss) =  size(session_vars_recall{aa}{ss}.Imaging_split{3}.trace_restricted,2);
    end
end

%% SCE/min
%all sessions
SCE_min_learn = SCE_learn_count./norun_learn_time;
SCE_min_recall= SCE_recall_count./norun_recall_time;

%A vs. B learn
SCE_min_learn_A = SCE_learn_count_A./norun_time_A_learn;
SCE_min_learn_B = SCE_learn_count_B./norun_time_B_learn;

%A vs. B recall
SCE_min_recall_A = SCE_recall_count_A./norun_time_A_recall;
SCE_min_recall_B = SCE_recall_count_B./norun_time_B_recall;

%% Normalized by neuron count
SCE_min_neuron_learn = SCE_min_learn./id_neurons_learn;
SCE_min_neuron_recall = SCE_min_recall./id_neurons_recall;

%% Do stats
ranksum(SCE_min_learn(:),SCE_min_recall(:))

%return x and for cdf plotting
[f1,x1] = ecdf(SCE_min_learn(:));
[f2,x2] = ecdf(SCE_min_recall(:));

figure
hold on
xlabel('SCEs/min')
ylabel('Cumulative probability')
s1 = stairs(x1,f1,'--','LineWidth',2,'Color',[139, 0, 139]/255);
s2 = stairs(x2,f2,'-','LineWidth',2,'Color',[139, 0, 139]/255);
set(gca,'FontSize',16)
set(gca,'LineWidth', 1.5)

legend([s1 s2],{'Learning','Recall'},'Location','southeast')

[h,p] = kstest2(SCE_min_learn(:),SCE_min_recall(:))


%% Scatter of learning diff as a function of SCE/min diff
for aa= 1:size(session_vars_learn,2)
    %all
    learn_perf_diff(aa,:) = diff(perf_learning{aa}.ses_perf(1,:));
    learn_SCE_min_diff(aa,:) = diff(SCE_min_learn(aa,:));
    %A
    learn_perf_diff_A(aa,:) = diff(perf_learning{aa}.ses_perf(2,:));
    learn_SCE_min_diff_A(aa,:) = diff(SCE_min_learn_A(aa,:));
    %B
    learn_perf_diff_B(aa,:) = diff(perf_learning{aa}.ses_perf(3,:));
    learn_SCE_min_diff_B(aa,:) = diff(SCE_min_learn_B(aa,:));
end

for aa= 1:size(session_vars_recall,2)
    %all
    recall_perf_diff(aa,:) = diff(perf_recall{aa}.ses_perf(1,:));
    recall_SCE_min_diff(aa,:) = diff(SCE_min_recall(aa,:));
    %A
    recall_perf_diff_A(aa,:) = diff(perf_recall{aa}.ses_perf(2,:));
    recall_SCE_min_diff_A(aa,:) = diff(SCE_min_recall_A(aa,:));
    %B
    recall_perf_diff_B(aa,:) = diff(perf_recall{aa}.ses_perf(3,:));
    recall_SCE_min_diff_B(aa,:) = diff(SCE_min_recall_B(aa,:));
end

%% Fit lines inti learning scatters - not clear b/c diff rate of learning for A vs. B
if 0
%fit line
x = learn_perf_diff(:);                             % Create Data
y = learn_SCE_min_diff(:);                          % Create Data
mdl = fitlm(x, y);                                  % Fit Data
B = mdl.Coefficients.Estimate;                      % Coefficients
CI = coefCI(mdl);                                   % Coefficient Confidence Intervals
[Ypred,YCI] = predict(mdl, x);                      % Fitted Regression Line & Confidence Intervals

Xnew = linspace(min(x), max(x), 1000);
[~,YCI] = predict(mdl, Xnew'); 

%should yield same results as corrcoef below
%[B,BINT,R,RINT,STATS] = regress( y,[x,ones(size(y',2),1)] )

%pearson R correlation
[r,p] = corrcoef(x,  y)



%% Scatterplot - change in perforance vs. change in SCE rate
figure
%subplot(1,2,1)
hold on
title('Learning')
axis square
xlabel((['\Delta','Performance [frac. correct]']))
ylabel((['\Delta','SCE rate [SCE/min]']))
xlim([-0.2 , 0.7])
ylim([-30 30])
scatter(learn_perf_diff(:), learn_SCE_min_diff(:),21,'filled','MarkerFaceColor',[139, 0, 139]/255)
%0 rate line
plot([-0.2 0.7], [0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.3);
%0 performance line
plot([0 0], [-30 30],'--','Color',[0.5 0.5 0.5],'LineWidth',0.3);

%plot confidence interval
plot(Xnew, YCI, '--k','LineWidth',1.5)
%plot predicted line
plot(x, Ypred,'-r','LineWidth',1.5)

set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)

% subplot(1,2,2)
% hold on
% title('Recall')
% axis square
% xlabel((['\Delta','Performance [frac. correct]']))
% ylabel((['\Delta','SCE rate [SCE/min]']))
% xlim([-0.4 , 0.7])
% ylim([-30 30])
% scatter(recall_perf_diff(:), recall_SCE_min_diff(:))
end

%% Fit lines inti learning scatters

%fit line
x = [learn_perf_diff_A(:);learn_perf_diff_B(:)] ;                             
y = [learn_SCE_min_diff_A(:); learn_SCE_min_diff_B(:)];

mdl = fitlm(x, y);                                  % Fit Data
B = mdl.Coefficients.Estimate;                      % Coefficients
CI = coefCI(mdl);                                   % Coefficient Confidence Intervals
[Ypred,YCI] = predict(mdl, x);                      % Fitted Regression Line & Confidence Intervals

%only within space of data
%Xnew = linspace(min(x), max(x), 1000);
%to the edges of plot 
Xnew = linspace(-0.5, 1, 1000);
[Ypred,YCI] = predict(mdl, Xnew'); 

%should yield same results as corrcoef below
%[B,BINT,R,RINT,STATS] = regress( y,[x,ones(size(y',2),1)] )

%pearson R correlation
[r,p] = corrcoef(x,  y)


%% Scatterplot - split A and B
figure
%subplot(2,1,1)
hold on
title('Learning ')
axis square
xlabel((['\Delta','Performance [frac. correct]']))
ylabel((['\Delta','SCE rate [SCE/min]']))
xlim([-0.5 , 1.0])
ylim([-20 35])
s1 = scatter(learn_perf_diff_A(:),learn_SCE_min_diff_A(:),24,'filled','MarkerFaceColor',[65,105,225]/255);
s2 = scatter(learn_perf_diff_B(:),learn_SCE_min_diff_B(:),24,'filled','MarkerFaceColor',[220,20,60]/255);
%0 rate line
plot([-0.5 1.0], [0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.3);
%0 performance line
plot([0 0], [-20 35],'--','Color',[0.5 0.5 0.5],'LineWidth',0.3);

%plot confidence interval
%plot(Xnew, YCI, '--k','LineWidth',1.5)

[hl,hp] = boundedline(Xnew, Ypred, abs(YCI-Ypred),'transparency', 0.08,'alpha','cmap',[139,0,139]./255);
%ob = outlinebounds(hl,hp)
%set(ob, 'linestyle', ':', 'color', [139,0,139]./255, 'marker', '.');
%plot predicted line
plot(Xnew, Ypred,'-','Color',[139,0,139]./255,'LineWidth',2)

set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)
legend([s1 s2],{'A','B'},'Location','southeast')



%% Plot
figure
subplot(1,2,1)
hold on
ylim([0 50])
plot(SCE_min_learn')
plot(mean(SCE_min_learn,1),'k--')

subplot(1,2,2)
hold on
ylim([0 50])
plot(SCE_min_recall')
plot(mean(SCE_min_recall,1),'k--')
end

