%% Load in PV and TC correlation data from recall animals

%cross session directories (recall experiments)
cross_dirs_recall = {'G:\OCGOL_stability_recall\I46\crossSession',...
    'G:\OCGOL_stability_recall\I45_RT\crossSession',...
    'G:\OCGOL_stability_recall\I42L_1\crossSession',...
    'G:\OCGOL_stability_recall\I42R_1\crossSession'};

%cross session directories (recall experiments)
cross_dirs_learning = {'G:\OCGOL_learning_short_term\I56_RTLS\crossSession',...
    'G:\OCGOL_learning_short_term\I57_RTLS\crossSession'};

%read in recall data
for ss=1:size(cross_dirs_recall,2)
    PV_TC_corr_recall(ss) = load(fullfile(cross_dirs_recall{ss},'PV_TC_corr.mat'));
end

%read in recall data
for ss=1:size(cross_dirs_learning,2)
 PV_TC_corr_learning(ss) = load(fullfile(cross_dirs_learning{ss},'PV_TC_corr.mat'));
end

%% Plot learning and recall TC correlation relative to day 1 on same plot

figure;
subplot(1,2,1)
hold on
title('TC A')
ylim([0 1])
%recall data
for ss = 1:4
    plot([2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanTC_rel_d1.ts.A,'g-')
end

%learning data
for ss = 1:2
    plot([2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanTC_rel_d1.ts.A,'r--')
end

subplot(1,2,2)
hold on
title('TC B')
ylim([0 1])
%recall data
for ss = 1:4
    plot([2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanTC_rel_d1.ts.B,'g-')
end

%learning data
for ss = 1:2
    plot([2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanTC_rel_d1.ts.B,'r--')
end

%% Plot learning and recall PV correlation relative to day 1 on same plot

figure;
subplot(1,2,1)
hold on
title('PV A')
ylim([0 1])
%recall data
for ss = 1:4
    plot([2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanPV_rel_d1.A,'g-')
end

%learning data
for ss = 1:2
    plot([2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanPV_rel_d1.A,'r--')
end

subplot(1,2,2)
hold on
title('PV B')
ylim([0 1])
%recall data
for ss = 1:4
    plot([2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanPV_rel_d1.B,'g-')
end

%learning data
for ss = 1:2
    plot([2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanPV_rel_d1.B,'r--')
end


%% Same day PV

figure;
hold on
title('Same day PV - at least 1 match to another session')
ylim([0 1])
%recall data
for ss = 1:4
    plot([1 2 3 6 7 8 9],PV_TC_corr_recall(ss).PV_TC_corr.meanPV_same_day.AB,'g-')
end

%learning data
for ss = 1:2
    plot([1 2 3 4 5 6],PV_TC_corr_learning(ss).PV_TC_corr.meanPV_same_day.AB,'r--')
end
