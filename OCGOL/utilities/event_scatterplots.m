function [] = event_scatterplots(AUC_min,tunedLogical)

%% Assign variables

onlyA_tuned = tunedLogical.si.onlyA_tuned;
onlyB_tuned = tunedLogical.si.onlyB_tuned;

AorB_tuned = tunedLogical.si.AorB_tuned;
AandB_tuned = tunedLogical.si.AandB_tuned;

%% Plot AUC/min scatterplots for different tuned neuron subclasses

figure;
%all neurons
subplot(1,4,1)
hold on
title('AUC/min - all neurons')
axis square
xlim([0 6]);
ylim([0 6]);
scatter(AUC_min{1},AUC_min{2}, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.7 0.7 0.7])
xlabel('A');
ylabel('B');
scatter(mean(AUC_min{1}),mean(AUC_min{2}),'MarkerFaceColor','r')
%plot center line (slope =1)
plot([0 6],[0 6],'Color',[0.5 0.5 0.5], 'LineStyle','--');

%A or B tuned
subplot(1,4,2)
hold on
title('AUC/min - A or B tuned')
axis square
xlim([0 6]);
ylim([0 6]);
scatter(AUC_min{1}(AorB_tuned),AUC_min{2}(AorB_tuned), 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.7 0.7 0.7])
xlabel('A');
ylabel('B');
scatter(mean(AUC_min{1}(AorB_tuned)),mean(AUC_min{2}(AorB_tuned)),'MarkerFaceColor','r')
%plot center line (slope =1)
plot([0 6],[0 6],'Color',[0.5 0.5 0.5], 'LineStyle','--');

%A tuned
subplot(1,4,3)
hold on
title('AUC/min - A tuned only')
axis square
xlim([0 6]);
ylim([0 6]);
scatter(AUC_min{1}(onlyA_tuned),AUC_min{2}(onlyA_tuned), 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.7 0.7 0.7])
xlabel('A');
ylabel('B');
scatter(mean(AUC_min{1}(onlyA_tuned)),mean(AUC_min{2}(onlyA_tuned)),'MarkerFaceColor','r')
%plot center line (slope =1)
plot([0 6],[0 6],'Color',[0.5 0.5 0.5], 'LineStyle','--');

%B tuned
subplot(1,4,4)
hold on
title('AUC/min - B tuned only')
axis square
xlim([0 6]);
ylim([0 6]);
scatter(AUC_min{1}(onlyB_tuned),AUC_min{2}(onlyB_tuned), 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.7 0.7 0.7])
xlabel('A');
ylabel('B');
scatter(mean(AUC_min{1}(onlyB_tuned)),mean(AUC_min{2}(onlyB_tuned)),'MarkerFaceColor','r')
%plot center line (slope =1)
plot([0 6],[0 6],'Color',[0.5 0.5 0.5], 'LineStyle','--');

end

