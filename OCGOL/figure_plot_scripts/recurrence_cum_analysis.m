function [outputArg1,outputArg2] = recurrence_cum_analysis(short_term_recall,long_term_recall)


%% Short term recall - recurrence since day 1
for aa=1:5
    %TS
    recurr_ex_short.ts.A(aa,:) = short_term_recall.recurr{aa}.recurr_ex.ts.A(1,:);
    recurr_ex_short.ts.B(aa,:) = short_term_recall.recurr{aa}.recurr_ex.ts.B(1,:);
    
    %SI
    recurr_ex_short.si.A(aa,:) = short_term_recall.recurr{aa}.recurr_ex.si.A(1,:);
    recurr_ex_short.si.B(aa,:) = short_term_recall.recurr{aa}.recurr_ex.si.B(1,:);

    %fraction active
    frac_active_short(aa,:) = short_term_recall.recurr{aa}.frac_active_ex.first_ses(1,:);
end

%Exclude sessions with less than or equal to 5 laps in A or B
%animal 2 , trial B, session 3 and 7
recurr_ex_short.ts.B(2,[3,7]) = nan;
recurr_ex_short.si.B(2,[3,7]) = nan;

%% Long term recall - recurrence since day 1
for aa=1:3
    %TS
    recurr_ex_long.ts.A(aa,:) = long_term_recall.recurr{aa}.recurr_ex.ts.A(1,:);
    recurr_ex_long.ts.B(aa,:) = long_term_recall.recurr{aa}.recurr_ex.ts.B(1,:);
    
    %SI
    recurr_ex_long.si.A(aa,:) = long_term_recall.recurr{aa}.recurr_ex.si.A(1,:);
    recurr_ex_long.si.B(aa,:) = long_term_recall.recurr{aa}.recurr_ex.si.B(1,:);

    %fraction active
    frac_active_long(aa,:) = long_term_recall.recurr{aa}.frac_active_ex.first_ses(1,:);
end

%Exclude sessions with less than or equal to 5 laps in A or B
%animal 3 , both A/B due to field shift, session 2 3 4 
recurr_ex_long.ts.A(3,[2,3,4]) = nan;
recurr_ex_long.si.A(3,[2,3,4]) = nan;

recurr_ex_long.ts.B(3,[2,3,4]) = nan;
recurr_ex_long.si.B(3,[2,3,4]) = nan;

%exclude long as well
frac_active_long(3,[2 3 4]) = nan;

%% Get mean and sem and plot - short term recall
recurr_short_mean.ts.A = nanmean(recurr_ex_short.ts.A,1);
recurr_short_mean.ts.B = nanmean(recurr_ex_short.ts.B,1);

recurr_short_mean.si.A = nanmean(recurr_ex_short.si.A,1);
recurr_short_mean.si.B = nanmean(recurr_ex_short.si.B,1);

%number of points for sem calculation - use either si or ts, but check both
%A and B trials
nb_sem_ST_recall.A = sum(~isnan(recurr_ex_short.ts.A),1);
nb_sem_ST_recall.B = sum(~isnan(recurr_ex_short.ts.B),1);

%get sem from std using numbers above
recurr_short_sem.ts.A = nanstd(recurr_ex_short.ts.A,0,1)./sqrt(nb_sem_ST_recall.A);
recurr_short_sem.si.A = nanstd(recurr_ex_short.si.A,0,1)./sqrt(nb_sem_ST_recall.A);

recurr_short_sem.ts.B = nanstd(recurr_ex_short.ts.B,0,1)./sqrt(nb_sem_ST_recall.B);
recurr_short_sem.si.B = nanstd(recurr_ex_short.si.B,0,1)./sqrt(nb_sem_ST_recall.B);

%% Get mean sem of frac active (short term recall)
mean_active.ST.recall = nanmean(frac_active_short,1);
nb_sem_active.ST.recall = sum(~isnan(frac_active_short),1);
sem_active.ST.recall = nanstd(frac_active_short,0,1)./sqrt(nb_sem_active.ST.recall);


%% Get mean and sem and plot
recurr_long_mean.ts.A = nanmean(recurr_ex_long.ts.A,1);
recurr_long_mean.ts.B = nanmean(recurr_ex_long.ts.B,1);

recurr_long_mean.si.A = nanmean(recurr_ex_long.si.A,1);
recurr_long_mean.si.B = nanmean(recurr_ex_long.si.B,1);

%number of points for sem calculation - use either si or ts, but check both
%A and B trials
nb_sem_LT_recall.A = sum(~isnan(recurr_ex_long.ts.A),1);
nb_sem_LT_recall.B = sum(~isnan(recurr_ex_long.ts.B),1);

%get sem from std using numbers above
recurr_long_sem.ts.A = nanstd(recurr_ex_long.ts.A,0,1)./sqrt(nb_sem_LT_recall.A);
recurr_long_sem.si.A = nanstd(recurr_ex_long.si.A,0,1)./sqrt(nb_sem_LT_recall.A);

recurr_long_sem.ts.B = nanstd(recurr_ex_long.ts.B,0,1)./sqrt(nb_sem_LT_recall.B);
recurr_long_sem.si.B = nanstd(recurr_ex_long.si.B,0,1)./sqrt(nb_sem_LT_recall.B);

%% Get mean sem of frac active (longterm recall)
mean_active.LT.recall = nanmean(frac_active_long,1);
nb_sem_active.LT.recall = sum(~isnan(frac_active_long),1);
sem_active.LT.recall = nanstd(frac_active_long,0,1)./sqrt(nb_sem_active.LT.recall);


%% Plot SHORT term recall recurrence and fraction active - TS
figure
hold on 
title('Recurrence - T.S.')
ylim([0 1.2])
xlim([0 8])
xticks(1:7)
xticklabels({'1','2','3','6','7','8','9'})
yticks(0:0.2:1)
xlabel('Days')
ylabel('Recurrence probability')
%Recall
rA = errorbar(1:7,recurr_short_mean.ts.A,recurr_short_sem.ts.A,'LineStyle','-','Linewidth',2,'Color',[65,105,225]/255);
rB = errorbar(1:7,recurr_short_mean.ts.B,recurr_short_sem.ts.B,'LineStyle','-','Linewidth',2,'Color',[220,20,60]/255);
%Active cells
ac = errorbar(1:7,mean_active.ST.recall,sem_active.ST.recall,'LineStyle','-','Linewidth',2,'Color',[34,139,34]/255);
%lA = errorbar(1:8,mean_mean_PV_learn_neighbor_all.A,sem_PV_learn_neighbor_all.A,'LineStyle','--','Linewidth',2,'Color',[65,105,225]/255);
%lB = errorbar(1:8,mean_mean_PV_learn_neighbor_all.B,sem_PV_learn_neighbor_all.B,'LineStyle','--','Linewidth',2,'Color',[220,20,60]/255);
set(gca,'FontSize',16)
set(gca,'Linewidth',2)

%legend([lA,lB,rA,rB],{'Learning A','Learning B','Recall A','Recall B'},'location','southeast')

figure
hold on 
title('Recurrence - S.I.')
ylim([0 1.2])
xlim([0 8])
xticks(1:7)
xticklabels({'1','2','3','6','7','8','9'})
yticks(0:0.2:1)
xlabel('Days')
ylabel('Recurrence probability')
%Recall
rA = errorbar(1:7,recurr_short_mean.si.A,recurr_short_sem.si.A,'LineStyle','-','Linewidth',2,'Color',[65,105,225]/255);
rB = errorbar(1:7,recurr_short_mean.si.B,recurr_short_sem.si.B,'LineStyle','-','Linewidth',2,'Color',[220,20,60]/255);
%Active cells
ac = errorbar(1:7,mean_active.ST.recall,sem_active.ST.recall,'LineStyle','-','Linewidth',2,'Color',[34,139,34]/255);
%lA = errorbar(1:8,mean_mean_PV_learn_neighbor_all.A,sem_PV_learn_neighbor_all.A,'LineStyle','--','Linewidth',2,'Color',[65,105,225]/255);
%lB = errorbar(1:8,mean_mean_PV_learn_neighbor_all.B,sem_PV_learn_neighbor_all.B,'LineStyle','--','Linewidth',2,'Color',[220,20,60]/255);
set(gca,'FontSize',16)
set(gca,'Linewidth',2)


%% Plot LONG term recall recurrence and fraction active - TS
figure
hold on 
title('Recurrence - T.S. - Long recall')
ylim([0 1.2])
xlim([0 8])
xticks(1:7)
xticklabels({'1','6','16','20','25','30'})
yticks(0:0.2:1)
xlabel('Days')
ylabel('Recurrence probability')
%Recall
rA = errorbar(1:6,recurr_long_mean.ts.A,recurr_long_sem.ts.A,'LineStyle','-','Linewidth',2,'Color',[65,105,225]/255);
rB = errorbar(1:6,recurr_long_mean.ts.B,recurr_long_sem.ts.B,'LineStyle','-','Linewidth',2,'Color',[220,20,60]/255);
%Active cells
ac = errorbar(1:6,mean_active.LT.recall,sem_active.LT.recall,'LineStyle','-','Linewidth',2,'Color',[34,139,34]/255);
%lA = errorbar(1:8,mean_mean_PV_learn_neighbor_all.A,sem_PV_learn_neighbor_all.A,'LineStyle','--','Linewidth',2,'Color',[65,105,225]/255);
%lB = errorbar(1:8,mean_mean_PV_learn_neighbor_all.B,sem_PV_learn_neighbor_all.B,'LineStyle','--','Linewidth',2,'Color',[220,20,60]/255);
set(gca,'FontSize',16)
set(gca,'Linewidth',2)

%legend([lA,lB,rA,rB],{'Learning A','Learning B','Recall A','Recall B'},'location','southeast')

figure
hold on 
title('Recurrence - S.I. - Long recall')
ylim([0 1.2])
xlim([0 8])
xticks(1:7)
xticklabels({'1','6','16','20','25','30'})
yticks(0:0.2:1)
xlabel('Days')
ylabel('Recurrence probability')
%Recall
rA = errorbar(1:6,recurr_long_mean.si.A,recurr_long_sem.si.A,'LineStyle','-','Linewidth',2,'Color',[65,105,225]/255);
rB = errorbar(1:6,recurr_long_mean.si.B,recurr_long_sem.si.B,'LineStyle','-','Linewidth',2,'Color',[220,20,60]/255);
%Active cells
ac = errorbar(1:6,mean_active.LT.recall,sem_active.LT.recall,'LineStyle','-','Linewidth',2,'Color',[34,139,34]/255);
%lA = errorbar(1:8,mean_mean_PV_learn_neighbor_all.A,sem_PV_learn_neighbor_all.A,'LineStyle','--','Linewidth',2,'Color',[65,105,225]/255);
%lB = errorbar(1:8,mean_mean_PV_learn_neighbor_all.B,sem_PV_learn_neighbor_all.B,'LineStyle','--','Linewidth',2,'Color',[220,20,60]/255);
set(gca,'FontSize',16)
set(gca,'Linewidth',2)

end

