function [outputArg1,outputArg2] = recurrence_cum_analysis_learn(short_term_learn)

%% Input parameters
animal_nb = 6;
for aa =1:6
    nb_ses_per_animal(aa) = size(short_term_learn.perf{aa}.ses_lap_ct,2);
end

%% Short term recall - recurrence since session 1 (session read)

%create blank matrix filled with nans - animal nb x max_ses dimension
recurr_ex_short.ts.A = nan(animal_nb, max(nb_ses_per_animal));
recurr_ex_short.ts.B = nan(animal_nb, max(nb_ses_per_animal));
recurr_ex_short.si.A = nan(animal_nb, max(nb_ses_per_animal));
recurr_ex_short.si.B = nan(animal_nb, max(nb_ses_per_animal));
frac_active_short = nan(animal_nb, max(nb_ses_per_animal));

for aa=1:6
    %TS
    recurr_ex_short.ts.A(aa,1:nb_ses_per_animal(aa)) = short_term_learn.recurr{aa}.recurr_ex.ts.A(1,1:nb_ses_per_animal(aa));
    recurr_ex_short.ts.B(aa,1:nb_ses_per_animal(aa)) = short_term_learn.recurr{aa}.recurr_ex.ts.B(1,1:nb_ses_per_animal(aa));
    
    %SI
    recurr_ex_short.si.A(aa,1:nb_ses_per_animal(aa)) = short_term_learn.recurr{aa}.recurr_ex.si.A(1,1:nb_ses_per_animal(aa));
    recurr_ex_short.si.B(aa,1:nb_ses_per_animal(aa)) = short_term_learn.recurr{aa}.recurr_ex.si.B(1,1:nb_ses_per_animal(aa));

    %fraction active
    frac_active_short(aa,1:nb_ses_per_animal(aa)) = short_term_learn.recurr{aa}.frac_active_ex.first_ses(1,1:nb_ses_per_animal(aa));
end

%% Exclude sessions with less than or equal to 5 laps in A or B or low laps
%overall; exclude those with field shift as well

%animal 1, session 3 and 4 trial A/B both trials - fewer cell, only slight
%field shift (consider leaving)
recurr_ex_short.ts.A(1,[3,4]) = nan;
recurr_ex_short.ts.B(1,[3,4]) = nan;
recurr_ex_short.si.A(1,[3,4]) = nan;
recurr_ex_short.si.B(1,[3,4]) = nan;
%animal 3,  session 2, both trials -  definite shift
recurr_ex_short.ts.A(3,[2]) = nan;
recurr_ex_short.ts.B(3,[2]) = nan;
recurr_ex_short.si.A(3,[2]) = nan;
recurr_ex_short.si.B(3,[2]) = nan;
%animal 5 , session 6 and 7, both trials - few laps overall
recurr_ex_short.ts.A(5,[6,7]) = nan;
recurr_ex_short.ts.B(5,[6,7]) = nan;
recurr_ex_short.si.A(5,[6,7]) = nan;
recurr_ex_short.si.B(5,[6,7]) = nan;

%remove from active sessions
frac_active_short(1,[3,4]) = nan;
frac_active_short(3,[2]) = nan;
frac_active_short(5,[6,7]) = nan;

%% Reassemble recurrence matrix and fraction active by days (rather than sessions)

%for each day
for dd=1:9
    switch dd
        case 1
            for aa=1:6
                recur_day_order.si.A(aa,1) = recurr_ex_short.si.A(aa,1);
                recur_day_order.si.B(aa,1) = recurr_ex_short.si.B(aa,1);
                
                recur_day_order.ts.A(aa,1) = recurr_ex_short.ts.A(aa,1);
                recur_day_order.ts.B(aa,1) = recurr_ex_short.ts.B(aa,1);
                
                active_day_order(aa,1) = frac_active_short(aa,1);
            end
        case 2
            for aa=1:6
                recur_day_order.si.A(aa,2) = recurr_ex_short.si.A(aa,2);
                recur_day_order.si.B(aa,2) = recurr_ex_short.si.B(aa,2);

                recur_day_order.ts.A(aa,2) = recurr_ex_short.ts.A(aa,2);
                recur_day_order.ts.B(aa,2) = recurr_ex_short.ts.B(aa,2);
                
                active_day_order(aa,2) = frac_active_short(aa,2);
            end
        case 3
            for aa=1:6
                recur_day_order.si.A(aa,3) = recurr_ex_short.si.A(aa,3);
                recur_day_order.si.B(aa,3) = recurr_ex_short.si.B(aa,3);

                recur_day_order.ts.A(aa,3) = recurr_ex_short.ts.A(aa,3);
                recur_day_order.ts.B(aa,3) = recurr_ex_short.ts.B(aa,3);
                
                active_day_order(aa,3) = frac_active_short(aa,3);
            end
        case 4
            for aa=1:6
                recur_day_order.si.A(aa,4) = recurr_ex_short.si.A(aa,4);
                recur_day_order.si.B(aa,4) = recurr_ex_short.si.B(aa,4);

                recur_day_order.ts.A(aa,4) = recurr_ex_short.ts.A(aa,4);
                recur_day_order.ts.B(aa,4) = recurr_ex_short.ts.B(aa,4);
                
                active_day_order(aa,4) = frac_active_short(aa,4);
            end            
        case 5
            for aa=1:6
                if aa==1
                    recur_day_order.si.A(aa,5) = nan;
                    recur_day_order.si.B(aa,5) = nan;
                    
                    recur_day_order.ts.A(aa,5) = nan;
                    recur_day_order.ts.B(aa,5) = nan;
                    
                    active_day_order(aa,5) = nan;
                    
                elseif aa==2
                    recur_day_order.si.A(aa,5) = nan;
                    recur_day_order.si.B(aa,5) = nan;
                    
                    recur_day_order.ts.A(aa,5) = nan;
                    recur_day_order.ts.B(aa,5) = nan;
                    
                    active_day_order(aa,5) = nan;
                    
                elseif aa==3
                    recur_day_order.si.A(aa,5) =  recurr_ex_short.si.A(aa,5);
                    recur_day_order.si.B(aa,5) =  recurr_ex_short.si.B(aa,5);
                    
                    recur_day_order.ts.A(aa,5) =  recurr_ex_short.ts.A(aa,5);
                    recur_day_order.ts.B(aa,5) =  recurr_ex_short.ts.B(aa,5);
                    
                    active_day_order(aa,5) = frac_active_short(aa,5);
                    
                elseif aa==4
                    recur_day_order.si.A(aa,5) =  recurr_ex_short.si.A(aa,5);
                    recur_day_order.si.B(aa,5) =  recurr_ex_short.si.B(aa,5);
                    
                    recur_day_order.ts.A(aa,5) =  recurr_ex_short.ts.A(aa,5);
                    recur_day_order.ts.B(aa,5) =  recurr_ex_short.ts.B(aa,5);
                    
                    active_day_order(aa,5) = frac_active_short(aa,5);
                    
                elseif aa ==5
                    recur_day_order.si.A(aa,5) =  recurr_ex_short.si.A(aa,5);
                    recur_day_order.si.B(aa,5) =  recurr_ex_short.si.B(aa,5);
                    
                    recur_day_order.ts.A(aa,5) =  recurr_ex_short.ts.A(aa,5);
                    recur_day_order.ts.B(aa,5) =  recurr_ex_short.ts.B(aa,5);
                    
                    active_day_order(aa,5) = frac_active_short(aa,5);
                    
                elseif aa==6
                    recur_day_order.si.A(aa,5) =  recurr_ex_short.si.A(aa,5);
                    recur_day_order.si.B(aa,5) =  recurr_ex_short.si.B(aa,5);
                    
                    recur_day_order.ts.A(aa,5) =  recurr_ex_short.ts.A(aa,5);
                    recur_day_order.ts.B(aa,5) =  recurr_ex_short.ts.B(aa,5);   
                    
                    active_day_order(aa,5) = frac_active_short(aa,5);
                end
            end
        case 6
            for aa=1:6
                %for each animal
                if aa==1
                    recur_day_order.si.A(aa,6) = recurr_ex_short.si.A(aa,5);
                    recur_day_order.si.B(aa,6) = recurr_ex_short.si.B(aa,5);
                    
                    recur_day_order.ts.A(aa,6) = recurr_ex_short.ts.A(aa,5);
                    recur_day_order.ts.B(aa,6) = recurr_ex_short.ts.B(aa,5);
                    
                    active_day_order(aa,6) = frac_active_short(aa,5);
                elseif aa==2
                    recur_day_order.si.A(aa,6) = recurr_ex_short.si.A(aa,5);
                    recur_day_order.si.B(aa,6) = recurr_ex_short.si.B(aa,5);
                    
                    recur_day_order.ts.A(aa,6) = recurr_ex_short.ts.A(aa,5);
                    recur_day_order.ts.B(aa,6) = recurr_ex_short.ts.B(aa,5);   
                    
                    active_day_order(aa,6) = frac_active_short(aa,5);
                    
                elseif aa==3
                    recur_day_order.si.A(aa,6) =  recurr_ex_short.si.A(aa,6);
                    recur_day_order.si.B(aa,6) =  recurr_ex_short.si.B(aa,6);
                    
                    recur_day_order.ts.A(aa,6) =  recurr_ex_short.ts.A(aa,6);
                    recur_day_order.ts.B(aa,6) =  recurr_ex_short.ts.B(aa,6); 
                    
                    active_day_order(aa,6) = frac_active_short(aa,6);
                    
                elseif aa==4
                    recur_day_order.si.A(aa,6) =  recurr_ex_short.si.A(aa,6);
                    recur_day_order.si.B(aa,6) =  recurr_ex_short.si.B(aa,6);
                    
                    recur_day_order.ts.A(aa,6) =  recurr_ex_short.ts.A(aa,6);
                    recur_day_order.ts.B(aa,6) =  recurr_ex_short.ts.B(aa,6);
                    
                    active_day_order(aa,6) = frac_active_short(aa,6);
                    
                elseif aa ==5
                    recur_day_order.si.A(aa,6) =  recurr_ex_short.si.A(aa,6);
                    recur_day_order.si.B(aa,6) =  recurr_ex_short.si.B(aa,6);
                    
                    recur_day_order.ts.A(aa,6) =  recurr_ex_short.ts.A(aa,6);
                    recur_day_order.ts.B(aa,6) =  recurr_ex_short.ts.B(aa,6);
                    
                    active_day_order(aa,6) = frac_active_short(aa,6);
                    
                elseif aa==6
                    recur_day_order.si.A(aa,6) =  recurr_ex_short.si.A(aa,6);
                    recur_day_order.si.B(aa,6) =  recurr_ex_short.si.B(aa,6);
                    
                    recur_day_order.ts.A(aa,6) =  recurr_ex_short.ts.A(aa,6);
                    recur_day_order.ts.B(aa,6) =  recurr_ex_short.ts.B(aa,6);
                    
                    active_day_order(aa,6) = frac_active_short(aa,6);
                end
            end
            
        case 7
            for aa=1:6
                %for each animal
                if aa==1
                    recur_day_order.si.A(aa,7) = recurr_ex_short.si.A(aa,6);
                    recur_day_order.si.B(aa,7) = recurr_ex_short.si.B(aa,6);
                    
                    recur_day_order.ts.A(aa,7) = recurr_ex_short.ts.A(aa,6);
                    recur_day_order.ts.B(aa,7) = recurr_ex_short.ts.B(aa,6);
                    
                    active_day_order(aa,7) = frac_active_short(aa,6);
                    
                elseif aa==2
                    recur_day_order.si.A(aa,7) = recurr_ex_short.si.A(aa,6);
                    recur_day_order.si.B(aa,7) = recurr_ex_short.si.B(aa,6);
                    
                    recur_day_order.ts.A(aa,7) = recurr_ex_short.ts.A(aa,6);
                    recur_day_order.ts.B(aa,7) = recurr_ex_short.ts.B(aa,6);
                    
                    active_day_order(aa,7) = frac_active_short(aa,6);
                elseif aa==3
                    recur_day_order.si.A(aa,7) = recurr_ex_short.si.A(aa,7);
                    recur_day_order.si.B(aa,7) = recurr_ex_short.si.B(aa,7);
                    
                    recur_day_order.ts.A(aa,7) = recurr_ex_short.ts.A(aa,7);
                    recur_day_order.ts.B(aa,7) = recurr_ex_short.ts.B(aa,7);
                    
                    active_day_order(aa,7) = frac_active_short(aa,7);
                elseif aa==4
                    recur_day_order.si.A(aa,7) = recurr_ex_short.si.A(aa,7);
                    recur_day_order.si.B(aa,7) = recurr_ex_short.si.B(aa,7);
                    
                    recur_day_order.ts.A(aa,7) = recurr_ex_short.ts.A(aa,7);
                    recur_day_order.ts.B(aa,7) = recurr_ex_short.ts.B(aa,7);
                    
                    active_day_order(aa,7) = frac_active_short(aa,7);
                elseif aa ==5
                    recur_day_order.si.A(aa,7) =  nan;
                    recur_day_order.si.B(aa,7) =  nan;
                    
                    recur_day_order.ts.A(aa,7) =  nan;
                    recur_day_order.ts.B(aa,7) =  nan;
                    
                    active_day_order(aa,7) = nan;
                    
                elseif aa==6
                    recur_day_order.si.A(aa,7) =  nan;
                    recur_day_order.si.B(aa,7) =  nan;
                    
                    recur_day_order.ts.A(aa,7) =  nan;
                    recur_day_order.ts.B(aa,7) =  nan;
                    
                    active_day_order(aa,7) = nan;
                end
            end            
        case 8
            for aa=1:6
                %for each animal
                if aa==1
                    recur_day_order.si.A(aa,8) = recurr_ex_short.si.A(aa,7);
                    recur_day_order.si.B(aa,8) = recurr_ex_short.si.B(aa,7);
                    
                    recur_day_order.ts.A(aa,8) = recurr_ex_short.ts.A(aa,7);
                    recur_day_order.ts.B(aa,8) = recurr_ex_short.ts.B(aa,7);
                    
                    active_day_order(aa,8) = frac_active_short(aa,7);
                    
                elseif aa==2
                    recur_day_order.si.A(aa,8) = recurr_ex_short.si.A(aa,7);
                    recur_day_order.si.B(aa,8) = recurr_ex_short.si.B(aa,7);
                    
                    recur_day_order.ts.A(aa,8) = recurr_ex_short.ts.A(aa,7);
                    recur_day_order.ts.B(aa,8) = recurr_ex_short.ts.B(aa,7);
                    
                    active_day_order(aa,8) = frac_active_short(aa,7);
                    
                elseif aa==3
                    recur_day_order.si.A(aa,8) = recurr_ex_short.si.A(aa,8);
                    recur_day_order.si.B(aa,8) = recurr_ex_short.si.B(aa,8);
                    
                    recur_day_order.ts.A(aa,8) = recurr_ex_short.ts.A(aa,8);
                    recur_day_order.ts.B(aa,8) = recurr_ex_short.ts.B(aa,8);
                    
                    active_day_order(aa,8) = frac_active_short(aa,8);
                    
                elseif aa==4
                    recur_day_order.si.A(aa,8) = recurr_ex_short.si.A(aa,8);
                    recur_day_order.si.B(aa,8) = recurr_ex_short.si.B(aa,8);
                    
                    recur_day_order.ts.A(aa,8) = recurr_ex_short.ts.A(aa,8);
                    recur_day_order.ts.B(aa,8) = recurr_ex_short.ts.B(aa,8);
                    
                    active_day_order(aa,8) = frac_active_short(aa,8);
                    
                elseif aa ==5
                    recur_day_order.si.A(aa,8) = recurr_ex_short.si.A(aa,7);
                    recur_day_order.si.B(aa,8) = recurr_ex_short.si.B(aa,7);
                    
                    recur_day_order.ts.A(aa,8) = recurr_ex_short.ts.A(aa,7);
                    recur_day_order.ts.B(aa,8) = recurr_ex_short.ts.B(aa,7);
                    
                    active_day_order(aa,8) = frac_active_short(aa,7);
                    
                elseif aa==6
                    recur_day_order.si.A(aa,8) =  nan;
                    recur_day_order.si.B(aa,8) =  nan;
                    
                    recur_day_order.ts.A(aa,8) =  nan;
                    recur_day_order.ts.B(aa,8) =  nan;
                    
                    active_day_order(aa,8) = nan;
                end
            end               
        case 9
            for aa=1:6
                %for each animal
                if aa==1
                    recur_day_order.si.A(aa,9) = nan;
                    recur_day_order.si.B(aa,9) = nan;
                    
                    recur_day_order.ts.A(aa,9) = nan;
                    recur_day_order.ts.B(aa,9) = nan;
                    
                    active_day_order(aa,9) = nan;
                    
                elseif aa==2
                    recur_day_order.si.A(aa,9) = nan;
                    recur_day_order.si.B(aa,9) = nan;
                   
                    recur_day_order.ts.A(aa,9) = nan;
                    recur_day_order.ts.B(aa,9) = nan;
                    
                    active_day_order(aa,9) = nan;
                    
                elseif aa==3
                    recur_day_order.si.A(aa,9) = nan;
                    recur_day_order.si.B(aa,9) = nan;
                    
                    recur_day_order.ts.A(aa,9) = nan;
                    recur_day_order.ts.B(aa,9) = nan;
                    
                    active_day_order(aa,9) = nan;
                elseif aa==4
                    recur_day_order.si.A(aa,9) = recurr_ex_short.si.A(aa,9);
                    recur_day_order.si.B(aa,9) = recurr_ex_short.si.B(aa,9);
                    
                    recur_day_order.ts.A(aa,9) = recurr_ex_short.ts.A(aa,9);
                    recur_day_order.ts.B(aa,9) = recurr_ex_short.ts.B(aa,9);
                    
                    active_day_order(aa,9) = frac_active_short(aa,9);
                    
                elseif aa ==5
                    recur_day_order.si.A(aa,9) = nan;
                    recur_day_order.si.B(aa,9) = nan;
                    
                    recur_day_order.ts.A(aa,9) = nan;
                    recur_day_order.ts.B(aa,9) = nan;       
                    
                    active_day_order(aa,9) = nan;
                elseif aa==6
                    recur_day_order.si.A(aa,9) = nan;
                    recur_day_order.si.B(aa,9) = nan;
                    
                    recur_day_order.ts.A(aa,9) = nan;
                    recur_day_order.ts.B(aa,9) = nan;
                    
                    active_day_order(aa,9) = nan;
                end
            end             
            
    end
end

%% Get mean and sem and plot - short term recall (day order matrix)
recurr_short_mean.ts.A = nanmean(recur_day_order.ts.A,1);
recurr_short_mean.ts.B = nanmean(recur_day_order.ts.B,1);

recurr_short_mean.si.A = nanmean(recur_day_order.si.A,1);
recurr_short_mean.si.B = nanmean(recur_day_order.si.B,1);

%number of points for sem calculation - use either si or ts, but check both
%A and B trials
nb_sem_ST_recall.A = sum(~isnan(recur_day_order.si.A),1);
nb_sem_ST_recall.B = sum(~isnan(recur_day_order.ts.B),1);

%get sem from std using numbers above
%A
recurr_short_sem.ts.A = nanstd(recur_day_order.ts.A,0,1)./sqrt(nb_sem_ST_recall.A);
recurr_short_sem.si.A = nanstd(recur_day_order.si.A,0,1)./sqrt(nb_sem_ST_recall.A);
%B 
recurr_short_sem.ts.B = nanstd(recur_day_order.ts.B,0,1)./sqrt(nb_sem_ST_recall.B);
recurr_short_sem.si.B = nanstd(recur_day_order.si.B,0,1)./sqrt(nb_sem_ST_recall.B);

%% Get mean sem of frac active (short term recall)
mean_active.ST.recall = nanmean(active_day_order,1);
nb_sem_active.ST.recall = sum(~isnan(active_day_order),1);
sem_active.ST.recall = nanstd(active_day_order,0,1)./sqrt(nb_sem_active.ST.recall);


%% Plot SHORT term recall recurrence and fraction active - TS

figure
hold on 
title('Recurrence Learning - S.I.')
ylim([0 1.2])
xlim([0 10])
xticks(1:9)
xticklabels({'1','2','3','4','5','6','7','8','9'})
yticks(0:0.2:1)
xlabel('Days')
ylabel('Recurrence probability')
%Recall
rA = errorbar(1:9,recurr_short_mean.si.A,recurr_short_sem.si.A,'LineStyle','--','Linewidth',2,'Color',[65,105,225]/255);
rB = errorbar(1:9,recurr_short_mean.si.B,recurr_short_sem.si.B,'LineStyle','--','Linewidth',2,'Color',[220,20,60]/255);
%Active cells
ac = errorbar(1:9,mean_active.ST.recall,sem_active.ST.recall,'LineStyle','-','Linewidth',2,'Color',[34,139,34]/255);
%lA = errorbar(1:8,mean_mean_PV_learn_neighbor_all.A,sem_PV_learn_neighbor_all.A,'LineStyle','--','Linewidth',2,'Color',[65,105,225]/255);
%lB = errorbar(1:8,mean_mean_PV_learn_neighbor_all.B,sem_PV_learn_neighbor_all.B,'LineStyle','--','Linewidth',2,'Color',[220,20,60]/255);
set(gca,'FontSize',16)
set(gca,'Linewidth',2)

legend([rA,rB, ac],{'Learning A','Learning B','Active cells'},'location','northeast')

%% Recurrence learning T.S.

figure
hold on 
title('Recurrence Learning - T.S.')
ylim([0 1.2])
xlim([0 10])
xticks(1:9)
xticklabels({'1','2','3','4','5','6','7','8','9'})
yticks(0:0.2:1)
xlabel('Days')
ylabel('Recurrence probability')
%Recall
rA = errorbar(1:9,recurr_short_mean.ts.A,recurr_short_sem.ts.A,'LineStyle','--','Linewidth',2,'Color',[65,105,225]/255);
rB = errorbar(1:9,recurr_short_mean.ts.B,recurr_short_sem.ts.B,'LineStyle','--','Linewidth',2,'Color',[220,20,60]/255);
%Active cells
ac = errorbar(1:9,mean_active.ST.recall,sem_active.ST.recall,'LineStyle','-','Linewidth',2,'Color',[34,139,34]/255);
%lA = errorbar(1:8,mean_mean_PV_learn_neighbor_all.A,sem_PV_learn_neighbor_all.A,'LineStyle','--','Linewidth',2,'Color',[65,105,225]/255);
%lB = errorbar(1:8,mean_mean_PV_learn_neighbor_all.B,sem_PV_learn_neighbor_all.B,'LineStyle','--','Linewidth',2,'Color',[220,20,60]/255);
set(gca,'FontSize',16)
set(gca,'Linewidth',2)

legend([rA,rB, ac],{'Learning A','Learning B','Active cells'},'location','northeast')


end

