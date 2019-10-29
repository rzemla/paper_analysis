function [outputArg1,outputArg2] = combine_STC_plot_multi_animal_short_vs_long_term(TC_corr_match_learning,TC_corr_match_recall)

%% Define number of animals and number of sessions for each animal
%define number of learning animals
nb_animal_learn = size(TC_corr_match_learning,2);
%define number of short term recall animals
nb_animal_recall = size(TC_corr_match_recall,2);

%number of sessions for learning
for aa=1:nb_animal_learn
    nb_ses_learn(aa) = size(TC_corr_match_learning{aa}.tc_corr_match.PV_corr_all_day,1);
end
%number of sessions for recall
for aa=1:nb_animal_recall
    nb_ses_recall(aa) = size(TC_corr_match_recall{aa}.tc_corr_match.PV_corr_all_day,1);
end

%% Combine relative to d1 STCs from multiple animals between any 2
%sessions/days
%% Short-term learning tuned TS STCs 
%for each learning session (define max for animal with max and leave all
%others
%9 sessions max for animal 4
for ss=1:9
    comb_STC_learning.A{1,ss} = [];
    comb_STC_learning.B{1,ss} = [];
end

%for each session

%for each animal - TS tuned neurons matching on both sessions
for aa=1:nb_animal_learn
    %check if session exists
    for ss=1:nb_ses_learn(aa)
        comb_STC_learning.A{1,ss} = [comb_STC_learning.A{1,ss} ; cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, ss}')'];
        comb_STC_learning.B{1,ss} = [comb_STC_learning.B{1,ss} ; cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, ss}')'];
    end
end


%sort by session 1 relative to session 1
for ss=1:9
    %sort for A matches
    [~,maxBin] = max(comb_STC_learning.A{1,ss}(:,1:100)', [], 1,'includenan');
    %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
    [~,sortOrder] = sort(maxBin,'ascend');
    
    comb_STC_learning_sort.A{1,ss} = comb_STC_learning.A{1,ss}(sortOrder,:);
    
    %sort for B matches
    [~,maxBin] = max(comb_STC_learning.B{1,ss}(:,1:100)', [], 1,'includenan');
    %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
    [~,sortOrder] = sort(maxBin,'ascend');
    
    comb_STC_learning_sort.B{1,ss} = comb_STC_learning.B{1,ss}(sortOrder,:);
end

%Plot TS tuned matching place cells for A and B trials
figure
for ss=2:9
    subplot(1,8,ss-1)
    imagesc(comb_STC_learning_sort.A{1,ss})
    hold on
    title('Learn - A')
    caxis([0 1])
    colormap('jet')
end

figure
for ss=2:9
    subplot(1,8,ss-1)
    imagesc(comb_STC_learning_sort.B{1,ss})
    hold on
    title('Learn - B')
    caxis([0 1])
    colormap('jet')
end

%% Recall - consecutive session alignment
for ss=1:6
    comb_STC_recall.A{1,ss} = [];
    comb_STC_recall.B{1,ss} = [];
end

for ss=1:6
    for aa=1:nb_animal_recall
        comb_STC_recall.A{1,ss} = [comb_STC_recall.A{1,ss} ; cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, ss}')'];
        comb_STC_recall.B{1,ss} = [comb_STC_recall.B{1,ss} ; cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, ss}')'];
    end
end

%sort by session 1 relative to session 1
for ss=1:6
    %sort for A matches
    [~,maxBin] = max(comb_STC_recall.A{1,ss}(:,1:100)', [], 1,'includenan');
    %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
    [~,sortOrder] = sort(maxBin,'ascend');
    
    comb_STC_recall_sort.A{1,ss} = comb_STC_recall.A{1,ss}(sortOrder,:);
    
    %sort for B matches
    [~,maxBin] = max(comb_STC_recall.B{1,ss}(:,1:100)', [], 1,'includenan');
    %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
    [~,sortOrder] = sort(maxBin,'ascend');
    
    comb_STC_recall_sort.B{1,ss} = comb_STC_recall.B{1,ss}(sortOrder,:);
end

figure
for ss=2:6
    subplot(1,5,ss-1)
    imagesc(comb_STC_recall_sort.A{1,ss})
    hold on
    title('Recall - A')
    caxis([0 1])
    colormap('jet')
end

figure
for ss=2:6
    subplot(1,5,ss-1)
    imagesc(comb_STC_recall_sort.B{1,ss})
    hold on
    title('Recall - B')
    caxis([0 1])
    colormap('jet')
end

%% Align relative to D1 (Figure 4F) - as a function of day distance - RECALL
%number of recall animals
%nb_animal_recall = 4;

%place into function
for dd=2:6
    %for each day relative to day 1 choose the corresponding imaging
    %sessions
    %for each day, ...choose sessions
    switch dd
        case 2
            for aa=1:nb_animal_recall
                comb_STC_recall_days.A{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 2}')'];
                comb_STC_recall_days.B{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 2}')'];
            end
        case 3
            for aa=1:nb_animal_recall
                comb_STC_recall_days.A{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 3}')'];
                comb_STC_recall_days.B{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 3}')'];
            end
        case 4 %equivalent day 1 distance, since not imaged on this day (use only 1 session rather than both)
            for aa=1:nb_animal_recall
                comb_STC_recall_days.A{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 4}')'];
                comb_STC_recall_days.B{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 4}')'];
            end
        case 5 %equivalent day 1 distance, since not imaged on this day (use only 1 session rather than both)
            for aa=1:nb_animal_recall
                comb_STC_recall_days.A{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 5}')'];
                comb_STC_recall_days.B{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 5}')'];
            end
        case 6
            for aa=1:nb_animal_recall
                comb_STC_recall_days.A{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 6}')'];
                comb_STC_recall_days.B{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 6}')'];
            end
            
            
    end
end

%% Align by NEIGHBORING days RECALL
%all animals treat the same
%D1 vs. D2 - ses 1 vs. 2
%D2 vs. D3 - ses 2 vs. 3
%D6 vs. D7 - ses 4 vs. 5
%D7 vs. D8 - ses 5 vs. 6
%D8 vs. D9 - ses 6 vs. 7

%run consecutive - 1 6 16 20 25 30 days (1 2 3 4 5 6 -session/pseudoday)

%place into function - number by previous day correlation
for dd=1:5
    %for each day relative to day 1 choose the corresponding imaging
    %sessions
    %for each day, ...choose sessions
    switch dd
        case 1
            for aa=1:nb_animal_recall
                neighbor_STC_recall_days.A{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 2}')'];
                neighbor_STC_recall_days.B{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 2}')'];
            end
        case 2
            for aa=1:nb_animal_recall
                neighbor_STC_recall_days.A{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{2, 3}')'];
                neighbor_STC_recall_days.B{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{2, 3}')'];
            end
        case 3 
            for aa=1:nb_animal_recall
                neighbor_STC_recall_days.A{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{3, 4}')'];
                neighbor_STC_recall_days.B{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{3, 4}')'];
            end
        case 4 
            for aa=1:nb_animal_recall
                neighbor_STC_recall_days.A{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{4, 5}')'];
                neighbor_STC_recall_days.B{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{4, 5}')'];
            end
        case 5
            for aa=1:nb_animal_recall
                neighbor_STC_recall_days.A{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{5, 6}')'];
                neighbor_STC_recall_days.B{aa,dd} = [cell2mat(TC_corr_match_recall{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{5, 6}')'];
            end
            
    end
end


%% Align relative to D1 (Figure 4F) - as a function of day distance - LEARN - updated to include all animals and sessions

%until day 9 - as far as the learning daygoes
for dd=2:9
    %for each day relative to day 1 choose the corresponding imaging
    %sessions
    %for each day, ...choose sessions
    switch dd
        case 2
            for aa=1:nb_animal_learn
                comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 2}')'];
                comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 2}')'];
            end
        case 3
            for aa=1:nb_animal_learn
                comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 3}')'];
                comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 3}')'];
            end
        case 4 %equivalent day 1 distance, since not imaged on this day (use only 1 session rather than both)
            for aa=1:nb_animal_learn
                comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 4}')'];
                comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 4}')'];
            end
        case 5 %equivalent day 1 distance, since not imaged on this day
            for aa=1:nb_animal_learn
                if aa==1
                    comb_STC_learn_days.A{aa,dd} = [];
                    comb_STC_learn_days.B{aa,dd} = [];
                elseif aa==2
                    comb_STC_learn_days.A{aa,dd} = [];
                    comb_STC_learn_days.B{aa,dd} = [];
                elseif aa >= 3 %all other animals day 5 corresponds to session 5
                comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 5}')'];
                comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 5}')'];                    
                    
                end
            end
        case 6
            for aa=1:nb_animal_learn
                if aa==1
                    comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 5}')'];
                    comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 5}')'];
                elseif aa==2
                    comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 5}')'];
                    comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 5}')'];
                elseif aa==3
                    comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 6}')'];
                    comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 6}')'];
                elseif aa>=4 %all other animals greater than or equal to 4
                    comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 6}')'];
                    comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 6}')'];
                end
            end
            
        case 7
            for aa=1:nb_animal_learn
                if aa==1
                    comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 6}')'];
                    comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 6}')'];
                elseif aa==2
                    comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 6}')'];
                    comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 6}')'];
                elseif aa==3
                    comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 7}')'];
                    comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 7}')'];
                elseif aa==4
                    comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 7}')'];
                    comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 7}')'];
                elseif aa==5 %no day 7
                    comb_STC_learn_days.A{aa,dd} = [];
                    comb_STC_learn_days.B{aa,dd} = [];
                elseif aa==6 %no day 7
                    comb_STC_learn_days.A{aa,dd} = [];
                    comb_STC_learn_days.B{aa,dd} = [];                    
                end
            end
            
        case 8
            for aa=1:nb_animal_learn
                if aa==1
                    comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 7}')'];
                    comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 7}')'];
                elseif aa==2
                    comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 7}')'];
                    comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 7}')'];
                elseif aa==3
                    comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 8}')'];
                    comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 8}')'];
                elseif aa==4
                    comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 8}')'];
                    comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 8}')'];
                elseif aa==5 %no day 7
                    comb_STC_learn_days.A{aa,dd} = [];
                    comb_STC_learn_days.B{aa,dd} = [];
                elseif aa==6 %no day 7
                    comb_STC_learn_days.A{aa,dd} = [];
                    comb_STC_learn_days.B{aa,dd} = [];                    
                end
            end
        case 9 %day 9
            for aa=1:nb_animal_learn
                if aa==1
                    comb_STC_learn_days.A{aa,dd} = [];
                    comb_STC_learn_days.B{aa,dd} = [];
                elseif aa==2
                    comb_STC_learn_days.A{aa,dd} = [];
                    comb_STC_learn_days.B{aa,dd} = [];
                elseif aa==3
                    comb_STC_learn_days.A{aa,dd} = [];
                    comb_STC_learn_days.B{aa,dd} = [];
                elseif aa==4
                    comb_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 8}')'];
                    comb_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 8}')'];
                elseif aa==5 %no day 7
                    comb_STC_learn_days.A{aa,dd} = [];
                    comb_STC_learn_days.B{aa,dd} = [];
                elseif aa==6 %no day 7
                    comb_STC_learn_days.A{aa,dd} = [];
                    comb_STC_learn_days.B{aa,dd} = [];                    
                end
            end            
    end
end

%% Remove sessions with mismatched/low-quality fields for cross-day - LEARN - animal 3/ session 2/day2 (same)
%animal 3, day/session 2
comb_STC_learn_days.A{3,2} = [];
comb_STC_learn_days.B{3,2} = [];

%% Align by neighboring days LEARN
%all animals treat the same
%D1 vs. D2 - ses 1 vs. 2
%D2 vs. D3 - ses 2 vs. 3
%D6 vs. D7 - ses 4 vs. 5
%D7 vs. D8 - ses 5 vs. 6
%D8 vs. D9 - ses 6 vs. 7

%construct neighboring day equivalent session match for learn (for current
%sets)

%place into function
for dd=[1 2 3 4 5 6 7 8]
    %for each day relative to day 1 choose the corresponding imaging
    %sessions
    %for each day, ...choose sessions
    switch dd
        case 1
            for aa=1:nb_animal_learn
                neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{1, 2}')'];
                neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{1, 2}')'];
            end
        case 2
            for aa=1:nb_animal_learn
                neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{2, 3}')'];
                neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{2, 3}')'];
            end
        case 3 
            for aa=1:nb_animal_learn
                neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{3, 4}')'];
                neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{3, 4}')'];
            end
        case 4
            for aa=1:nb_animal_learn
                if aa==1
                    neighbor_STC_learn_days.A{aa,dd} = [];
                    neighbor_STC_learn_days.B{aa,dd} = [];
                elseif aa==2
                    neighbor_STC_learn_days.A{aa,dd} = [];
                    neighbor_STC_learn_days.B{aa,dd} = [];
                elseif aa==3
                    neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{4, 5}')'];
                    neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{4, 5}')'];
                elseif aa>=4
                    neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{4, 5}')'];
                    neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{4, 5}')'];
                end
            end
        case 5 %day 5 vs day 6
            for aa=1:nb_animal_learn
                if aa==1
                    neighbor_STC_learn_days.A{aa,dd} = [];
                    neighbor_STC_learn_days.B{aa,dd} = [];
                elseif aa==2
                    neighbor_STC_learn_days.A{aa,dd} = [];
                    neighbor_STC_learn_days.B{aa,dd} = [];
                elseif aa==3
                    neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{5, 6}')'];
                    neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{5, 6}')'];
                elseif aa>=4
                     neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{5, 6}')'];
                    neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{5, 6}')'];
                end
            end
            
        case 6
            for aa=1:nb_animal_learn
                if aa==1
                    neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{5, 6}')'];
                    neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{5, 6}')'];

                elseif aa==2
                    neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{5, 6}')'];
                    neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{5, 6}')'];
                elseif aa==3
                    neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{6, 7}')'];
                    neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{6, 7}')'];
                elseif aa==4
                    neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{6, 7}')'];
                    neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{6, 7}')'];
                elseif aa==5
                    neighbor_STC_learn_days.A{aa,dd} = [];
                    neighbor_STC_learn_days.B{aa,dd} = [];    
                elseif aa==6
                    neighbor_STC_learn_days.A{aa,dd} = [];
                    neighbor_STC_learn_days.B{aa,dd} = [];                    
                end
            end
        case 7
            for aa=1:nb_animal_learn
                if aa==1
                    neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{6, 7}')'];
                    neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{6, 7}')'];
                elseif aa==2
                    neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{6, 7}')'];
                    neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{6, 7}')'];
                elseif aa==3
                    neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{7, 8}')'];
                    neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{7, 8}')'];
                elseif aa==4
                    neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{7, 8}')'];
                    neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{7, 8}')'];
                elseif aa==5
                    neighbor_STC_learn_days.A{aa,dd} = [];
                    neighbor_STC_learn_days.B{aa,dd} = [];    
                elseif aa==6
                    neighbor_STC_learn_days.A{aa,dd} = [];
                    neighbor_STC_learn_days.B{aa,dd} = [];                    
                end
            end
        case 8 %8 vs day 9
            for aa=1:nb_animal_learn
                if aa==1
                    neighbor_STC_learn_days.A{aa,dd} = [];
                    neighbor_STC_learn_days.B{aa,dd} = []; 
                elseif aa==2
                    neighbor_STC_learn_days.A{aa,dd} = [];
                    neighbor_STC_learn_days.B{aa,dd} = []; 
                elseif aa==3
                    neighbor_STC_learn_days.A{aa,dd} = [];
                    neighbor_STC_learn_days.B{aa,dd} = []; 
                elseif aa==4
                    neighbor_STC_learn_days.A{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.A{8, 9}')'];
                    neighbor_STC_learn_days.B{aa,dd} = [cell2mat(TC_corr_match_learning{aa}.tc_corr_match.matching_ROI_all_day_STC.ts.B{8, 9}')'];
                elseif aa==5
                    neighbor_STC_learn_days.A{aa,dd} = [];
                    neighbor_STC_learn_days.B{aa,dd} = [];    
                elseif aa==6
                    neighbor_STC_learn_days.A{aa,dd} = [];
                    neighbor_STC_learn_days.B{aa,dd} = [];                    
                end
            end             
    end
end

%% Make correction for misaligned day 2 session for animal 3
%1 vs. 2
neighbor_STC_learn_days.A{3,1} = [];
%2 vs. 3
neighbor_STC_learn_days.A{3,2} = [];

%% Collapse day-matched STCs for TS tuned neurons across days neighbor - LEARN

%for each day relative to dd index (dd=1 --> 1 vs. 2; dd=2 --> 
for dd=1:8 %1 day short of max for neighbor
    %A trials
    neighbor_STC_learn_days_all.A{dd} = cell2mat(neighbor_STC_learn_days.A(:,dd));
    %B trials
    neighbor_STC_learn_days_all.B{dd} = cell2mat(neighbor_STC_learn_days.B(:,dd));
end


%% Get mean and sem for TC on neighboring LEARNING days (uses names defined for learn - works the same)
%across all pooled neurons
[mean_TC_learn_neighbor_all,sem_TC_learn_neighbor_all] = TC_corr_learn_neighbor_all_neurons(neighbor_STC_learn_days_all);

%do per animal (uncollapsed input - need separate out from each animal) -
[mean_mean_PV_learn_neighbor_all,sem_PV_learn_neighbor_all] = PV_corr_learn_neighbor_all_neurons(neighbor_STC_learn_days);

%% Collapse day-matched STCs for TS tuned neurons across days (rel D1) - RECALL

%for each day relative to 1 (1 vs. 2. 1 vs. 3...)
for dd=2:6
    %A trials
    comb_STC_recall_days_all.A{dd} = cell2mat(comb_STC_recall_days.A(:,dd));
    %B trials
    comb_STC_recall_days_all.B{dd} = cell2mat(comb_STC_recall_days.B(:,dd));
end

%% Calculate TC correlation relative to day 1

%correlation coefcient for each neuron
for dd=2:6
    %A
    TC_comb_STC_recall_days_all.A{dd} = diag(corr(comb_STC_recall_days_all.A{dd}(:,1:100)',comb_STC_recall_days_all.A{dd}(:,101:200)'));
    %B
    TC_comb_STC_recall_days_all.B{dd} = diag(corr(comb_STC_recall_days_all.B{dd}(:,1:100)',comb_STC_recall_days_all.B{dd}(:,101:200)'));
end

%get mean for each day (across all animals)
%A
mean_TC_recall_all.A = cellfun(@mean,TC_comb_STC_recall_days_all.A);
%B
mean_TC_recall_all.B = cellfun(@mean,TC_comb_STC_recall_days_all.B);

%% Get sem

sem_TC_recall_all.A = [nan, cellfun(@(x) std(x,0,1), TC_comb_STC_recall_days_all.A(2:6),'UniformOutput',true)./...
                        sqrt(cellfun(@(x) size(x,1), TC_comb_STC_recall_days_all.A(2:6),'UniformOutput',true))];
                    
sem_TC_recall_all.B = [nan, cellfun(@(x) std(x,0,1), TC_comb_STC_recall_days_all.B(2:6),'UniformOutput',true)./...
                        sqrt(cellfun(@(x) size(x,1), TC_comb_STC_recall_days_all.B(2:6),'UniformOutput',true))];
                    
                    
%% Collapse day-matched STCs for TS tuned neurons across days neighbor - RECALL

%for each day relative to dd index (dd=1 --> 1 vs. 2; dd=2 --> 
for dd=1:5
    %A trials
    neighbor_STC_recall_days_all.A{dd} = cell2mat(neighbor_STC_recall_days.A(:,dd));
    %B trials
    neighbor_STC_recall_days_all.B{dd} = cell2mat(neighbor_STC_recall_days.B(:,dd));
end

%% Get mean and sem for TC on neighboring RECALL days
%1 less for long-term recall sessions 
nb_ses = 5;
[mean_TC_recall_neighbor_all,sem_TC_recall_neighbor_all] = TC_corr_recall_neighbor_all_neurons(neighbor_STC_recall_days_all,nb_ses);

%do per animal (uncollapsed input - need separate out from each animal) -
[mean_mean_PV_recall_neighbor_all,sem_PV_recall_neighbor_all] = PV_corr_recall_neighbor_all_neurons(neighbor_STC_recall_days,nb_ses);
                  
                    
%% Split TC and PV calculation by ANIMAL (rather than neuron pool) - seems necessary for PV - RECALL

%calculate PV/TC correlation for each relative session for each animal
for aa=1:nb_animal_recall
    for dd=2:6
        %A
        %get diagonal of correlation
        TC_recall_days_ind.A{aa,dd} = diag(corr(comb_STC_recall_days.A{aa,dd}(:,1:100)',comb_STC_recall_days.A{aa,dd}(:,101:200)'));
        PV_recall_days_ind.A{aa,dd} = diag(corr(comb_STC_recall_days.A{aa,dd}(:,1:100),comb_STC_recall_days.A{aa,dd}(:,101:200)));
        
        %get mean TC score for each animal/day
        mean_TC_recall_days_ind.A(aa,dd) = nanmean(TC_recall_days_ind.A{aa,dd});
        mean_PV_recall_days_ind.A(aa,dd) = nanmean(PV_recall_days_ind.A{aa,dd});
        
        %B
        %get diagonal of correlation
        TC_recall_days_ind.B{aa,dd} = diag(corr(comb_STC_recall_days.B{aa,dd}(:,1:100)',comb_STC_recall_days.B{aa,dd}(:,101:200)'));
        PV_recall_days_ind.B{aa,dd} = diag(corr(comb_STC_recall_days.B{aa,dd}(:,1:100),comb_STC_recall_days.B{aa,dd}(:,101:200)));
        
        %get mean TC score for each animal/day
        mean_TC_recall_days_ind.B(aa,dd) = nanmean(TC_recall_days_ind.B{aa,dd});
        mean_PV_recall_days_ind.B(aa,dd) = nanmean(PV_recall_days_ind.B{aa,dd});
    end
end

%use same number of animals for here b/c have same number of animal on each
%session
for dd=2:6
    %A
    %get sem TC score for each animal/day
    sem_TC_recall_days_ind.A(:,dd) = std(mean_TC_recall_days_ind.A(:,dd))./sqrt(nb_animal_recall);
    sem_PV_recall_days_ind.A(:,dd) = nanstd(mean_PV_recall_days_ind.A(:,dd))./sqrt(nb_animal_recall);
    %B
    sem_TC_recall_days_ind.B(:,dd) = std(mean_TC_recall_days_ind.B(:,dd))./sqrt(nb_animal_recall);
    sem_PV_recall_days_ind.B(:,dd) = nanstd(mean_PV_recall_days_ind.B(:,dd))./sqrt(nb_animal_recall);
    %A
    %get sem TC score for each animal/day
    mean_mean_TC_recall_days_ind.A(:,dd) = mean(mean_TC_recall_days_ind.A(:,dd));
    mean_mean_PV_recall_days_ind.A(:,dd) = mean(mean_PV_recall_days_ind.A(:,dd));
    %B
    mean_mean_TC_recall_days_ind.B(:,dd) = mean(mean_TC_recall_days_ind.B(:,dd));
    mean_mean_PV_recall_days_ind.B(:,dd) = mean(mean_PV_recall_days_ind.B(:,dd));
end

%% Split TC and PV calculation by animal (rather than neuron pool) - seems necessary for PV - LEARN

%calculate PV/TC correlation for each relative session for each animal
for aa=1:nb_animal_learn
    for dd=2:9
        %A
        if ~isempty(comb_STC_learn_days.A{aa,dd})
            %get diagonal of correlation
            TC_learn_days_ind.A{aa,dd} = diag(corr(comb_STC_learn_days.A{aa,dd}(:,1:100)',comb_STC_learn_days.A{aa,dd}(:,101:200)'));
            PV_learn_days_ind.A{aa,dd} = diag(corr(comb_STC_learn_days.A{aa,dd}(:,1:100),comb_STC_learn_days.A{aa,dd}(:,101:200)));
            
            %get mean TC score for each animal/day
            mean_TC_learn_days_ind.A(aa,dd) = nanmean(TC_learn_days_ind.A{aa,dd});
            mean_PV_learn_days_ind.A(aa,dd) = nanmean(PV_learn_days_ind.A{aa,dd});
            
        else %set to empty
            TC_learn_days_ind.A{aa,dd} = [];
            PV_learn_days_ind.A{aa,dd} = [];
            mean_TC_learn_days_ind.A(aa,dd) = nan;
            mean_PV_learn_days_ind.A(aa,dd) = nan;
        end
            
        %B
        if ~isempty(comb_STC_learn_days.B{aa,dd})
        %get diagonal of correlation
        TC_learn_days_ind.B{aa,dd} = diag(corr(comb_STC_learn_days.B{aa,dd}(:,1:100)',comb_STC_learn_days.B{aa,dd}(:,101:200)'));
        PV_learn_days_ind.B{aa,dd} = diag(corr(comb_STC_learn_days.B{aa,dd}(:,1:100),comb_STC_learn_days.B{aa,dd}(:,101:200)));
        
        %get mean TC score for each animal/day
        mean_TC_learn_days_ind.B(aa,dd) = nanmean(TC_learn_days_ind.B{aa,dd});
        mean_PV_learn_days_ind.B(aa,dd) = nanmean(PV_learn_days_ind.B{aa,dd});
        else
            TC_learn_days_ind.B{aa,dd} = [];
            PV_learn_days_ind.B{aa,dd} = [];
            mean_TC_learn_days_ind.B(aa,dd) = nan;
            mean_PV_learn_days_ind.B(aa,dd) = nan;
            
        end
    end
end

%number of animals adjusted for by day (can use A only)
nb_animals_per_day_adj_learn = sum(~cellfun(@isempty,comb_STC_learn_days.A),1);

for dd=2:9
    %A
    %get sem TC score for each animal/day
    sem_TC_learn_days_ind.A(:,dd) = nanstd(mean_TC_learn_days_ind.A(:,dd))./sqrt(nb_animals_per_day_adj_learn(dd));
    sem_PV_learn_days_ind.A(:,dd) = nanstd(mean_PV_learn_days_ind.A(:,dd))./sqrt(nb_animals_per_day_adj_learn(dd));
    %B
    sem_TC_learn_days_ind.B(:,dd) = nanstd(mean_TC_learn_days_ind.B(:,dd))./sqrt(nb_animals_per_day_adj_learn(dd));
    sem_PV_learn_days_ind.B(:,dd) = nanstd(mean_PV_learn_days_ind.B(:,dd))./sqrt(nb_animals_per_day_adj_learn(dd));
    %A
    %get sem TC score for each animal/day
    mean_mean_TC_learn_days_ind.A(:,dd) = nanmean(mean_TC_learn_days_ind.A(:,dd));
    mean_mean_PV_learn_days_ind.A(:,dd) = nanmean(mean_PV_learn_days_ind.A(:,dd));
    %B
    mean_mean_TC_learn_days_ind.B(:,dd) = nanmean(mean_TC_learn_days_ind.B(:,dd));
    mean_mean_PV_learn_days_ind.B(:,dd) = nanmean(mean_PV_learn_days_ind.B(:,dd));
end



%% Collapse day-matched STCs for TS tuned neurons across days (rel D1)

%for each day relative to 1 (1 vs. 2. 1 vs. 3...)
for dd=2:9
    %A trials
    comb_STC_learn_days_all.A{dd} = cell2mat(comb_STC_learn_days.A(:,dd));
    %B trials
    comb_STC_learn_days_all.B{dd} = cell2mat(comb_STC_learn_days.B(:,dd));
end

%% Calculate TC/PV correlation relative to day 1

%correlation coefcient for each neuron
for dd=2:9
    %A
    TC_comb_STC_learn_days_all.A{dd} = diag(corr(comb_STC_learn_days_all.A{dd}(:,1:100)',comb_STC_learn_days_all.A{dd}(:,101:200)'));
    %B
    TC_comb_STC_learn_days_all.B{dd} = diag(corr(comb_STC_learn_days_all.B{dd}(:,1:100)',comb_STC_learn_days_all.B{dd}(:,101:200)'));
    %A
    PV_comb_STC_learn_days_all.A{dd} = diag(corr(comb_STC_learn_days_all.A{dd}(:,1:100),comb_STC_learn_days_all.A{dd}(:,101:200)));
    %B
    PV_comb_STC_learn_days_all.B{dd} = diag(corr(comb_STC_learn_days_all.B{dd}(:,1:100),comb_STC_learn_days_all.B{dd}(:,101:200)));
end

%get mean for each day (across all animals)
%A
mean_TC_learn_all.A = cellfun(@mean,TC_comb_STC_learn_days_all.A);
%B
mean_TC_learn_all.B = cellfun(@mean,TC_comb_STC_learn_days_all.B);

%PV correlation
%A
mean_PV_learn_all.A = cellfun(@mean,PV_comb_STC_learn_days_all.A);
%B
mean_PV_learn_all.B = cellfun(@mean,PV_comb_STC_learn_days_all.B);

%% Get sem
%TC
sem_TC_learn_all.A = [nan, cellfun(@(x) std(x,0,1), TC_comb_STC_learn_days_all.A(2:9),'UniformOutput',true)./...
    sqrt(cellfun(@(x) size(x,1), TC_comb_STC_learn_days_all.A(2:9),'UniformOutput',true))];

sem_TC_learn_all.B = [nan, cellfun(@(x) std(x,0,1), TC_comb_STC_learn_days_all.B(2:9),'UniformOutput',true)./...
    sqrt(cellfun(@(x) size(x,1), TC_comb_STC_learn_days_all.B(2:9),'UniformOutput',true))];
%PV
sem_PV_learn_all.A = [nan, cellfun(@(x) nanstd(x,0,1), PV_comb_STC_learn_days_all.A(2:9),'UniformOutput',true)./...
    sqrt(cellfun(@(x) size(x,1), PV_comb_STC_learn_days_all.A(2:9),'UniformOutput',true))];

sem_PV_learn_all.B = [nan, cellfun(@(x) nanstd(x,0,1), PV_comb_STC_learn_days_all.B(2:9),'UniformOutput',true)./...
    sqrt(cellfun(@(x) size(x,1), PV_comb_STC_learn_days_all.B(2:9),'UniformOutput',true))];
                    
%% Plot line plot with errorbars - all neurons - T.C.

figure
hold on 
title('TC correlation - T.S.')
ylim([0 1])
xlim([1 7])
xticks(2:6)
xticklabels({'6','16','20','25','30'})
yticks(0:0.2:1)
xlabel('Days')
ylabel('Correlation score')
%Recall
rA = errorbar(2:6,mean_TC_recall_all.A(2:6),sem_TC_recall_all.A(2:6),'LineStyle','-','Linewidth',2,'Color',[65,105,225]/255);
rB = errorbar(2:6,mean_TC_recall_all.B(2:6),sem_TC_recall_all.B(2:6),'LineStyle','-','Linewidth',2,'Color',[220,20,60]/255);
%Learn
%lA = errorbar(2:9,mean_TC_learn_all.A(2:9),sem_TC_learn_all.A(2:9),'LineStyle','--','Linewidth',2,'Color',[65,105,225]/255);
%lB = errorbar(2:9,mean_TC_learn_all.B(2:9),sem_TC_learn_all.B(2:9),'LineStyle','--','Linewidth',2,'Color',[220,20,60]/255);
set(gca,'FontSize',16)
set(gca,'Linewidth',2)

%legend([lA,lB,rA,rB],{'Learning A','Learning B','Recall A','Recall B'},'location','northeast')

%% Plot line plot with errorbars - by animal - T.C.
figure
hold on 
title('TC correlation - T.S.')
ylim([0 1])
xlim([1 7])
xticks(2:6)
xticklabels({'6','16','20','25','30'})
yticks(0:0.2:1)
xlabel('Days')
ylabel('Correlation score')
%Recall
rA = errorbar(2:6,mean_mean_TC_recall_days_ind.A(2:6),sem_TC_recall_days_ind.A(2:6),'LineStyle','-','Linewidth',2,'Color',[65,105,225]/255);
rB = errorbar(2:6,mean_mean_TC_recall_days_ind.B(2:6),sem_TC_recall_days_ind.B(2:6),'LineStyle','-','Linewidth',2,'Color',[220,20,60]/255);
%Learn
% lA = errorbar(2:9,mean_mean_TC_learn_days_ind.A(2:9),sem_TC_learn_days_ind.A(2:9),'LineStyle','--','Linewidth',2,'Color',[65,105,225]/255);
% rB = errorbar(2:9,mean_mean_TC_learn_days_ind.B(2:9),sem_TC_learn_days_ind.B(2:9),'LineStyle','--','Linewidth',2,'Color',[220,20,60]/255);

set(gca,'FontSize',16)
set(gca,'Linewidth',2)

%legend([lA,lB,rA,rB],{'Learning A','Learning B','Recall A','Recall B'},'location','northeast')

%% Plot line plot with errorbars - by animal - P.V.
figure
hold on 
title('PV correlation - T.S.')
ylim([0 1])
xlim([1 7])
xticks(2:6)
xticklabels({'6','16','20','25','30'})
yticks(0:0.2:1)
xlabel('Days')
ylabel('Correlation score')
%Recall
rA = errorbar(2:6,mean_mean_PV_recall_days_ind.A(2:6),sem_PV_recall_days_ind.A(2:6),'LineStyle','-','Linewidth',2,'Color',[65,105,225]/255);
rB = errorbar(2:6,mean_mean_PV_recall_days_ind.B(2:6),sem_PV_recall_days_ind.B(2:6),'LineStyle','-','Linewidth',2,'Color',[220,20,60]/255);
%Learn
%lA = errorbar(2:9,mean_mean_PV_learn_days_ind.A(2:9),sem_PV_learn_days_ind.A(2:9),'LineStyle','--','Linewidth',2,'Color',[65,105,225]/255);
%lB = errorbar(2:9,mean_mean_PV_learn_days_ind.B(2:9),sem_PV_learn_days_ind.B(2:9),'LineStyle','--','Linewidth',2,'Color',[220,20,60]/255);

set(gca,'FontSize',16)
set(gca,'Linewidth',2)
%legend([lA,lB,rA,rB],{'Learning A','Learning B','Recall A','Recall B'},'location','northeast')

%% Plot NEIGHBOR T.C correlation for learning/recall - individual neurons pooled

figure
hold on 
title('TC correlation - T.S. - neighbors')
ylim([0 1])
xlim([0 6])
xticks(1:5)
xticklabels({'1-6','6-16','16-20','20-25','25-30'})
%yticks(0:0.2:1)
xlabel('Neighboring days')
ylabel('Correlation score')
%Recall
rA = errorbar(1:5,mean_TC_recall_neighbor_all.A,sem_TC_recall_neighbor_all.A,'LineStyle','-','Linewidth',2,'Color',[65,105,225]/255);
rB = errorbar(1:5,mean_TC_recall_neighbor_all.B,sem_TC_recall_neighbor_all.B,'LineStyle','-','Linewidth',2,'Color',[220,20,60]/255);
%Learn
%lA = errorbar(1:8,mean_TC_learn_neighbor_all.A,sem_TC_learn_neighbor_all.A,'LineStyle','--','Linewidth',2,'Color',[65,105,225]/255);
%lB = errorbar(1:8,mean_TC_learn_neighbor_all.B,sem_TC_learn_neighbor_all.B,'LineStyle','--','Linewidth',2,'Color',[220,20,60]/255);
set(gca,'FontSize',16)
set(gca,'Linewidth',2)

%legend([lA,lB,rA,rB],{'Learning A','Learning B','Recall A','Recall B'},'location','southeast')

%% Plot NEIGHBOR P.V. correlation for learning/recall - by animals

figure
hold on 
title('PV correlation - T.S. - neighbors')
ylim([0 1])
xlim([0 6])
xticks(1:5)
xticklabels({'1-6','6-16','16-20','20-25','25-30'})
%yticks(0:0.2:1)
xlabel('Neighboring days')
ylabel('Correlation score')
%Recall
rA = errorbar(1:5,mean_mean_PV_recall_neighbor_all.A,sem_PV_recall_neighbor_all.A,'LineStyle','-','Linewidth',2,'Color',[65,105,225]/255);
rB = errorbar(1:5,mean_mean_PV_recall_neighbor_all.B,sem_PV_recall_neighbor_all.B,'LineStyle','-','Linewidth',2,'Color',[220,20,60]/255);
%Learn
%lA = errorbar(1:8,mean_mean_PV_learn_neighbor_all.A,sem_PV_learn_neighbor_all.A,'LineStyle','--','Linewidth',2,'Color',[65,105,225]/255);
%lB = errorbar(1:8,mean_mean_PV_learn_neighbor_all.B,sem_PV_learn_neighbor_all.B,'LineStyle','--','Linewidth',2,'Color',[220,20,60]/255);
set(gca,'FontSize',16)
set(gca,'Linewidth',2)

%legend([lA,lB,rA,rB],{'Learning A','Learning B','Recall A','Recall B'},'location','southeast')


%% Plots - OLD CODE BELOW
%{
%Plot TS tuned matching place cells for A and B trials
figure
for ss=2:7
    subplot(1,6,ss-1)
    imagesc(comb_STC_recall_sort.A{1,ss})
    hold on
    caxis([0 1])
    colormap('jet')
end

figure
for ss=2:7
    subplot(1,6,ss-1)
    imagesc(comb_STC_recall_sort.B{1,ss})
    hold on
    caxis([0 1])
    colormap('jet')
end

comb_STC_recall_sort_7 = comb_STC_recall{1,7}(sortOrder,:);


diag(corr(comb_STC_recall_sort_7(:,1:100)',comb_STC_recall_sort_7(:,101:200)'))

mean(diag(corr(comb_STC_recall_sort_3(:,1:100)',comb_STC_recall_sort_3(:,101:200)')))

mean(diag(corr(comb_STC_learning_sort_3(:,1:100)',comb_STC_learning_sort_3(:,101:200)')))

diag(corr(comb_STC_learning_sort_6(:,1:100)',comb_STC_learning_sort_6(:,101:200)'))

[~,maxBin] = max(comb_STC_recall{1,3}(:,1:100)', [], 1,'includenan');
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,sortOrder] = sort(maxBin,'ascend');

comb_STC_recall_sort_3 = comb_STC_recall{1,3}(sortOrder,:);

%% PV correlate the matching maps

%find pv diagnonals for learning A
for ss=1:6
    PV_diag.learn.A(ss,:) = diag(corr(comb_STC_learning_sort.A{1,ss}(:,1:100),comb_STC_learning_sort.A{1,ss}(:,101:200)))
    PV_diag.learn.B(ss,:) = diag(corr(comb_STC_learning_sort.B{1,ss}(:,1:100),comb_STC_learning_sort.B{1,ss}(:,101:200)))
end

for ss=1:7
    PV_diag.recall.A(ss,:) = diag(corr(comb_STC_recall_sort.A{1,ss}(:,1:100),comb_STC_recall_sort.A{1,ss}(:,101:200)))
    PV_diag.recall.B(ss,:) = diag(corr(comb_STC_recall_sort.B{1,ss}(:,1:100),comb_STC_recall_sort.B{1,ss}(:,101:200)))
end

%% TC correlate the matching maps
for ss=1:6
    TC_diag.learn.A{ss} = diag(corr(comb_STC_learning_sort.A{1,ss}(:,1:100)',comb_STC_learning_sort.A{1,ss}(:,101:200)'))
    TC_diag.learn.B{ss} = diag(corr(comb_STC_learning_sort.B{1,ss}(:,1:100)',comb_STC_learning_sort.B{1,ss}(:,101:200)'))
end

for ss=1:7
    TC_diag.recall.A{ss} = diag(corr(comb_STC_recall_sort.A{1,ss}(:,1:100)',comb_STC_recall_sort.A{1,ss}(:,101:200)'))
    TC_diag.recall.B{ss} = diag(corr(comb_STC_recall_sort.B{1,ss}(:,1:100)',comb_STC_recall_sort.B{1,ss}(:,101:200)'))
end


%learning
%pick day
day_sel =6

figure
subplot(1,2,1)
imagesc(corr(comb_STC_learning_sort.A{1,day_sel}(:,1:100),comb_STC_learning_sort.A{1,day_sel}(:,101:200)))
hold on
title('A')
caxis([0 1])
colormap('jet')

subplot(1,2,2)
imagesc(corr(comb_STC_learning_sort.B{1,day_sel}(:,1:100),comb_STC_learning_sort.B{1,day_sel}(:,101:200)))
hold on
title('B')
caxis([0 1])
colormap('jet')
colorbar


%recall
day_sel =2

figure
subplot(1,2,1)
imagesc(corr(comb_STC_recall_sort.A{1,day_sel}(:,1:100),comb_STC_recall_sort.A{1,day_sel}(:,101:200)))
hold on
title('A')
caxis([0 1])
colormap('jet')

subplot(1,2,2)
imagesc(corr(comb_STC_recall_sort.B{1,day_sel}(:,1:100),comb_STC_recall_sort.B{1,day_sel}(:,101:200)))
hold on
title('B')
caxis([0 1])
colormap('jet')
colorbar

corr(comb_STC_recall_sort.A{1,7}(:,1:100),comb_STC_recall_sort.A{1,7}(:,101:200))

%% Extract neurons with high TC correlation scores (>0.6)
ss=6
%high_TC_idx.learn.A{ss} = find(TC_diag.learn.A{1, ss} < 0.3);
%high_TC_idx.learn.A{ss} = find(TC_diag.learn.A{1, ss} < 0.6 & TC_diag.learn.A{1, ss} >= 0.3);
high_TC_idx.learn.A{ss} = find(TC_diag.learn.A{1, ss} >= 0.6);

figure
imagesc(comb_STC_learning_sort.A{1, ss}(high_TC_idx.learn.A{ss},:))
hold on
caxis([0 1])
colormap('jet')

ss=6
%high_TC_idx.learn.B{ss} = find(TC_diag.learn.B{1, ss} < 0.3);
%high_TC_idx.learn.B{ss} = find(TC_diag.learn.B{1, ss} < 0.6 & TC_diag.learn.B{1, ss} >= 0.3);
high_TC_idx.learn.B{ss} = find(TC_diag.learn.B{1, ss} >= 0.6);

figure
imagesc(comb_STC_learning_sort.B{1, ss}(high_TC_idx.learn.B{ss},:))
hold on
caxis([0 1])
colormap('jet')

% ss=6
% %high_TC_idx.recall.A{ss} = find(TC_diag.recall.A{1, ss} < 0.3);
% %high_TC_idx.recall.A{ss} = find(TC_diag.recall.A{1, ss} < 0.6 & TC_diag.recall.A{1, ss} >= 0.3);
% high_TC_idx.recall.A{ss} = find(TC_diag.recall.A{1, ss} >= 0.6);
% 
% figure
% imagesc(comb_STC_recall_sort.A{1, ss}(high_TC_idx.recall.A{ss},:))
% hold on
% caxis([0 1])
% colormap('jet')

ss=6
%high_TC_idx.recall.B{ss} = find(TC_diag.recall.B{1, ss} < 0.3);
%high_TC_idx.recall.B{ss} = find(TC_diag.recall.B{1, ss} < 0.6 & TC_diag.recall.B{1, ss} >= 0.3);
high_TC_idx.recall.B{ss} = find(TC_diag.recall.B{1, ss} >= 0.6);

figure
imagesc(comb_STC_recall_sort.B{1, ss}(high_TC_idx.recall.B{ss},:))
hold on
caxis([0 1])
colormap('jet')
%}
end

