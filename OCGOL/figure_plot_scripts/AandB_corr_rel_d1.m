function [outputArg1,outputArg2] = AandB_corr_rel_d1(TC_corr_match_learning,TC_corr_match_recall)

%% Get number of animals

nb_recall_animals = size(TC_corr_match_recall,2);
nb_learn_animals = size(TC_corr_match_learning,2);

%% Get number of sessions for each animal

%recall
for aa=1:nb_recall_animals
    nb_ses_recall(aa) = size(TC_corr_match_recall{aa}.tc_corr_match.STC_mat_AB_A,2);
end

%learn
for aa=1:nb_learn_animals
    nb_ses_learn(aa) = size(TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_A,2);
end

%% Relative to D1 comparison - LEARN - relative to session only (not sorted by day yet or low trial count/field shift filtered)

%preallocate with nans
TC_corr_AB_matching_learning_early = nan(nb_learn_animals,max(nb_ses_learn));
TC_corr_AB_matching_learning_late = nan(nb_learn_animals,max(nb_ses_learn));

%for all learning animals
for aa=1:nb_learn_animals
    
    %session relative to session 1 (d1)
    %for all sessions for given animal
    for day_sel=2:nb_ses_learn(aa)
        
        %for A&B tuned neurons, get STC (100 bin) on A trials and STC (100
        %bin) on B trials
        %session 1 matches
        STC_AB_1 = [TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_A{1, day_sel}(:,1:100), TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_B{1, day_sel}(:,1:100)];
        %comparison session matched (to session 1 matched neurons)
        STC_AB_later = [TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_A{1, day_sel}(:,101:200), TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_B{1, day_sel}(:,101:200)];
        
        %% Sort D1 - sort all the STCs by session 1, A trials
        %A sort
        [~,test_maxBin] = max(STC_AB_1(:,1:100)', [], 1,'includenan');
        %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
        [~,test_sortOrder] = sort(test_maxBin,'ascend');
        
        STC_AB_sort_A = STC_AB_1(test_sortOrder,1:100);
        
        %B
        % [~,test_maxBin] = max(STC_AB_1(:,101:200)', [], 1,'includenan');
        % %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
        % [~,test_sortOrder] = sort(test_maxBin,'ascend');
        
        STC_AB_sort_B = STC_AB_1(test_sortOrder,101:200);
        
        %% Sort later day matching
        %A sort
        % [~,test_maxBin] = max(STC_AB_later(:,1:100)', [], 1,'includenan');
        % %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
        % [~,test_sortOrder] = sort(test_maxBin,'ascend');
        STC_AB_later_sort_A = STC_AB_later(test_sortOrder,1:100);
        
        %B
        % [~,test_maxBin] = max(STC_AB_later(:,101:200)', [], 1,'includenan');
        % %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
        % [~,test_sortOrder] = sort(test_maxBin,'ascend');
        STC_AB_later_sort_B = STC_AB_later(test_sortOrder,101:200);
        
        %% Correlate tuning curves each session against each other on session 1 vs later session
        %if matrix not empty - mean correlation score (across all neurons)
        if ~isempty(STC_AB_sort_A)
            TC_corr_AB_matching_learning_early(aa,day_sel) = mean(diag(corr(STC_AB_sort_A', STC_AB_sort_B','Type','Pearson','rows','all')));
            TC_corr_AB_matching_learning_late(aa,day_sel) = mean(diag(corr(STC_AB_later_sort_A', STC_AB_later_sort_B','Type','Pearson','rows','all')));
        else
            TC_corr_AB_matching_learning_early(aa,day_sel) = nan;
            TC_corr_AB_matching_learning_late(aa,day_sel) = nan;
            
        end
        %get diagonal correlation values (correlation value for each
        %neuron)
        if ~isempty(STC_AB_sort_A)
            TC_corr_AB_matching_learning_early_diag{aa,day_sel} = diag(corr(STC_AB_sort_A', STC_AB_sort_B','Type','Pearson','rows','all'));
            TC_corr_AB_matching_learning_late_diag{aa,day_sel} = diag(corr(STC_AB_later_sort_A', STC_AB_later_sort_B','Type','Pearson','rows','all'));
        else
            TC_corr_AB_matching_learning_early_diag{aa,day_sel} = nan;
            TC_corr_AB_matching_learning_late_diag{aa,day_sel} = nan;
            
        end
        
        %clear the starting matrices
        clear STC_AB_1 STC_AB_later
        
    end
end

%% Filter sessions relative to session 1 based on field shift/trial count - LEARN
%animal 1, session 3,4 - A/B both trials - means
TC_corr_AB_matching_learning_early(1,[3,4])= nan;
TC_corr_AB_matching_learning_late(1,[3,4])= nan;
%each neuron correlation
TC_corr_AB_matching_learning_early_diag(1,[3,4]) = {nan};
TC_corr_AB_matching_learning_late_diag(1,[3,4]) = {nan};

%animal 3, session 2 - A/B both trials - means
TC_corr_AB_matching_learning_early(3,[2])= nan;
TC_corr_AB_matching_learning_late(3,[2])= nan;
%each neuron correlation
TC_corr_AB_matching_learning_early_diag(3,[2]) = {nan};
TC_corr_AB_matching_learning_late_diag(3,[2]) = {nan};

%animal 5, session 6 and 7 - A/B both trials - means
TC_corr_AB_matching_learning_early(5,[6,7])= nan;
TC_corr_AB_matching_learning_late(5,[6,7])= nan;
%each neuron correlation
TC_corr_AB_matching_learning_early_diag(5,[6,7]) = {nan};
TC_corr_AB_matching_learning_late_diag(5,[6,7]) = {nan};



%% Cross-comparison to all other sessions - LEARN - DEBUG THIS TO WORK ACROSS ALL ANIMALS/SESSION
for aa=1:nb_learn_animals
    
    %session relative to session 1st
    for day_sel = 1:6
        
        %session relative to 2nd session - CONTINUE HERE
        for day_sel_2 = 1:6
            
            STC_AB_1 = [TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_A{day_sel, day_sel_2}(:,1:100), TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_B{day_sel, day_sel_2}(:,1:100)];
            STC_AB_later = [TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_A{day_sel, day_sel_2}(:,101:200), TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_B{day_sel, day_sel_2}(:,101:200)];
            
            %% Sort D1
            %A sort
            [~,test_maxBin] = max(STC_AB_1(:,1:100)', [], 1,'includenan');
            %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
            [~,test_sortOrder] = sort(test_maxBin,'ascend');
            
            STC_AB_sort_A = STC_AB_1(test_sortOrder,1:100);
            
            %B
            % [~,test_maxBin] = max(STC_AB_1(:,101:200)', [], 1,'includenan');
            % %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
            % [~,test_sortOrder] = sort(test_maxBin,'ascend');
            
            STC_AB_sort_B = STC_AB_1(test_sortOrder,101:200);
            
            
            %% Sort later day matching
            %A sort
            % [~,test_maxBin] = max(STC_AB_later(:,1:100)', [], 1,'includenan');
            % %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
            % [~,test_sortOrder] = sort(test_maxBin,'ascend');
            
            STC_AB_later_sort_A = STC_AB_later(test_sortOrder,1:100);
            
            %B
            % [~,test_maxBin] = max(STC_AB_later(:,101:200)', [], 1,'includenan');
            % %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
            % [~,test_sortOrder] = sort(test_maxBin,'ascend');
            
            STC_AB_later_sort_B = STC_AB_later(test_sortOrder,101:200);
            
            %% Correlate each session against each other
            %if nay matrix empty
%             if ~isempty(STC_AB_sort_A)
%                 TC_corr_AB_matching_learning_early(aa,day_sel) = mean(diag(corr(STC_AB_sort_A', STC_AB_sort_B','Type','Pearson','rows','all')));
%                 TC_corr_AB_matching_learning_late(aa,day_sel) = mean(diag(corr(STC_AB_later_sort_A', STC_AB_later_sort_B','Type','Pearson','rows','all')));
%             else
%                 TC_corr_AB_matching_learning_early(aa,day_sel) = nan;
%                 TC_corr_AB_matching_learning_late(aa,day_sel) = nan;
%                 
%             end
            %get diagonal correlation values
            if ~isempty(STC_AB_sort_A)
                TC_corr_AB_matching_learning_early_diag_all{aa,day_sel, day_sel_2}= diag(corr(STC_AB_sort_A', STC_AB_sort_B','Type','Pearson','rows','all'));
                TC_corr_AB_matching_learning_late_diag_all{aa,day_sel, day_sel_2} = diag(corr(STC_AB_later_sort_A', STC_AB_later_sort_B','Type','Pearson','rows','all'));
            else
                TC_corr_AB_matching_learning_early_diag_all{aa,day_sel, day_sel_2} = nan;
                TC_corr_AB_matching_learning_late_diag_all{aa,day_sel, day_sel_2} = nan;
                
            end
            
        end
    end
end


%% Relative to D1 comparison - recall

%preallocate with nans
TC_corr_AB_matching_recall_early = nan(nb_recall_animals,max(nb_ses_recall));
TC_corr_AB_matching_recall_late = nan(nb_recall_animals,max(nb_ses_recall));

for aa=1:nb_recall_animals
    
    for day_sel=2:7
        
%         STC_AB_1 = [TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_A{1, day_sel}(:,1:100), TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_B{1, day_sel}(:,1:100)];
%         STC_AB_later = [TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_A{1, day_sel}(:,101:200), TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_B{1, day_sel}(:,101:200)];

        STC_AB_1 = [TC_corr_match_recall{aa}.tc_corr_match.STC_mat_AB_A{1, day_sel}(:,1:100), TC_corr_match_recall{aa}.tc_corr_match.STC_mat_AB_B{1, day_sel}(:,1:100)];
        STC_AB_later = [TC_corr_match_recall{aa}.tc_corr_match.STC_mat_AB_A{1, day_sel}(:,101:200), TC_corr_match_recall{aa}.tc_corr_match.STC_mat_AB_B{1, day_sel}(:,101:200)];
        
        %% Sort D1
        %A sort
        [~,test_maxBin] = max(STC_AB_1(:,1:100)', [], 1,'includenan');
        %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
        [~,test_sortOrder] = sort(test_maxBin,'ascend');
        
        STC_AB_sort_A = STC_AB_1(test_sortOrder,1:100);
        
        %B
        % [~,test_maxBin] = max(STC_AB_1(:,101:200)', [], 1,'includenan');
        % %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
        % [~,test_sortOrder] = sort(test_maxBin,'ascend');
        
        STC_AB_sort_B = STC_AB_1(test_sortOrder,101:200);
        
        
        %% Sort later day matching
        %A sort
        % [~,test_maxBin] = max(STC_AB_later(:,1:100)', [], 1,'includenan');
        % %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
        % [~,test_sortOrder] = sort(test_maxBin,'ascend');
        
        STC_AB_later_sort_A = STC_AB_later(test_sortOrder,1:100);
        
        %B
        % [~,test_maxBin] = max(STC_AB_later(:,101:200)', [], 1,'includenan');
        % %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
        % [~,test_sortOrder] = sort(test_maxBin,'ascend');
        
        STC_AB_later_sort_B = STC_AB_later(test_sortOrder,101:200);
        
        %% Correlate each session against each other
        %if nay matrix empty
        if ~isempty(STC_AB_sort_A)
            TC_corr_AB_matching_recall_early(aa,day_sel) = mean(diag(corr(STC_AB_sort_A', STC_AB_sort_B','Type','Pearson','rows','all')));
            TC_corr_AB_matching_recall_late(aa,day_sel) = mean(diag(corr(STC_AB_later_sort_A', STC_AB_later_sort_B','Type','Pearson','rows','all')));
            
        else
            TC_corr_AB_matching_recall_early(aa,day_sel) = nan;
            TC_corr_AB_matching_recall_late(aa,day_sel) = nan;
            
        end
        %get the diagnoal correlation values
        if ~isempty(STC_AB_sort_A)
            TC_corr_AB_matching_recall_early_diag{aa,day_sel} = diag(corr(STC_AB_sort_A', STC_AB_sort_B','Type','Pearson','rows','all'));
            TC_corr_AB_matching_recall_late_diag{aa,day_sel} = diag(corr(STC_AB_later_sort_A', STC_AB_later_sort_B','Type','Pearson','rows','all'));
            
        else
            TC_corr_AB_matching_recall_early_diag{aa,day_sel} = nan;
            TC_corr_AB_matching_recall_late_diag{aa,day_sel} = nan;
            
        end
        
        
    end
end

%% Get sem for learning and recall for TC correations -  LEARNING AND RECALL BY ANIMAL AVERAGE

%change this to reflect number of animals per session (by animal)
nb_animal_per_ses_learn = sum(~isnan(TC_corr_AB_matching_learning_early),1)
nb_animal_per_ses_recall = sum(~isnan(TC_corr_AB_matching_recall_early),1)

TC_learning_early_sem = nanstd(TC_corr_AB_matching_learning_early,[],1)./sqrt(nb_animal_per_ses_learn);
TC_learning_late_sem = nanstd(TC_corr_AB_matching_learning_late,[],1)./sqrt(nb_animal_per_ses_learn);

TC_recall_early_sem =nanstd(TC_corr_AB_matching_recall_early,[],1)./sqrt(nb_animal_per_ses_recall);
TC_recall_late_sem = nanstd(TC_corr_AB_matching_recall_late,[],1)./sqrt(nb_animal_per_ses_recall);

%% Get mean and sem by neuron; not animal -LEARNING

%get total neuron count for sessions - learning nb ROI comparison 
for ss=2:9
    %collect all the TC values
    corr_vals_TC_learn{ss} = cell2mat(TC_corr_AB_matching_learning_early_diag(:,ss));
    %remove nans
    corr_vals_TC_learn{ss}(isnan(corr_vals_TC_learn{ss})) = [];
    %get count
    nbROI_TC_learn(ss) = size(corr_vals_TC_learn{ss},1);
end

%get mean across all matching neurons
for ss=2:9
    mean_learn_early_all_ROIs(ss) = nanmean(cell2mat(TC_corr_AB_matching_learning_early_diag(:,ss)));
    mean_learn_late_all_ROIs(ss) = nanmean(cell2mat(TC_corr_AB_matching_learning_late_diag(:,ss)));
end

%get sem across all matching neurons
for ss=2:9
    sem_learn_early_all_ROIs(ss) = nanstd(cell2mat(TC_corr_AB_matching_learning_early_diag(:,ss)),0,1)/sqrt(nbROI_TC_learn(ss));
    sem_learn_late_all_ROIs(ss) = nanstd(cell2mat(TC_corr_AB_matching_learning_late_diag(:,ss)),0,1)/sqrt(nbROI_TC_learn(ss));
end

%% Get mean and sem by neuron; not animal -RECALL

%learning nb ROI comparison
for ss=2:7
    %collect all the TC values
    corr_vals_TC_recall{ss} = cell2mat(TC_corr_AB_matching_recall_early_diag(:,ss));
    %remove nans
    corr_vals_TC_recall{ss}(isnan(corr_vals_TC_recall{ss})) = [];
    %get count
    nbROI_TC_recall(ss) = size(corr_vals_TC_recall{ss},1);
end

%get mean across all matching neurons
for ss=2:7
    mean_recall_early_all_ROIs(ss) = nanmean(cell2mat(TC_corr_AB_matching_recall_early_diag(:,ss)));
    mean_recall_late_all_ROIs(ss) = nanmean(cell2mat(TC_corr_AB_matching_recall_late_diag(:,ss)));
end

%get sem across all matching neurons
for ss=2:7
    sem_recall_early_all_ROIs(ss) = nanstd(cell2mat(TC_corr_AB_matching_recall_early_diag(:,ss)),0,1)/sqrt(nbROI_TC_recall(ss));
    sem_recall_late_all_ROIs(ss) = nanstd(cell2mat(TC_corr_AB_matching_recall_late_diag(:,ss)),0,1)/sqrt(nbROI_TC_recall(ss));
end

%% Construct correlation matrix for learning vs recall - mean of all values

%learning
% for ii=1:6
%     for jj=1:6
%         mean_corr_learn_early(ii,jj) = nanmean(cell2mat(TC_corr_AB_matching_learning_early_diag_all(:,ii,jj)))
%         mean_corr_learn_late(ii,jj) = nanmean(cell2mat(TC_corr_AB_matching_learning_late_diag_all(:,ii,jj)))
%     end
% end
% 
% 
% figure;
% subplot(2,2,1)
%  imagesc(mean_corr_learn_early)
%  hold on
%  caxis([0.2 0.7]) 
%  colormap('jet')
%  subplot(2,2,2)
%  imagesc(mean_corr_learn_late)
%  hold on
%  caxis([0.2 0.7]) 
%  colormap('jet')
 
%% Statistics - do Paired Mann Whitnney U comp on collapsed cells (all neurons)

%learning
for ss=2:9
    [p_learn(ss),~] = signrank(cell2mat(TC_corr_AB_matching_learning_early_diag(:,ss)),cell2mat(TC_corr_AB_matching_learning_late_diag(:,ss)));
end

%recall
for ss = 2:7
    [p_recall(ss),~] = signrank(cell2mat(TC_corr_AB_matching_recall_early_diag(:,ss)),cell2mat(TC_corr_AB_matching_recall_late_diag(:,ss)));
end

%%  Plot mean and sem -  neuron mean

figure('Position',[2058 307 1008 481]);
subplot(1,2,1)
hold on
title('Learning')
xlabel('Session #')
ylabel('Correlation score')
axis square
ylim([0.2 0.9])
%p1 = plot(nanmean(TC_corr_AB_matching_learning_early,1),'k-')
p1 = errorbar(mean_learn_early_all_ROIs(:,2:end),sem_learn_early_all_ROIs(2:end),'k-','LineWidth',2);
%p2 = plot(nanmean(TC_corr_AB_matching_learning_late,1),'k--')
p2 = errorbar(mean_learn_late_all_ROIs(:,2:end),sem_learn_late_all_ROIs(2:end),'Color', [139,0,139]/255,'LineWidth',2);
legend([p1 p2], {'Day 1','Session'},'Location','southwest')
xticks([1:8])
xtickangle(45)
for ii=2:6
    x_labels_learn{ii-1} = [num2str(ii)];
end
xticklabels(x_labels_learn)
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
xlim([0,9])

subplot(1,2,2)
hold on
title('Recall')
xlabel('Session #')
ylabel('Correlation score')
axis square
ylim([0.2 0.9])
%p1 = plot(nanmean(TC_corr_AB_matching_recall_early,1),'k-')
p1 = errorbar(mean_recall_early_all_ROIs(:,2:end),sem_recall_early_all_ROIs(2:end),'k-','LineWidth',2);
%p2 = plot(nanmean(TC_corr_AB_matching_learning_late,1),'k--')
p2 = errorbar(mean_recall_late_all_ROIs(:,2:end),sem_recall_late_all_ROIs(2:end),'Color', [139,0,139]/255,'LineWidth',2);
legend([p1 p2], {'Day 1','Session'},'Location','southeast')
xticks([1:6])
xtickangle(45)
for ii=2:7
    x_labels_learn{ii-1} = [ num2str(ii)];
end
xticklabels(x_labels_learn)
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
xlim([0,7])

%% Plot mean and sem - animal mean (all neurons for given animal matching)

figure('Position',[2058 307 1008 481]);
subplot(1,2,1)
hold on
title('Learning')
xlabel('Session #')
ylabel('Correlation score')
ylim([0 0.8])
%p1 = plot(nanmean(TC_corr_AB_matching_learning_early,1),'k-')
p1 = errorbar(nanmean(TC_corr_AB_matching_learning_early(:,2:end),1),TC_learning_early_sem(2:end),'k-','LineWidth',2);
%p2 = plot(nanmean(TC_corr_AB_matching_learning_late,1),'k--')
p2 = errorbar(nanmean(TC_corr_AB_matching_learning_late(:,2:end),1),TC_learning_late_sem(2:end),'Color', [139,0,139]/255,'LineWidth',2);
legend([p1 p2], {'Day 1','Session'},'Location','southeast')
xticks([1:9])
xtickangle(45)
for ii=2:9
    x_labels_learn{ii-1} = [num2str(ii)];
end
xticklabels(x_labels_learn)
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
xlim([0,10])

subplot(1,2,2)
hold on
title('Recall')
xlabel('Session #')
ylabel('Correlation score')
ylim([0 0.8])
%p1 = plot(nanmean(TC_corr_AB_matching_recall_early,1),'k-')
p1 = errorbar(nanmean(TC_corr_AB_matching_recall_early(:,2:end),1),TC_recall_early_sem(2:end),'k-','LineWidth',2)
%p2 = plot(nanmean(TC_corr_AB_matching_recall_late,1),'k--')
p2 = errorbar(nanmean(TC_corr_AB_matching_recall_late(:,2:end),1),TC_recall_late_sem(2:end),'Color', [139,0,139]/255,'LineWidth',2)
legend([p1 p2], {'Day 1','Session'},'Location','southeast')
xticks([1:6])
xtickangle(45)

for ii=2:7
    x_labels_learn{ii-1} = [ num2str(ii)];
end
xticklabels(x_labels_learn)
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
xlim([0,7])

%% Plot - later STCs
if 0
figure('Position',[2002 44  364  945])
subplot(2,2,1)
imagesc(STC_AB_sort_A)
hold on
title('A')
subplot(2,2,2)
imagesc(STC_AB_sort_B)
hold on
title('B')

subplot(2,2,3)
imagesc(STC_AB_later_sort_A)
hold on
title('Later day')
subplot(2,2,4)
imagesc(STC_AB_later_sort_B)
end

%% Merge D1 STCS for learning - sample matching neurons for visualization
if 0
    d1_STC_AB = [TC_corr_match_learning{1, 1}.tc_corr_match.STC_mat_AB_A{1, 1}(:,1:100), TC_corr_match_learning{1, 1}.tc_corr_match.STC_mat_AB_B{1, 1}(:,1:100)];
    
    %A sort
    [~,test_maxBin] = max(d1_STC_AB(:,1:100)', [], 1,'includenan');
    %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
    [~,test_sortOrder] = sort(test_maxBin,'ascend');
    
    d1_STC_AB_sort_A = d1_STC_AB(test_sortOrder,1:100);
    
    %B
    [~,test_maxBin] = max(d1_STC_AB(:,101:200)', [], 1,'includenan');
    %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
    [~,test_sortOrder] = sort(test_maxBin,'ascend');
    
    d1_STC_AB_sort_B = d1_STC_AB(test_sortOrder,101:200);
end

end

