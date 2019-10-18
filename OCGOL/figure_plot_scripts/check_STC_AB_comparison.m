function [outputArg1,outputArg2] = check_STC_AB_comparison(TC_corr_match_learning,TC_corr_match_recall)


%day 1
%TC_corr_match_learning{1, 1}.tc_corr_match.STC_mat_AB_A{1, 1}
%TC_corr_match_learning{1, 1}.tc_corr_match.STC_mat_AB_B{1, 1}

nb_recall_animals = size(TC_corr_match_recall,2);
nb_learn_animals = size(TC_corr_match_learning,2);

%% Merge D1 STCS for learning

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

%% Relative to D1 comparison - learning - relative to day 1 only
for aa=1:3
    
    %session relative to session 1 (d1)
    for day_sel=2:6
        
        STC_AB_1 = [TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_A{1, day_sel}(:,1:100), TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_B{1, day_sel}(:,1:100)];
        STC_AB_later = [TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_A{1, day_sel}(:,101:200), TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_B{1, day_sel}(:,101:200)];
        
        %         STC_AB_1 = [TC_corr_match_recall{aa}.tc_corr_match.STC_mat_AB_A{1, day_sel}(:,1:100), TC_corr_match_recall{aa}.tc_corr_match.STC_mat_AB_B{1, day_sel}(:,1:100)];
        %         STC_AB_later = [TC_corr_match_recall{aa}.tc_corr_match.STC_mat_AB_A{1, day_sel}(:,101:200), TC_corr_match_recall{aa}.tc_corr_match.STC_mat_AB_B{1, day_sel}(:,101:200)];
        %
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
            TC_corr_AB_matching_learning_early(aa,day_sel) = mean(diag(corr(STC_AB_sort_A', STC_AB_sort_B','Type','Pearson','rows','all')));
            TC_corr_AB_matching_learning_late(aa,day_sel) = mean(diag(corr(STC_AB_later_sort_A', STC_AB_later_sort_B','Type','Pearson','rows','all')));
        else
            TC_corr_AB_matching_learning_early(aa,day_sel) = nan;
            TC_corr_AB_matching_learning_late(aa,day_sel) = nan;
            
        end
        %get diagonal correlation values
        if ~isempty(STC_AB_sort_A)
            TC_corr_AB_matching_learning_early_diag{aa,day_sel} = diag(corr(STC_AB_sort_A', STC_AB_sort_B','Type','Pearson','rows','all'));
            TC_corr_AB_matching_learning_late_diag{aa,day_sel} = diag(corr(STC_AB_later_sort_A', STC_AB_later_sort_B','Type','Pearson','rows','all'));
        else
            TC_corr_AB_matching_learning_early_diag{aa,day_sel} = nan;
            TC_corr_AB_matching_learning_late_diag{aa,day_sel} = nan;
            
        end
        
    end
end

%% Cross-comparison to all other sessions -learning 
for aa=1:3
    
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

%% Get sem for learning and recall for TC correations

TC_learning_early_sem =nanstd(TC_corr_AB_matching_learning_early,[],1)./sqrt(size(TC_corr_AB_matching_learning_early,1));
TC_learning_late_sem = nanstd(TC_corr_AB_matching_learning_late,[],1)./sqrt(size(TC_corr_AB_matching_learning_late,1));

TC_recall_early_sem =nanstd(TC_corr_AB_matching_recall_early,[],1)./sqrt(size(TC_corr_AB_matching_recall_early,1));
TC_recall_late_sem = nanstd(TC_corr_AB_matching_recall_late,[],1)./sqrt(size(TC_corr_AB_matching_recall_late,1));


%% Construct correlation matrix for learning vs recall - mean of all values

%learning
for ii=1:6
    for jj=1:6
        mean_corr_learn_early(ii,jj) = nanmean(cell2mat(TC_corr_AB_matching_learning_early_diag_all(:,ii,jj)))
         mean_corr_learn_late(ii,jj) = nanmean(cell2mat(TC_corr_AB_matching_learning_late_diag_all(:,ii,jj)))
    end
end


figure;
subplot(2,2,1)
 imagesc(mean_corr_learn_early)
 hold on
 caxis([0.2 0.7]) 
 colormap('jet')
 subplot(2,2,2)
 imagesc(mean_corr_learn_late)
 hold on
 caxis([0.2 0.7]) 
 colormap('jet')
 
%% Do Paired Mann Whitnney U comp on collapsed cells (all neurons

%learning
for test_ses=2:6
    [p_learn(test_ses),~] = signrank(cell2mat(TC_corr_AB_matching_learning_early_diag(:,test_ses)),cell2mat(TC_corr_AB_matching_learning_late_diag(:,test_ses)));
end

%learning nb ROi comparison
for test_ses=2:6
    nbROI_TC_learn(test_ses) = length(cell2mat(TC_corr_AB_matching_learning_early_diag(:,test_ses)))
end

%recall
for test_ses = 2:7
    [p_recall(test_ses),~] = signrank(cell2mat(TC_corr_AB_matching_recall_early_diag(:,test_ses)),cell2mat(TC_corr_AB_matching_recall_late_diag(:,test_ses)));
end

%learning nb ROi comparison
for test_ses=2:7
    nbROI_TC_recall(test_ses) = length(cell2mat(TC_corr_AB_matching_recall_early_diag(:,test_ses)))
end
%% Get mean and sem by neuron; not animal -LEARNING
%get mean across all matching neurons
for test_ses=2:6
    mean_learn_early_all_ROIs(test_ses) = nanmean(cell2mat(TC_corr_AB_matching_learning_early_diag(:,test_ses)));
    mean_learn_late_all_ROIs(test_ses) = nanmean(cell2mat(TC_corr_AB_matching_learning_late_diag(:,test_ses)));
end

%get sem across all matching neurons
for test_ses=2:6
    sem_learn_early_all_ROIs(test_ses) = nanstd(cell2mat(TC_corr_AB_matching_learning_early_diag(:,test_ses)),0,1)/sqrt(nbROI_TC_learn(test_ses));
    sem_learn_late_all_ROIs(test_ses) = nanstd(cell2mat(TC_corr_AB_matching_learning_late_diag(:,test_ses)),0,1)/sqrt(nbROI_TC_learn(test_ses));
end

%% Get mean and sem by neuron; not animal -RECALL
%get mean across all matching neurons
for test_ses=2:7
    mean_recall_early_all_ROIs(test_ses) = nanmean(cell2mat(TC_corr_AB_matching_recall_early_diag(:,test_ses)));
    mean_recall_late_all_ROIs(test_ses) = nanmean(cell2mat(TC_corr_AB_matching_recall_late_diag(:,test_ses)));
end

%get sem across all matching neurons
for test_ses=2:7
    sem_recall_early_all_ROIs(test_ses) = nanstd(cell2mat(TC_corr_AB_matching_recall_early_diag(:,test_ses)),0,1)/sqrt(nbROI_TC_recall(test_ses));
    sem_recall_late_all_ROIs(test_ses) = nanstd(cell2mat(TC_corr_AB_matching_recall_late_diag(:,test_ses)),0,1)/sqrt(nbROI_TC_recall(test_ses));
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
legend([p1 p2], {'Reference','Session'},'Location','southeast')
xticks([1:5])
xtickangle(45)
for ii=2:6
    x_labels_learn{ii-1} = [num2str(ii)];
end
xticklabels(x_labels_learn)
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
xlim([0,6])

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
legend([p1 p2], {'Reference','Session'},'Location','southeast')
xticks([1:6])
xtickangle(45)
for ii=2:7
    x_labels_learn{ii-1} = [ num2str(ii)];
end
xticklabels(x_labels_learn)
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
xlim([0,7])

%% Plot mean and sem - animal mean across all neurons

figure('Position',[2058 307 1008 481]);
subplot(1,2,1)
hold on
title('Learning')
xlabel('Session #')
ylabel('Correlation score')
ylim([0.2 0.8])
%p1 = plot(nanmean(TC_corr_AB_matching_learning_early,1),'k-')
p1 = errorbar(nanmean(TC_corr_AB_matching_learning_early(:,2:end),1),TC_learning_early_sem(2:end),'k-','LineWidth',2);
%p2 = plot(nanmean(TC_corr_AB_matching_learning_late,1),'k--')
p2 = errorbar(nanmean(TC_corr_AB_matching_learning_late(:,2:end),1),TC_learning_late_sem(2:end),'Color', [139,0,139]/255,'LineWidth',2);
legend([p1 p2], {'Reference','Session'},'Location','southeast')
xticks([1:5])
xtickangle(45)
for ii=2:6
    x_labels_learn{ii-1} = [num2str(ii)];
end
xticklabels(x_labels_learn)
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
xlim([0,6])

subplot(1,2,2)
hold on
title('Recall')
xlabel('Session #')
ylabel('Correlation score')
ylim([0.2 0.8])
%p1 = plot(nanmean(TC_corr_AB_matching_recall_early,1),'k-')
p1 = errorbar(nanmean(TC_corr_AB_matching_recall_early(:,2:end),1),TC_recall_early_sem(2:end),'k-','LineWidth',2)
%p2 = plot(nanmean(TC_corr_AB_matching_recall_late,1),'k--')
p2 = errorbar(nanmean(TC_corr_AB_matching_recall_late(:,2:end),1),TC_recall_late_sem(2:end),'Color', [139,0,139]/255,'LineWidth',2)
legend([p1 p2], {'Reference','Session'},'Location','southeast')
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

