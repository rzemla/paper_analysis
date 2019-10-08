function [outputArg1,outputArg2] = trajectory_analysis(TC_corr_match_learning,TC_corr_match_recall)


%% Extract matching STCs

for aa=1:3
    
    %session relative to session 1st
    for day_sel = 1:6
        
        %session relative to 2nd session - CONTINUE HERE
        for day_sel_2 = 1:6
            
            STC_AB_1{aa,day_sel,day_sel_2} = [TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_A{day_sel, day_sel_2}(:,1:100), TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_B{day_sel, day_sel_2}(:,1:100)];
            STC_AB_later{aa,day_sel,day_sel_2} = [TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_A{day_sel, day_sel_2}(:,101:200), TC_corr_match_learning{aa}.tc_corr_match.STC_mat_AB_B{day_sel, day_sel_2}(:,101:200)];
            
            %% Sort D1
            %A sort
            %[~,test_maxBin] = max(STC_AB_1(:,1:100)', [], 1,'includenan');
            %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
            %[~,test_sortOrder] = sort(test_maxBin,'ascend');
            
            %STC_AB_sort_A = STC_AB_1(test_sortOrder,1:100);
            
            %B
            % [~,test_maxBin] = max(STC_AB_1(:,101:200)', [], 1,'includenan');
            % %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
            % [~,test_sortOrder] = sort(test_maxBin,'ascend');
            
            %STC_AB_sort_B = STC_AB_1(test_sortOrder,101:200);
            
            
            %% Sort later day matching
            %A sort
            % [~,test_maxBin] = max(STC_AB_later(:,1:100)', [], 1,'includenan');
            % %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
            % [~,test_sortOrder] = sort(test_maxBin,'ascend');
            
            %STC_AB_later_sort_A = STC_AB_later(test_sortOrder,1:100);
            
            %B
            % [~,test_maxBin] = max(STC_AB_later(:,101:200)', [], 1,'includenan');
            % %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
            % [~,test_sortOrder] = sort(test_maxBin,'ascend');
            
            %STC_AB_later_sort_B = STC_AB_later(test_sortOrder,101:200);
            

            
        end
    end
end





%% Sort STCs by A and B trials

aa=2
ii= 1
jj=6

%A sort
[~,maxBin_A] = max(STC_AB_1{aa,ii,jj}(:,1:100)', [], 1,'includenan');
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,sortOrder_A] = sort(maxBin_A,'ascend');
%B sort
[~,maxBin_B] = max(STC_AB_1{aa,ii,jj}(:,101:200)', [], 1,'includenan');
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,sortOrder_B] = sort(maxBin_B,'ascend');


STC_AB_1_sort_A = STC_AB_1{aa,ii,jj}(sortOrder_A,1:100);
STC_AB_1_sort_B = STC_AB_1{aa,ii,jj}(sortOrder_A,101:200);


[~,maxBin_A] = max(STC_AB_later{aa,ii,jj}(:,1:100)', [], 1,'includenan');
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,sortOrder_A] = sort(maxBin_A,'ascend');
%B sort
[~,maxBin_B] = max(STC_AB_later{aa,ii,jj}(:,101:200)', [], 1,'includenan');
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,sortOrder_B] = sort(maxBin_B,'ascend');


STC_AB_later_sort_A = STC_AB_later{aa,ii,jj}(sortOrder_A,1:100);
STC_AB_later_sort_B = STC_AB_later{aa,ii,jj}(sortOrder_A,101:200);

%% Plot side by side 

figure
subplot(2,2,1)
imagesc(STC_AB_1_sort_A)
caxis([0 1])
colormap('jet')

subplot(2,2,2)
imagesc(STC_AB_1_sort_B)
caxis([0 1])
colormap('jet')

subplot(2,2,3)
imagesc(STC_AB_later_sort_A)
caxis([0 1])
colormap('jet')

subplot(2,2,4) 
imagesc(STC_AB_later_sort_B)
caxis([0 1])
colormap('jet')


end








