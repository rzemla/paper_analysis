function [PV_TC_corr] = PV_TC_corr_across_days(animal_data, tunedLogical,registered,options)


%% Get list of matching ROIs across sessions

%matched on all days
%matching_list = registered.multi.assigned_all;
%matched on any day (including nans)
%matching_list = registered.multi.assigned;

%these neurons are just matched (no other criteria applied)
%use matches that were additionally manually filtered for mismatches
matching_list = registered.multi.assigned_filtered;

%these neurons have above list filtered to have a significant place field
%and at least 5 distinct events within the place field
%event,PF, and criteria filtered
matching_list_filtered = registered.multi.matching_list_filtered;

%number of sessions
nb_ses = size(options.sessionSelect,2);
%should yield the same number of correct matrix alignment
%nb_ses = size(matching_list,2);

%assign the selected input trials
options.selectSes = options.selectTrial;


%% Extract mean STC map in each spatial bin (not normalized and not occupancy divided) (100 bins)
%for each session
for ss = options.sessionSelect%1:size(animal_data,2)
    %normalized to A or B trials independently
    A_STC{ss} = animal_data{ss}.Place_cell{options.selectSes(1)}.Spatial_tuning_curve;
    B_STC{ss} = animal_data{ss}.Place_cell{options.selectSes(2)}.Spatial_tuning_curve;
    
    %Gs smoothed, but not normalized (nn) to itself
    %animal_data{ss}.Place_cell{1}.Spatial_Info.rate_map_smooth{8}
    %same as the following
    A_STC_noNorm{ss} = animal_data{ss}.Place_cell{options.selectSes(1)}.Spatial_tuning_curve_no_norm;
    B_STC_noNorm{ss} = animal_data{ss}.Place_cell{options.selectSes(2)}.Spatial_tuning_curve_no_norm;
  
    %dF/F (like STC - normalized) 
    %A_df{ss} = animal_data{ss}.Place_cell{1}.Spatial_tuning_dF;
    %B_df{ss} = animal_data{ss}.Place_cell{2}.Spatial_tuning_dF;
    
    %not occupancy normalized mean dF values across respective A and B
    %trials (not Gaussian smoothed)
    A_df_non_oc{ss} = animal_data{ss}.Place_cell{options.selectSes(1)}.Spatial_Info.mean_dF_map{8};
    B_df_non_oc{ss} = animal_data{ss}.Place_cell{options.selectSes(2)}.Spatial_Info.mean_dF_map{8};
    
end


%% Normalized to max STC value across A and B trials (for each ROI)
%for each session, take max value for each ROI
for ss=options.sessionSelect
    %make cumulative matrix for that session
    comb_STC = [A_STC_noNorm{ss}; B_STC_noNorm{ss}];
    %get min and max value for each ROI (min should all be 0 for STC based on event
    %map)
    min_STC = min(comb_STC,[],1);
    max_STC = max(comb_STC,[],1);
    %get max - min difference for each ROI for [0-1] normalization below
    diff_max_min_STC = max_STC - min_STC;
    
    %feature scale/normalize (0-1 range)
    A_STC_norm{ss} = (A_STC_noNorm{ss} - min_STC)./(diff_max_min_STC);
    B_STC_norm{ss} = (B_STC_noNorm{ss} - min_STC)./(diff_max_min_STC);
end


%% Create combined STC matrix of neurons matching across all sessions
%changed here from fixed 7 value to actual number of sessions 
%generate sample map for 1 vs. 2 day with NaN filled with -1 value
%generate blank ROI by bins matrix

%normalized to A/self trials
matchSTCs_A = zeros(size(matching_list,1),nb_ses*100);
matchSTCs_B = zeros(size(matching_list,1),nb_ses*100);

%normalized across both A and B trials
matchSTCs_A_norm = zeros(size(matching_list,1),nb_ses*100);
matchSTCs_B_norm = zeros(size(matching_list,1),nb_ses*100);

%not self normalized both A and B trials
matchSTCs_A_noNorm = zeros(size(matching_list,1),nb_ses*100);
matchSTCs_B_Nonorm = zeros(size(matching_list,1),nb_ses*100);


%fill the nan assignments to STC to be NaN
for ss=options.sessionSelect
    nan_log = isnan(matching_list(:,ss));
    %norm self
    matchSTCs_A(nan_log,1+(ss-1)*100:ss*100) = NaN;
    matchSTCs_B(nan_log,1+(ss-1)*100:ss*100) = NaN;
    %norm both trials
    matchSTCs_A_norm(nan_log,1+(ss-1)*100:ss*100) = NaN;
    matchSTCs_B_norm(nan_log,1+(ss-1)*100:ss*100) = NaN;
    
    %no norm both trials
    matchSTCs_A_noNorm(nan_log,1+(ss-1)*100:ss*100) = NaN;
    matchSTCs_B_noNorm(nan_log,1+(ss-1)*100:ss*100) = NaN;
end

%fill the matching assignments with trial-normalized STC
for ss=options.sessionSelect
    nan_log = isnan(matching_list(:,ss));
    %normalized original output - each cell normalized to itself for each
    %trial
    matchSTCs_A(~nan_log,1+(ss-1)*100:ss*100) = A_STC{ss}(:,matching_list(~nan_log,ss))';
    matchSTCs_B(~nan_log,1+(ss-1)*100:ss*100) = B_STC{ss}(:,matching_list(~nan_log,ss))';
    %norm both trials
    matchSTCs_A_norm(~nan_log,1+(ss-1)*100:ss*100) = A_STC_norm{ss}(:,matching_list(~nan_log,ss))';
    matchSTCs_B_norm(~nan_log,1+(ss-1)*100:ss*100) = B_STC_norm{ss}(:,matching_list(~nan_log,ss))';
    
    %no norm both trials
    matchSTCs_A_noNorm(~nan_log,1+(ss-1)*100:ss*100) = A_STC_noNorm{ss}(:,matching_list(~nan_log,ss))';
    matchSTCs_B_noNorm(~nan_log,1+(ss-1)*100:ss*100) = B_STC_noNorm{ss}(:,matching_list(~nan_log,ss))';    
end

%sort the values according to day 1 (A)
%maxBin - spatial bin where activity is greatest for each ROI
[~,maxBin_all_A] = max(matchSTCs_A(:,1:100)', [], 1,'includenan');
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,sortOrder_all_A] = sort(maxBin_all_A,'ascend');

%sort the values according to day 1 (B)
%maxBin - spatial bin where activity is greatest for each ROI
[~,maxBin_all_B] = max(matchSTCs_B(:,1:100)', [], 1,'includenan');
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,sortOrder_all_B] = sort(maxBin_all_B,'ascend');

%sorted matrix by day 1 (independent sorting)
matchSTCs_sorted.A = matchSTCs_A(sortOrder_all_A,:);
matchSTCs_sorted.B = matchSTCs_B(sortOrder_all_B,:);

%remove nan values that originate from STC calculation or non-match NaN
%split NaN values from A and B trials (based on 1st session and remove rest) - split in 2 matrices
nan_d1_log.A = isnan(matchSTCs_sorted.A(:,1));
nan_d1_log.B = isnan(matchSTCs_sorted.B(:,1));

matchSTC_sorted_nan.A = matchSTCs_sorted.A(nan_d1_log.A,:);
matchSTC_sorted_nonan.A = matchSTCs_sorted.A(~nan_d1_log.A,:);

matchSTC_sorted_nan.B = matchSTCs_sorted.B(nan_d1_log.B,:);
matchSTC_sorted_nonan.B = matchSTCs_sorted.B(~nan_d1_log.B,:);

%combine matrices with d1 nans below rest of sorted neurons
matchSTC_nan_sorted.A = [matchSTC_sorted_nonan.A; matchSTC_sorted_nan.A];
matchSTC_nan_sorted.B = [matchSTC_sorted_nonan.B; matchSTC_sorted_nan.B];

%% PV correlation analysis across days (relative to D1) for A and B trials

%split matched neuron STC (non_norm) into cell of matrices
for ss = options.sessionSelect
    matchSTCs_nn.A{ss} = matchSTCs_A_noNorm(:,1+(ss-1)*100:ss*100);
    matchSTCs_nn.B{ss} = matchSTCs_B_noNorm(:,1+(ss-1)*100:ss*100);
end

%run population correlation between non-norm event based STCs
for ss = options.sessionSelect(2:end)
    %A corr across days
    PVcorr_rel_d1.A{ss-1} = corr(matchSTCs_nn.A{1},matchSTCs_nn.A{ss}, 'type','Pearson','rows','complete');
    %B corr across days
    PVcorr_rel_d1.B{ss-1} = corr(matchSTCs_nn.B{1},matchSTCs_nn.B{ss}, 'type','Pearson','rows','complete');   
end

%get mean PV score on each day relative to day 1
for ss = options.sessionSelect(2:end)
    %A corr across days
    meanPV_rel_d1.A(ss-1) = nanmean(diag(PVcorr_rel_d1.A{ss-1}));
    %B corr across days
    meanPV_rel_d1.B(ss-1) = nanmean(diag(PVcorr_rel_d1.B{ss-1}));
end

%% PV correlation for all matching neurons (each session against each session) 

%remove all nan values from the matrices before doing the correlation
%allows you to count the number of neurons being correlated

%find those that are nans in both subsets
%this is the logical that says which neurons should be selected
%these indices correspond to the indices in the match matrix
%ROI x bins
for ii = options.sessionSelect
    for jj = options.sessionSelect
        %for A
        keep_nonan_ROIs.A = ~(logical(sum(isnan(matchSTCs_nn.A{ii}),2)) | logical(sum(isnan(matchSTCs_nn.A{jj}),2)));
        %for B
        keep_nonan_ROIs.B = ~(logical(sum(isnan(matchSTCs_nn.B{ii}),2)) | logical(sum(isnan(matchSTCs_nn.B{jj}),2)));
        
        %these are the indices in the matching_list matrix (should be the
        %same b/c these are matching)
        PV_corr_ROI_count.A(ii,jj) = size(find(keep_nonan_ROIs.A == 1),1);
        PV_corr_ROI_count.B(ii,jj) = size(find(keep_nonan_ROIs.B == 1),1);
        
        %get PV correlation of neurons being matched across sesssions
        PV_corr_idx_match.A{ii,jj} = matching_list(keep_nonan_ROIs.A,[ii,jj]);
        PV_corr_idx_match.B{ii,jj} = matching_list(keep_nonan_ROIs.B,[ii,jj]);
        
        %A corr across days
        PVcorr_all_ses_no_nan.A{ii,jj} = corr(matchSTCs_nn.A{ii}(keep_nonan_ROIs.A,:),matchSTCs_nn.A{jj}(keep_nonan_ROIs.A,:), 'type','Pearson');
        %B corr across days
        PVcorr_all_ses_no_nan.B{ii,jj} = corr(matchSTCs_nn.B{ii}(keep_nonan_ROIs.B,:),matchSTCs_nn.B{jj}(keep_nonan_ROIs.B,:), 'type','Pearson');
        
        %save the associated STCs in cell format - for A matches
        PVcorr_all_ses_no_nan_STCs.A{ii,jj}{1} = matchSTCs_nn.A{ii}(keep_nonan_ROIs.A,:);
        PVcorr_all_ses_no_nan_STCs.A{ii,jj}{2} = matchSTCs_nn.A{jj}(keep_nonan_ROIs.A,:);
        
        %save the associated STCs in cell format - for A matches
        PVcorr_all_ses_no_nan_STCs.B{ii,jj}{1} = matchSTCs_nn.B{ii}(keep_nonan_ROIs.B,:);
        PVcorr_all_ses_no_nan_STCs.B{ii,jj}{2} = matchSTCs_nn.B{jj}(keep_nonan_ROIs.B,:);        
    end
end

%QC checked
%check that when manually running correlation between 2 sessions, get the
%same result (same matrix
% test_corr_mat = corr(PVcorr_all_ses_no_nan_STCs.A{3, 2}{1},PVcorr_all_ses_no_nan_STCs.A{3, 2}{2})
% generated_corr_mat = PVcorr_all_ses_no_nan.A{3,2}
% 
% isequal(test_corr_mat,generated_corr_mat)


%run population correlation between non-norm event based STCs
for ii = options.sessionSelect
    for jj = options.sessionSelect
        %A corr across days
        PVcorr_all_ses.A{ii,jj} = corr(matchSTCs_nn.A{ii},matchSTCs_nn.A{jj}, 'type','Pearson','rows','complete');
        %B corr across days
        PVcorr_all_ses.B{ii,jj} = corr(matchSTCs_nn.B{ii},matchSTCs_nn.B{jj}, 'type','Pearson','rows','complete');
    end
end

%% TESTING CODE
%check if you remove nans you get the same result - try one if you remove
%nans
%ses 1 vs ses 2 
% PVcorr_all_ses.A{1, 2}  
% 
% %ROI x bin
% x = matchSTCs_nn.A{1};
% y = matchSTCs_nn.A{2};
% 
% %find indices of nans in both arrays, merge, remove from both, and run the
% %correlation
% 
% %find those that are nans in both subsets
% remove_subset_ROIs = logical(sum(isnan(x),2)) | logical(sum(isnan(y),2));
% 
% %both matrices with nans removed
% x_nonan = x(~remove_subset_ROIs,:);
% y_nonan = y(~remove_subset_ROIs,:);
% 
% PVcorr_test_1_2 = corr(x_nonan,y_nonan, 'type','Pearson');
% 
% %same result as with running corr with complete parameter 
% %use this approach for getting the number of neurons being correlated
% %across days
% isequal(PVcorr_all_ses.A{1, 2}, PVcorr_test_1_2)

%QC check against day 1 generated data - checks out
% isequal(PVcorr_all_ses.A{1,5}, PVcorr_rel_d1.A{4})
% 
% figure
% subplot(1,2,1) 
% imagesc(PVcorr_all_ses.A{1,4})
% subplot(1,2,2)
% imagesc(PVcorr_rel_d1.A{3})

%% TESTING CODE
%compare all the matrices - same result
% for ii = options.sessionSelect
%     for jj = options.sessionSelect
%         %do for A
%         equal_matrix.A(ii,jj) = isequal(PVcorr_all_ses_no_nan.A{ii,jj},PVcorr_all_ses.A{ii,jj});
%         %do for B
%         equal_matrix.B(ii,jj) = isequal(PVcorr_all_ses_no_nan.B{ii,jj},PVcorr_all_ses.B{ii,jj});
%     end
% end

%% Same day correlation for all neurons (non-matching)
%same day correlation between A and B trials
for ss = options.sessionSelect
    %between A and B trials
    PVcorr_same_day.AB{ss} = corr(matchSTCs_nn.A{ss},matchSTCs_nn.B{ss}, 'type','Pearson','rows','complete');
end

%take mean of diagnonal of same day PV correlation
for ss = options.sessionSelect
    %between A and B trials
    meanPV_same_day.AB(ss) = nanmean(diag(PVcorr_same_day.AB{ss}));
end

%% PV correlation (all neurons) - SAME DAY for each session
%expect the value to to down with learning and remain constant with recall
%no filter at the moment
for ss = options.sessionSelect
    PVcorr_all_same_day{ss} = corr(A_STC_noNorm{ss}',B_STC_noNorm{ss}', 'type','Pearson','rows','complete');
    PVcorr_all_same_day_diag(ss,:) = diag(PVcorr_all_same_day{ss});
end


%plot A vs B PV vectors as fxn of track/bin position across sessions
cmap_purple = cbrewer('seq','YlGn',nb_ses);
figure;
hold on
for ss=1:nb_ses
    plot(PVcorr_all_same_day_diag(ss,:),'Color',cmap_purple(ss,:))
end


%plot as raster for the animal
figure
imagesc(PVcorr_all_same_day_diag)
hold on
caxis([0 1])
colormap('jet')

%plot diagonal curves across sessions


%% TC correlation between session of all neurons (non norm STC)

for ss = options.sessionSelect
    TCcorr_all_same_day{ss} = corr(A_STC_noNorm{ss},B_STC_noNorm{ss}, 'type','Pearson','rows','complete');
    TCcorr_all_same_day_diag{ss} = diag(TCcorr_all_same_day{ss});
end

%% Plot the mean PV correlation (same day) as a line plot

figure('Position',[2845 320 430 330]);
hold on
title('Population vector correlation to each other')
ylim([0 1])
ylabel('Mean correlation')
p1 = plot(meanPV_same_day.AB,'m');
legend([p1], 'A vs. B')
xticks(options.sessionSelect)
if options.learning_data ==1
    xticklabels({'1 ', '2', '3','4','5', '6'})
else
    xticklabels({'1 ', '2', '3','6','7', '8', '9'})
end

%% Plot the mean PV correlation (rel to D1) as a line plot

figure('Position',[2845 320 430 330]);
hold on
title('Population vector correlation relative to D1')
ylim([0 1])
ylabel('Mean correlation')
p1 = plot(meanPV_rel_d1.A,'b');
p2 = plot(meanPV_rel_d1.B,'r');
legend([p1 p2], 'A','B')
xticks(options.sessionSelect)
if options.learning_data ==1
    xticklabels({'1 vs. 2', '1 vs.3', '1 vs. 4','1 vs. 5','1 vs. 6'})
else
    xticklabels({'1 vs. 2', '1 vs.3', '1 vs. 6','1 vs. 7','1 vs. 8', '1 vs. 9'})
end

%% TC correlation for neurons (relative to day1)

%SI
%for A tuned matching day 2 day (relative to D1)
%for each session (relative to D1
for ss=options.sessionSelect(2:end)
    %get index of matches for filtered match matrix
    non_nan_idx = find(sum(~isnan(matching_list_filtered.si_Aall_filt_event_filt(:,[1,ss])),2) == 2);
    
    %generate day2 day match matrix
    d2d_match_list = matching_list_filtered.si_Aall_filt_event_filt(non_nan_idx,[1,ss]);
        
    %A corr across days
    TCcorr_rel_d1.si.A{ss-1} = corr(A_STC_noNorm{1}(:,d2d_match_list(:,1)),A_STC_noNorm{ss}(:,d2d_match_list(:,2)), 'type','Pearson','rows','complete');

end

%SI
%for B tuned matching day 2 day (relative to D1)
%for each session (relative to D1
for ss=options.sessionSelect(2:end)
    %get index of matches for filtered match matrix
    non_nan_idx = find(sum(~isnan(matching_list_filtered.si_Ball_filt_event_filt(:,[1,ss])),2) == 2);
    
    %generate day2 day match matrix
    d2d_match_list = matching_list_filtered.si_Ball_filt_event_filt(non_nan_idx,[1,ss]);
        
    %A corr across days
    TCcorr_rel_d1.si.B{ss-1} = corr(B_STC_noNorm{1}(:,d2d_match_list(:,1)),B_STC_noNorm{ss}(:,d2d_match_list(:,2)), 'type','Pearson','rows','complete');

end

%TS
%for A tuned matching day 2 day (relative to D1)
%for each session (relative to D1
for ss=options.sessionSelect(2:end)
    %get index of matches for filtered match matrix
    non_nan_idx = find(sum(~isnan(matching_list_filtered.ts_Aall_filt_event_filt(:,[1,ss])),2) == 2);
    
    %generate day2 day match matrix
    d2d_match_list = matching_list_filtered.ts_Aall_filt_event_filt(non_nan_idx,[1,ss]);
    
    %if not matched
    if ~isempty(d2d_match_list)
        %A corr across days
        TCcorr_rel_d1.ts.A{ss-1} = corr(A_STC_noNorm{1}(:,d2d_match_list(:,1)),A_STC_noNorm{ss}(:,d2d_match_list(:,2)), 'type','Pearson','rows','complete');
    else
        TCcorr_rel_d1.ts.A{ss-1} = nan;
    end
end

%TS
%for B tuned matching day 2 day (relative to D1)
%for each session (relative to D1
for ss=options.sessionSelect(2:end)
    %get index of matches for filtered match matrix
    non_nan_idx = find(sum(~isnan(matching_list_filtered.ts_Ball_filt_event_filt(:,[1,ss])),2) == 2);
    
    %generate day2 day match matrix
    d2d_match_list = matching_list_filtered.ts_Ball_filt_event_filt(non_nan_idx,[1,ss]);
       %if not matched
       if ~isempty(d2d_match_list)
           %A corr across days
           TCcorr_rel_d1.ts.B{ss-1} = corr(B_STC_noNorm{1}(:,d2d_match_list(:,1)),B_STC_noNorm{ss}(:,d2d_match_list(:,2)), 'type','Pearson','rows','complete');
       else
           TCcorr_rel_d1.ts.B{ss-1} = nan;
       end
end

%take the diagonal and mean for TC correlation relative to D1
for ss = options.sessionSelect(2:end)
    %SI
    %A corr across days
    meanTC_rel_d1.si.A(ss-1) = nanmean(diag(TCcorr_rel_d1.si.A{ss-1}));
    %B corr across days
    meanTC_rel_d1.si.B(ss-1) = nanmean(diag(TCcorr_rel_d1.si.B{ss-1}));
    %TS
    %A corr across days
    meanTC_rel_d1.ts.A(ss-1) = nanmean(diag(TCcorr_rel_d1.ts.A{ss-1}));
    %B corr across days
    meanTC_rel_d1.ts.B(ss-1) = nanmean(diag(TCcorr_rel_d1.ts.B{ss-1}));    
end

%% TC correlation for neurons (do this session against asession and compare) - modification of code above - SI

%SI tuned - A trials only
%for A tuned matching day 2 day (relative to D1)
%for each session (relative to D1
for ii = options.sessionSelect
    for ss=options.sessionSelect
        %get index of matches for filtered match matrix
        non_nan_idx = find(sum(~isnan(matching_list_filtered.si_Aall_filt_event_filt(:,[ii,ss])),2) == 2);
        
        %generate day2 day match matrix
        d2d_match_list = matching_list_filtered.si_Aall_filt_event_filt(non_nan_idx,[ii,ss]);
        
        %indexes that match between any 2 days
        TCcorr_idx_match.si.A{ii,ss} = d2d_match_list;
        
    %if not matched
    if ~isempty(d2d_match_list)
        %A corr across days
        TCcorr_all_ses.si.A{ii,ss} = corr(A_STC_noNorm{ii}(:,d2d_match_list(:,1)),A_STC_noNorm{ss}(:,d2d_match_list(:,2)), 'type','Pearson','rows','complete');
                
        %save the associated STCs in cell format - for A matches
        TCcorr_all_ses_STCs.si.A{ii,ss}{1} = A_STC_noNorm{ii}(:,d2d_match_list(:,1));
        TCcorr_all_ses_STCs.si.A{ii,ss}{2} = A_STC_noNorm{ss}(:,d2d_match_list(:,2));
    else
        TCcorr_all_ses.si.A{ii,ss} = nan;
        TCcorr_all_ses_STCs.si.A{ii,ss}{1} = nan;
        TCcorr_all_ses_STCs.si.A{ii,ss}{2} = nan;       
        
    end
    
    end
end

%do sample check here
%QC checked
%check that when manually running correlation between 2 sessions, get the
%same result (same matrix
% test_corr_mat = corr(TCcorr_all_ses_STCs.si.A{3,2}{1},TCcorr_all_ses_STCs.si.A{3,2}{2})
% generated_corr_mat = TCcorr_all_ses.si.A{3,2}
% % 
% isequal(test_corr_mat,generated_corr_mat)

%SI - B trials only
%for B tuned matching day 2 day (relative to D1)
%for each session (relative to D1
for ii = options.sessionSelect
    for ss = options.sessionSelect
        %get index of matches for filtered match matrix
        non_nan_idx = find(sum(~isnan(matching_list_filtered.si_Ball_filt_event_filt(:,[ii,ss])),2) == 2);
        
        %generate day2 day match matrix
        d2d_match_list = matching_list_filtered.si_Ball_filt_event_filt(non_nan_idx,[ii,ss]);
        
        %indexes that match between any 2 days
        TCcorr_idx_match.si.B{ii,ss} = d2d_match_list;
        
        %if not matched
        if ~isempty(d2d_match_list)
            %A corr across days
            TCcorr_all_ses.si.B{ii,ss} = corr(B_STC_noNorm{ii}(:,d2d_match_list(:,1)),B_STC_noNorm{ss}(:,d2d_match_list(:,2)), 'type','Pearson','rows','complete');
            
            %save the associated STCs in cell format - for B matches
            TCcorr_all_ses_STCs.si.B{ii,ss}{1} = B_STC_noNorm{ii}(:,d2d_match_list(:,1));
            TCcorr_all_ses_STCs.si.B{ii,ss}{2} = B_STC_noNorm{ss}(:,d2d_match_list(:,2));
        else
            TCcorr_all_ses.si.B{ii,ss} = nan;
            TCcorr_all_ses_STCs.si.B{ii,ss}{1} = nan;
            TCcorr_all_ses_STCs.si.B{ii,ss}{2} = nan;
            
        end
    end
end


%do sample check here
%QC checked
%check that when manually running correlation between 2 sessions, get the
%same result (same matrix
% test_corr_mat = corr(TCcorr_all_ses_STCs.si.B{3,2}{1},TCcorr_all_ses_STCs.si.B{3,2}{2})
% generated_corr_mat = TCcorr_all_ses.si.B{3,2}
% % % 
% isequal(test_corr_mat,generated_corr_mat)

%extract the neurons counts for each TC correlation by taking the size of
%each matrix
%SI - A trials and B trials
for ii = options.sessionSelect
    for ss = options.sessionSelect
        TCcorr_all_ses_neuron_count.si.A(ii,ss) = size(TCcorr_all_ses.si.A{ii,ss},1);
        TCcorr_all_ses_neuron_count.si.B(ii,ss) = size(TCcorr_all_ses.si.B{ii,ss},1);
    end
end

%QC check - CORRECT FOR SI TUNED - THESE ARE ALREADY min 5 event filtered
%and 1 PF filtered

% %then check that these values are the same as generated from D1 comparisons
% isequal(TCcorr_rel_d1.si.A{4-1},TCcorr_all_ses.si.A{1,4})
% 
% %first check the correlation values are different across time (diagonal)
% figure
% hold on
% title('Diagonals for TC correlations between different sessions - S.I. tuned neurons')
% %plot self
% plot(diag(TCcorr_all_ses.si.A{1,1}))
% %plot 1 vs. 4
% plot(diag(TCcorr_all_ses.si.A{1,4}))
% %plot 3 vs. 4
% plot(diag(TCcorr_all_ses.si.A{3,4}))
% %plot 3 vs. 5
% plot(diag(TCcorr_all_ses.si.A{3,5}))

%% TC correlation for neurons (do this session against every session and compare) - modification of code above - TS

%TS - A
%for each session against each session
for ii=options.sessionSelect
    for ss=options.sessionSelect
        %get index of matches for filtered match matrix
        non_nan_idx = find(sum(~isnan(matching_list_filtered.ts_Aall_filt_event_filt(:,[ii,ss])),2) == 2);
        
        %generate day2 day match matrix
        d2d_match_list = matching_list_filtered.ts_Aall_filt_event_filt(non_nan_idx,[ii,ss]);
        
        %indexes that match between any 2 days
        TCcorr_idx_match.ts.A{ii,ss} = d2d_match_list;        
        
        %if not matched
        if ~isempty(d2d_match_list)
            %A corr across days
            TCcorr_all_ses.ts.A{ii,ss} = corr(A_STC_noNorm{ii}(:,d2d_match_list(:,1)),A_STC_noNorm{ss}(:,d2d_match_list(:,2)), 'type','Pearson','rows','complete');
            
            %save the associated STCs in cell format - for B matches
            TCcorr_all_ses_STCs.ts.A{ii,ss}{1} = A_STC_noNorm{ii}(:,d2d_match_list(:,1));
            TCcorr_all_ses_STCs.ts.A{ii,ss}{2} = A_STC_noNorm{ss}(:,d2d_match_list(:,2));
        else
            TCcorr_all_ses.ts.A{ii,ss} = nan;
            TCcorr_all_ses_STCs.ts.A{ii,ss}{1} = nan;
            TCcorr_all_ses_STCs.ts.A{ii,ss}{2} = nan;
        end
    end
end

%QC checked
%check that when manually running correlation between 2 sessions, get the
%same result (same matrix
% test_corr_mat = corr(TCcorr_all_ses_STCs.ts.A{3,2}{1},TCcorr_all_ses_STCs.ts.A{3,2}{2})
% generated_corr_mat = TCcorr_all_ses.ts.A{3,2}
% % % 
% isequal(test_corr_mat,generated_corr_mat)

%TS - B
%for each session against each session
for ii=options.sessionSelect
    for ss=options.sessionSelect
        %get index of matches for filtered match matrix
        non_nan_idx = find(sum(~isnan(matching_list_filtered.ts_Ball_filt_event_filt(:,[ii,ss])),2) == 2);
        
        %generate day2 day match matrix
        d2d_match_list = matching_list_filtered.ts_Ball_filt_event_filt(non_nan_idx,[ii,ss]);
        
        %matching idxs
        TCcorr_idx_match.ts.B{ii,ss} = d2d_match_list;
        
        %if not matched
        if ~isempty(d2d_match_list)
            %A corr across days
            TCcorr_all_ses.ts.B{ii,ss} = corr(B_STC_noNorm{ii}(:,d2d_match_list(:,1)),B_STC_noNorm{ss}(:,d2d_match_list(:,2)), 'type','Pearson','rows','complete');

            %save the associated STCs in cell format - for B matches
            TCcorr_all_ses_STCs.ts.B{ii,ss}{1} = B_STC_noNorm{ii}(:,d2d_match_list(:,1));
            TCcorr_all_ses_STCs.ts.B{ii,ss}{2} = B_STC_noNorm{ss}(:,d2d_match_list(:,2));
            
        else
            TCcorr_all_ses.ts.B{ii,ss} = nan;
            TCcorr_all_ses_STCs.ts.B{ii,ss}{1} = nan;
            TCcorr_all_ses_STCs.ts.B{ii,ss}{2} = nan;            
            
        end
    end
end

%QC check
% test_corr_mat = corr(TCcorr_all_ses_STCs.ts.B{3,2}{1},TCcorr_all_ses_STCs.ts.B{3,2}{2})
% generated_corr_mat = TCcorr_all_ses.ts.B{3,2}
% % % 
% isequal(test_corr_mat,generated_corr_mat)


%SI - A trials and B trials
for ii = options.sessionSelect
    for ss = options.sessionSelect
        TCcorr_all_ses_neuron_count.ts.A(ii,ss) = size(TCcorr_all_ses.ts.A{ii,ss},1);
        TCcorr_all_ses_neuron_count.ts.B(ii,ss) = size(TCcorr_all_ses.ts.B{ii,ss},1);
    end
end

%% QC check - CORRECT FOR SI TUNED - THESE ARE ALREADY min 5 event filtered
%and 1 PF filtered

%then check that these values are the same as generated from D1 comparisons
% isequal(TCcorr_rel_d1.ts.B{6-1},TCcorr_all_ses.ts.B{1,6})
% 
% %first check the correlation values are different across time (diagonal)
% figure
% hold on
% title('Diagonals for TC correlations between different sessions - T.S. tuned neurons')
% %plot self
% plot(diag(TCcorr_all_ses.ts.A{1,1}))
% %plot 1 vs. 4
% plot(diag(TCcorr_all_ses.ts.A{1,4}))
% %plot 3 vs. 4
% plot(diag(TCcorr_all_ses.ts.A{3,4}))
% %plot 3 vs. 5
% plot(diag(TCcorr_all_ses.ts.A{3,5}))

%% A&B SI tuned neurons  - SI
% CORRELATION FOR A&B neurons will serve as a control to the code that
%is already written for Figure 4H - implemented below

for ii=options.sessionSelect
    for ss=options.sessionSelect
        %get index of matches for filtered match matrix
        non_nan_idx = find(sum(~isnan(matching_list_filtered.si_AB_filt_event_filt(:,[ii,ss])),2) == 2);
        
        %generate day2 day match matrix
        d2d_match_list = matching_list_filtered.si_AB_filt_event_filt(non_nan_idx,[ii,ss]);
        
        %matching idxs
        TCcorr_idx_match.si.AB{ii,ss} = d2d_match_list;
        
        %if there are matching neurons
        if ~isempty(d2d_match_list)
            %generate the matching nonnormalized STCs
            AB_combined_nonNorm_STCs.si{ii,ss} = [A_STC_noNorm{ii}(:,d2d_match_list(:,1))', B_STC_noNorm{ii}(:,d2d_match_list(:,1))';
                A_STC_noNorm{ss}(:,d2d_match_list(:,2))', B_STC_noNorm{ss}(:,d2d_match_list(:,2))'];
            
            %split into separate cells to make correlations easier (same format
            %as above just split into cells
            AB_combined_nonNorm_STCs_cell_split.si{ii,ss}{1,1} = A_STC_noNorm{ii}(:,d2d_match_list(:,1))';
            AB_combined_nonNorm_STCs_cell_split.si{ii,ss}{1,2} = B_STC_noNorm{ii}(:,d2d_match_list(:,1))';
            AB_combined_nonNorm_STCs_cell_split.si{ii,ss}{2,1} = A_STC_noNorm{ss}(:,d2d_match_list(:,2))';
            AB_combined_nonNorm_STCs_cell_split.si{ii,ss}{2,2} = B_STC_noNorm{ss}(:,d2d_match_list(:,2))';
        
            %A vs B TC correlation matrix on first session across days for
            %A&B neurons
            TCcorr_all_ses.si.AB.first{ii,ss} = corr(A_STC_noNorm{ii}(:,d2d_match_list(:,1)),B_STC_noNorm{ii}(:,d2d_match_list(:,1)), 'type','Pearson','rows','complete');
            TCcorr_all_ses.si.AB.second{ii,ss} = corr(A_STC_noNorm{ss}(:,d2d_match_list(:,2)),B_STC_noNorm{ss}(:,d2d_match_list(:,2)), 'type','Pearson','rows','complete');
       
        else %fill this up with nans
            TCcorr_all_ses.si.AB.first{ii,ss} = nan;
            TCcorr_all_ses.si.AB.second{ii,ss} = nan;
            AB_combined_nonNorm_STCs.si{ii,ss} = nan;
            AB_combined_nonNorm_STCs_cell_split.si{ii,ss} = nan;
        end
    end
end

%get neuron number counts for each neuron
%SI A&B neurons - A trials and B trials
for ii = options.sessionSelect
    for ss = options.sessionSelect
        TCcorr_all_ses_neuron_count.si.AB(ii,ss) = size(TCcorr_all_ses.si.AB.first{ii,ss},1);
    end
end

%% QC Check SI
%generate list of day 1 correlations
% for ii=options.sessionSelect
%     first_corr(ii) = mean(diag(TCcorr_all_ses.si.AB.first{1, ii}));
%     second_corr(ii) = mean(diag(TCcorr_all_ses.si.AB.second{1, ii}));
% end
% 
% figure
% hold on
% plot(first_corr)
% plot(second_corr)
% 
% %checked - works
% %QC check above split
% %isequal(AB_combined_nonNorm_STCs{2,4}, cell2mat( AB_combined_nonNorm_STCs_cell_split{2,4}))
% 
% 
% [~,maxBin_all_A] = max(AB_combined_nonNorm_STCs_cell_split.si{ii,ss}{1,1}', [], 1,'includenan');
% %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
% [~,sortOrder_all_A] = sort(maxBin_all_A,'ascend');
% 
% figure
% subplot(2,2,1)
% imagesc(AB_combined_nonNorm_STCs_cell_split.si{ii,ss}{1,1}(sortOrder_all_A,:))
% subplot(2,2,2)
% imagesc(AB_combined_nonNorm_STCs_cell_split.si{ii,ss}{1,2}(sortOrder_all_A,:))
% subplot(2,2,3)
% imagesc(AB_combined_nonNorm_STCs_cell_split.si{ii,ss}{2,1}(sortOrder_all_A,:))
% subplot(2,2,4)
% imagesc(AB_combined_nonNorm_STCs_cell_split.si{ii,ss}{2,2}(sortOrder_all_A,:))


%% A&B TS tuned neurons  
% CORRELATION FOR A&B neurons will serve as a control to the code that
%is already written for Figure 4H - implemented below

for ii=options.sessionSelect
    for ss=options.sessionSelect
        %get index of matches for filtered match matrix
        non_nan_idx = find(sum(~isnan(matching_list_filtered.ts_AB_filt_event_filt(:,[ii,ss])),2) == 2);
        
        %generate day2 day match matrix
        d2d_match_list = matching_list_filtered.ts_AB_filt_event_filt(non_nan_idx,[ii,ss]);
        
        %matching idxs
        TCcorr_idx_match.ts.AB{ii,ss} = d2d_match_list;
        
        
        %if there are matching neurons
        if ~isempty(d2d_match_list)
            %generate the matching nonnormalized STCs
            AB_combined_nonNorm_STCs.ts{ii,ss} = [A_STC_noNorm{ii}(:,d2d_match_list(:,1))', B_STC_noNorm{ii}(:,d2d_match_list(:,1))';
                A_STC_noNorm{ss}(:,d2d_match_list(:,2))', B_STC_noNorm{ss}(:,d2d_match_list(:,2))'];
            
            %split into separate cells to make correlations easier (same format
            %as above just split into cells
            AB_combined_nonNorm_STCs_cell_split.ts{ii,ss}{1,1} = A_STC_noNorm{ii}(:,d2d_match_list(:,1))';
            AB_combined_nonNorm_STCs_cell_split.ts{ii,ss}{1,2} = B_STC_noNorm{ii}(:,d2d_match_list(:,1))';
            AB_combined_nonNorm_STCs_cell_split.ts{ii,ss}{2,1} = A_STC_noNorm{ss}(:,d2d_match_list(:,2))';
            AB_combined_nonNorm_STCs_cell_split.ts{ii,ss}{2,2} = B_STC_noNorm{ss}(:,d2d_match_list(:,2))';
        
            %A vs B TC correlation matrix on first session across days for
            %A&B neurons
            TCcorr_all_ses.ts.AB.first{ii,ss} = corr(A_STC_noNorm{ii}(:,d2d_match_list(:,1)),B_STC_noNorm{ii}(:,d2d_match_list(:,1)), 'type','Pearson','rows','complete');
            TCcorr_all_ses.ts.AB.second{ii,ss} = corr(A_STC_noNorm{ss}(:,d2d_match_list(:,2)),B_STC_noNorm{ss}(:,d2d_match_list(:,2)), 'type','Pearson','rows','complete');
       
        else %fill this up with nans
            TCcorr_all_ses.ts.AB.first{ii,ss} = nan;
            TCcorr_all_ses.ts.AB.second{ii,ss} = nan;
            AB_combined_nonNorm_STCs.ts{ii,ss} = nan;
            AB_combined_nonNorm_STCs_cell_split.ts{ii,ss} = nan;
        end
    end
end


%get neuron number counts for each neuron
%SI A&B neurons - A trials and B trials
for ii = options.sessionSelect
    for ss = options.sessionSelect
        TCcorr_all_ses_neuron_count.ts.AB(ii,ss) = size(TCcorr_all_ses.ts.AB.first{ii,ss},1);
    end
end

%% QC Check TS
%generate list of day 1 correlations
% for ii=options.sessionSelect
%     first_corr(ii) = mean(diag(TCcorr_all_ses.ts.AB.first{1, ii}));
%     second_corr(ii) = mean(diag(TCcorr_all_ses.ts.AB.second{1, ii}));
% end
% 
% figure
% hold on
% plot(first_corr)
% plot(second_corr)
% 
% %checked - works
% %QC check above split
% %isequal(AB_combined_nonNorm_STCs{2,4}, cell2mat( AB_combined_nonNorm_STCs_cell_split{2,4}))
% 
% 
% [~,maxBin_all_A] = max(AB_combined_nonNorm_STCs_cell_split.ts{ii,ss}{1,1}', [], 1,'includenan');
% %sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
% [~,sortOrder_all_A] = sort(maxBin_all_A,'ascend');
% 
% figure
% subplot(2,2,1)
% imagesc(AB_combined_nonNorm_STCs_cell_split.ts{ii,ss}{1,1}(sortOrder_all_A,:))
% subplot(2,2,2)
% imagesc(AB_combined_nonNorm_STCs_cell_split.ts{ii,ss}{1,2}(sortOrder_all_A,:))
% subplot(2,2,3)
% imagesc(AB_combined_nonNorm_STCs_cell_split.ts{ii,ss}{2,1}(sortOrder_all_A,:))
% subplot(2,2,4)
% imagesc(AB_combined_nonNorm_STCs_cell_split.ts{ii,ss}{2,2}(sortOrder_all_A,:))

%% Plot - tuning curve correlations

figure('Position',[2040 320 800 330]);
subplot(1,2,1)
hold on
title('Tuning curve correlation relative to D1 (S.I.)')
ylim([0 1])
ylabel('Mean correlation')
p1 = plot(meanTC_rel_d1.si.A,'b');
p2 = plot(meanTC_rel_d1.si.B,'r');
legend([p1 p2], 'A','B')
xticks(options.sessionSelect)
if options.learning_data ==1
    xticklabels({'1 vs. 2', '1 vs.3', '1 vs. 4','1 vs. 5','1 vs. 6'})

else
    xticklabels({'1 vs. 2', '1 vs.3', '1 vs. 6','1 vs. 7','1 vs. 8', '1 vs. 9'})
end

subplot(1,2,2)
hold on
title('Tuning curve correlation relative to D1 (T.S.)')
ylim([0 1])
ylabel('Mean correlation')
p1 = plot(meanTC_rel_d1.ts.A,'b');
p2 = plot(meanTC_rel_d1.ts.B,'r');
legend([p1 p2], 'A','B')
xticks(options.sessionSelect)
if options.learning_data ==1
    xticklabels({'1 vs. 2', '1 vs.3', '1 vs. 4','1 vs. 5','1 vs. 6'})
else
    xticklabels({'1 vs. 2', '1 vs.3', '1 vs. 6','1 vs. 7','1 vs. 8', '1 vs. 9'})
end

%% Export data for cumulative analysis for all animals

%%%% THIS IS IMPORTANT %%%%
%TC CORRELATION DATA FOR EXPORT (BOTH S.I and T.S in substruct)
PV_TC_corr.TCcorr_all_ses = TCcorr_all_ses;
PV_TC_corr.TCcorr_all_ses_neuron_count = TCcorr_all_ses_neuron_count;

%STCs used to generate TC correlations for A and B cells across days 
PV_TC_corr.TCcorr_all_ses_STCs = TCcorr_all_ses_STCs;

%A&B tuned neurons compared against each session for sis and ts
PV_TC_corr.AB_combined_nonNorm_STCs = AB_combined_nonNorm_STCs;
PV_TC_corr.AB_combined_nonNorm_STCs_cell_split = AB_combined_nonNorm_STCs_cell_split;

%idxs of TC correlated neurons
PV_TC_corr.TCcorr_idx_match = TCcorr_idx_match;

%PV CORRELATION DATA FOR EXPORT (ALL NEURONS)
%PV correlation matrices - same results
PV_TC_corr.PVcorr_all_ses_no_nan = PVcorr_all_ses_no_nan;
PV_TC_corr.PVcorr_all_ses = PVcorr_all_ses;

%numbers of ROIs PV correlated against any 2 sessions
PV_TC_corr.PV_corr_ROI_count = PV_corr_ROI_count;
%idxs of neurons being matched across days
PV_TC_corr.PV_corr_idx_match = PV_corr_idx_match;

%STCs used to generate PV correlations
PV_TC_corr.PVcorr_all_ses_no_nan_STCs = PVcorr_all_ses_no_nan_STCs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%previous data

%PV correlations relative to day 1
PV_TC_corr.meanPV_rel_d1 = meanPV_rel_d1;

%PV correlations same day
PV_TC_corr.meanPV_same_day = meanPV_same_day;
%PV correlation - same day - all neurons
PV_TC_corr.PVcorr_all_same_day = PVcorr_all_same_day;
%PV correlation - same day - all neurons - diagnonal values (bin by bin)
PV_TC_corr.PVcorr_all_same_day_diag = PVcorr_all_same_day_diag;
%
%TC correlations relative to day 1
PV_TC_corr.meanTC_rel_d1 = meanTC_rel_d1;
%contains the correlation matrices for each day relative to day1 for thos
%neurons that matched (both si and ts criteria)
PV_TC_corr_rel_d1_mat_nonnon = TCcorr_rel_d1;

%Same day A vs B. TC correlations for each neuron (all neurons)
PV_TC_corr.TC_same_day_mat = TCcorr_all_same_day;
PV_TC_corr.TC_same_day_diag = TCcorr_all_same_day_diag;

end

