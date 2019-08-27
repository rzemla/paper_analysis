function [outputArg1,outputArg2] = non_norm_matching_STC_rasters_learning(animal_data, tunedLogical,registered,options)


%% Get list of matching ROIs across sessions

%matched on all days
%matching_list = registered.multi.assigned_all;
%matched on any day (including nans)
%matching_list = registered.multi.assigned;
%use matches that were additionally manually filtered for mismatches
matching_list = registered.multi.assigned_filtered;


%% Extract mean STC map in each spatial bin (not normalized and not occupancy divided) (100 bins)
%for each session
for ss = options.sessionSelect%1:size(animal_data,2)
    %normalized to A or B trials independently
    A_STC{ss} = animal_data{ss}.Place_cell{4}.Spatial_tuning_curve;
    B_STC{ss} = animal_data{ss}.Place_cell{5}.Spatial_tuning_curve;
    
    %Gs smoothed, but not normalized (nn) to itself
    %animal_data{ss}.Place_cell{1}.Spatial_Info.rate_map_smooth{8}
    %same as the following
    A_STC_noNorm{ss} = animal_data{ss}.Place_cell{4}.Spatial_tuning_curve_no_norm;
    B_STC_noNorm{ss} = animal_data{ss}.Place_cell{5}.Spatial_tuning_curve_no_norm;
  
    %dF/F (like STC - normalized) 
    %A_df{ss} = animal_data{ss}.Place_cell{1}.Spatial_tuning_dF;
    %B_df{ss} = animal_data{ss}.Place_cell{2}.Spatial_tuning_dF;
    
    %not occupancy normalized mean dF values across respective A and B
    %trials (not Gaussian smoothed)
    A_df_non_oc{ss} = animal_data{ss}.Place_cell{4}.Spatial_Info.mean_dF_map{8};
    B_df_non_oc{ss} = animal_data{ss}.Place_cell{5}.Spatial_Info.mean_dF_map{8};
    
end

%take raw dF/F values, Gaussian smoothed and occupancy normalize

%normalized occupancy
%animal_data{1, 1}.Place_cell{1, 1}.Spatial_Info.proba_bin  


%% Normalized to max STC value across A and B trials (for each ROI)
%for each session, take max value for each ROI
for ss=1:6
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

%% Debug/qc related
%isequal(animal_data{1, 1}.Place_cell{1, 1}.Spatial_Info.rate_map_smooth{1, 8}(:,35),...
 %        animal_data{1, 1}.Place_cell{1, 1}.Spatial_tuning_curve_no_norm(:,35))



%% Plot STC of neurons matching across all sessions - black STC for NaNs

%colormap with black values below 0;
%default call returns 64 values for cmap range
% cmap = colormap('jet');
% %black map
% cmap_black = repmat([0 0 0],63,1);
% %combined map
% cmap_extended = [cmap_black; cmap];

%generate sample map for 1 vs. 2 day with NaN filled with -1 value
%generate blank ROI by bins matrix

%normalized to A/self trials
matchSTCs_A = zeros(size(matching_list,1),7*100);
matchSTCs_B = zeros(size(matching_list,1),7*100);

%normalized across both A and B trials
matchSTCs_A_norm = zeros(size(matching_list,1),7*100);
matchSTCs_B_norm = zeros(size(matching_list,1),7*100);

%not self normalized both A and B trials
matchSTCs_A_noNorm = zeros(size(matching_list,1),7*100);
matchSTCs_B_Nonorm = zeros(size(matching_list,1),7*100);


%fill the nan assignments to STC to be NaN
for ss=1:6
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
for ss=1:6
    nan_log = isnan(matching_list(:,ss));
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
%maxBin - spatial bin where activity is greatest for each ROI - returns
%bin index with max value
[~,maxBin_all_A] = max(matchSTCs_A(:,1:100)', [], 1,'includenan');
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
%sort according to max bin values for all ROIs
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


%% Sort B trials relative to A trials
matchSTCs_sorted.BrelA = matchSTCs_B(sortOrder_all_A,:);
matchSTC_sorted_nan.BrelA = matchSTCs_sorted.BrelA(nan_d1_log.A,:);
matchSTC_sorted_nonan.BrelA = matchSTCs_sorted.BrelA(~nan_d1_log.A,:);
matchSTC_nan_sorted.BrelA = [matchSTC_sorted_nonan.BrelA; matchSTC_sorted_nan.BrelA];


%% Plot the raster sorted independently

figure
%A
subplot(1,2,1)
%create blank alpha shading matrix where 
imAlpha=ones(size(matchSTC_nan_sorted.A));
imAlpha(isnan(matchSTC_nan_sorted.A))=0;
imagesc(matchSTC_nan_sorted.A,'AlphaData',imAlpha);
%set background axis color to black
set(gca,'color',0*[1 1 1]);
%set colormap to 
colormap(gca,'jet');
%B
subplot(1,2,2)
%create blank alpha shading matrix where 
imAlpha=ones(size(matchSTC_nan_sorted.B));
imAlpha(isnan(matchSTC_nan_sorted.B))=0;
imagesc(matchSTC_nan_sorted.B,'AlphaData',imAlpha);
%set background axis color to black
set(gca,'color',0*[1 1 1]);
%set colormap to 
colormap(gca,'jet');

%% Plot rasters of B ROIs/trial sorted relative to A trials

figure
%A
subplot(1,2,1)
%create blank alpha shading matrix where
%set equal (max) transparency across the matrix
imAlpha=ones(size(matchSTC_nan_sorted.A));
%set transparency of nan values to 0 (non transparency/min)
imAlpha(isnan(matchSTC_nan_sorted.A))=0;
%plot raster with transparency matrix set
imagesc(matchSTC_nan_sorted.A,'AlphaData',imAlpha);
%set background axis color to black
set(gca,'color',0*[1 1 1]);
%set colormap to 
colormap(gca,'jet');

%B
subplot(1,2,2)
%create blank alpha shading matrix where 
%set equal (max) transparency across the matrix
imAlpha=ones(size(matchSTC_nan_sorted.BrelA));
%set transparency of nan values to 0 (non transparency/min)
imAlpha(isnan(matchSTC_nan_sorted.BrelA))=0;
%plot raster with transparency matrix set
imagesc(matchSTC_nan_sorted.BrelA,'AlphaData',imAlpha);
%set background axis color to black
set(gca,'color',0*[1 1 1]);
%set colormap to 
colormap(gca,'jet');

%% Plot rasters normalized for each ROI across both sessions
%use sort order for self-norm A trials
matchSTCs_norm_sorted.A = matchSTCs_A_norm(sortOrder_all_A,:);
matchSTCs_norm_sorted.B = matchSTCs_B_norm(sortOrder_all_A,:);

%push nans to bottom of raster
matchSTC_norm_sorted_nan.A = matchSTCs_norm_sorted.A(nan_d1_log.A,:);
matchSTC_norm_sorted_nonan.A = matchSTCs_norm_sorted.A(~nan_d1_log.A,:);

matchSTC_norm_sorted_nan.BrelA = matchSTCs_norm_sorted.B(nan_d1_log.A,:);
matchSTC_norm_sorted_nonan.BrelA = matchSTCs_norm_sorted.B(~nan_d1_log.A,:);

%combine matrices with d1 nans below rest of sorted neurons
matchSTC_norm_nan_sorted.A = [matchSTC_norm_sorted_nonan.A; matchSTC_norm_sorted_nan.A];
matchSTC_norm_nan_sorted.BrelA = [matchSTC_norm_sorted_nonan.BrelA; matchSTC_norm_sorted_nan.BrelA];

%% Plot normalized and A sorted STC rasters
figure
%A
subplot(1,2,1)
%create blank alpha shading matrix where
%set equal (max) transparency across the matrix
imAlpha=ones(size(matchSTC_norm_nan_sorted.A));
%set transparency of nan values to 0 (non transparency/min)
imAlpha(isnan(matchSTC_norm_nan_sorted.A))=0;
%plot raster with transparency matrix set
imagesc(matchSTC_norm_nan_sorted.A,'AlphaData',imAlpha);
%set background axis color to black
set(gca,'color',0*[1 1 1]);
%set colormap to 
colormap(gca,'jet');

%B
subplot(1,2,2)
%create blank alpha shading matrix where 
%set equal (max) transparency across the matrix
imAlpha=ones(size(matchSTC_norm_nan_sorted.BrelA));
%set transparency of nan values to 0 (non transparency/min)
imAlpha(isnan(matchSTC_norm_nan_sorted.BrelA))=0;
%plot raster with transparency matrix set
imagesc(matchSTC_norm_nan_sorted.BrelA,'AlphaData',imAlpha);
%set background axis color to black
set(gca,'color',0*[1 1 1]);
%set colormap to 
colormap(gca,'jet');

%% PV correlation analysis across days (relative to D1) for A and B trials

%split matched neuron STC (non_norm) into cell of matrices
for ss = 1:6
    matchSTCs_nn.A{ss} = matchSTCs_A_noNorm(:,1+(ss-1)*100:ss*100);
    matchSTCs_nn.B{ss} = matchSTCs_B_noNorm(:,1+(ss-1)*100:ss*100);
end

%run population correlation between non-norm event based STCs
for ss = 2:6
    %A corr across days
    PVcorr_rel_d1.A{ss-1} = corr(matchSTCs_nn.A{1},matchSTCs_nn.A{ss}, 'type','Pearson','rows','complete');
    %B corr across days
    PVcorr_rel_d1.B{ss-1} = corr(matchSTCs_nn.B{1},matchSTCs_nn.B{ss}, 'type','Pearson','rows','complete');   
end

%get mean PV score on each day relative to day 1
for ss = 2:6
    %A corr across days
    meanPV_rel_d1.A(ss-1) = nanmean(diag(PVcorr_rel_d1.A{ss-1}));
    %B corr across days
    meanPV_rel_d1.B(ss-1) = nanmean(diag(PVcorr_rel_d1.B{ss-1}));
end

%same day correlation between A and B trials
for ss = 1:6
    %between A and B trials
    PVcorr_same_day.AB{ss} = corr(matchSTCs_nn.A{ss},matchSTCs_nn.B{ss}, 'type','Pearson','rows','complete');
end

%% Plot the mean PV correlation as a line plot

figure;
hold on
title('Population vector correlation relative to D1')
ylim([0 1])
ylabel('Mean correlation')
p1 = plot(meanPV_rel_d1.A,'b');
p2 = plot(meanPV_rel_d1.B,'r');
legend([p1 p2], 'A','B')
xticks(1:5)
xticklabels({'1 vs. 2', '1 vs.3', '1 vs. 6','1 vs. 7','1 vs. 8'})


%% %% PV correlation analysis across days (function of time)

%% OLD CODE
%{
%% Plot STCs from ROIs matching 1 and any chosen session thereafter
ses_comp = 7;

%get the matching ROIs from 2 ses
matching_ses_ROI_idxs = matching_list(:,[1,ses_comp]);
%get rid of nan values and get non-nan ROIs on both ses
matching_ses_ROI_idxs_nonan = matching_ses_ROI_idxs(find(sum(isnan(matching_ses_ROI_idxs),2)==0),:);

%generate matching A STC
A_comp_STCs = [A_STC{1}(:,matching_ses_ROI_idxs_nonan(:,1))', A_STC{ses_comp}(:,matching_ses_ROI_idxs_nonan(:,2))'];
%generate matching B STC
B_comp_STCs = [B_STC{1}(:,matching_ses_ROI_idxs_nonan(:,1))', B_STC{ses_comp}(:,matching_ses_ROI_idxs_nonan(:,2))'];

%sort by A on ses 1 and apply sort downstream
%sort by A trials on session 1
[~,matched_maxBin_A] = max(A_comp_STCs(:,1:100)', [], 1);
%[~,matched_maxBin_A] = max(A_comp_STCs(:,101:200)', [], 1);
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,matched_sortOrder_A] = sort(matched_maxBin_A,'ascend');
%find max bin for each ROI and return index
[~,matched_maxBin_B] = max(B_comp_STCs(:,1:100)', [], 1);
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,matched_sortOrder_B] = sort(matched_maxBin_B,'ascend');

%plot
figure('Position',[2630 50 490 890])
subplot(2,1,1)
imagesc(A_comp_STCs(matched_sortOrder_A,:))
hold on
colormap('jet')
caxis([0 1]);

subplot(2,1,2)
imagesc(B_comp_STCs(matched_sortOrder_B,:))
hold on
colormap('jet')
caxis([0 1]);

%}

end

