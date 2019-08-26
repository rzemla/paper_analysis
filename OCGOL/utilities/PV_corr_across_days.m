function [outputArg1,outputArg2] = PV_corr_across_days(animal_data, tunedLogical,registered,options)


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
    A_STC{ss} = animal_data{ss}.Place_cell{1}.Spatial_tuning_curve;
    B_STC{ss} = animal_data{ss}.Place_cell{2}.Spatial_tuning_curve;
    
    %A_STC{ss} =
    %B_STC{ss} = 
    
end

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
matchSTCs_A = zeros(size(matching_list,1),7*100);
matchSTCs_B = zeros(size(matching_list,1),7*100);

%fill the nan assignments to STC to be NaN
for ss=1:7
   nan_log = isnan(matching_list(:,ss));
   matchSTCs_A(nan_log,1+(ss-1)*100:ss*100) = NaN;
   matchSTCs_B(nan_log,1+(ss-1)*100:ss*100) = NaN;
end

%fill the matching assignments with trial-normalized STC
for ss=1:7
   nan_log = isnan(matching_list(:,ss));
   matchSTCs_A(~nan_log,1+(ss-1)*100:ss*100) = A_STC{ss}(:,matching_list(~nan_log,ss))';
   matchSTCs_B(~nan_log,1+(ss-1)*100:ss*100) = B_STC{ss}(:,matching_list(~nan_log,ss))';
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

%sorted matrix by day 1
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


end
