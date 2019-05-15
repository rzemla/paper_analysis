function [outputArg1,outputArg2] = plotSTC_OCGOL_training(animal_data, tunedLogical)

%% Import variables

%% Define tuned cells across trials

%spatial information criterion
Atuned = tunedLogical.si.Atuned;
Btuned = tunedLogical.si.Btuned;

AandB_tuned =  tunedLogical.si.AandB_tuned;
AorB_tuned = tunedLogical.si.AorB_tuned;
onlyA_tuned = tunedLogical.si.onlyA_tuned;
onlyB_tuned = tunedLogical.si.onlyB_tuned;

%% Extract mean dF map in each spatial bin (not normalized and not occupancy divided) (100 bins)
A_df = animal_data{1}.Place_cell{1}.Spatial_Info.mean_dF_map{8}; 
B_df = animal_data{1}.Place_cell{2}.Spatial_Info.mean_dF_map{8}; 

A_df_both = A_df(:,AandB_tuned);
B_df_both = B_df(:,AandB_tuned);

A_df_onlyA = A_df(:,onlyA_tuned);
B_df_onlyA = B_df(:,onlyA_tuned);


%% Sort A dF/F STCs by maximum rate

%maxBin - spatial bin where activity is greatest for each ROI
[~,maxBin] = max(A_df_both, [], 1);

%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,sortOrderBoth] = sort(maxBin,'ascend');

clear maxBin

[~,maxBin] = max(A_df_onlyA, [], 1);
[~,sortOrderA_only] = sort(maxBin,'ascend');

%sort activity from both trials
% A_STC_both_sorted = A_STC_both(:,sortOrder);
% B_STC_both_sorted = B_STC_both(:,sortOrder);
% 
% A_STC_onlyA_sorted = A_STC_onlyA(:,sortOrder);
% B_STC_onlyA_sorted = B_STC_onlyA(:,sortOrder);

A_df_both_sorted = A_df_both(:,sortOrderBoth);
B_df_both_sorted = B_df_both(:,sortOrderBoth);

A_df_onlyA_sorted = A_df_onlyA(:,sortOrderA_only);
B_df_onlyA_sorted = B_df_onlyA(:,sortOrderA_only);

% STC_sorted_nonan=STC_sorted(~any(isnan(STC_sorted),2),:);
% STC_dF_sorted_nonan=STC_dF_sorted(~any(isnan(STC_sorted),2),:);

%find(isnan(animal_data{1}.Place_cell{1}.Spatial_tuning_curve) ==1)

%% Plot STCs side-by-side

% figure;
% subplot(1,2,1)
% imagesc(A_STC_onlyA_sorted')
% hold on;
% ylabel('Neuron #');
% title('A trials')
% colormap('jet');
% 
% subplot(1,2,2)
% imagesc(B_STC_onlyA_sorted')
% hold on;
% title('B trials')
% colormap('jet')
% %colorbar

%% Plot dF/F maps side-by-side

figure;
subplot(1,2,1)
imagesc(A_df_both_sorted')
hold on;
ylabel('Neuron #');
title('A trials')
caxis([0 2])
colormap('jet');
colorbar

subplot(1,2,2)
imagesc(B_df_both_sorted')
hold on;
title('B trials')
caxis([0 2])
colormap('jet')
colorbar


%% Plot dF/F maps side-by-side

figure;
subplot(1,2,1)
imagesc(A_df_onlyA_sorted')
hold on;
ylabel('Neuron #');
title('A trials - only A tuned')
caxis([0 2])
colormap('jet');
colorbar

subplot(1,2,2)
imagesc(B_df_onlyA_sorted')
hold on;
title('B trials - only A tuned')
caxis([0 2])
colormap('jet')
colorbar


%% Extract STCs with tuned ROIs - in nontuned neurons will scale the weakest signal to 1 regardless of tuning

%definition of spatial tuning curve
%Gaussian smoothed onset rate map / spatial bin occupancy time (sec)
%Normalization from (0-1) for each ROI (ROI-by-ROI)

%these contain NaNs
STC_A = animal_data{1}.Place_cell{1}.Spatial_tuning_curve;
STC_B = animal_data{1}.Place_cell{2}.Spatial_tuning_curve;

A_STC_both = STC_A(:,AandB_tuned);
B_STC_both = STC_B(:,AandB_tuned);

A_STC_onlyA = STC_A(:,onlyA_tuned);
B_STC_onlyA = STC_B(:,onlyA_tuned);



end

