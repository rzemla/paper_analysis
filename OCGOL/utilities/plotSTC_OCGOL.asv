function [outputArg1,outputArg2] = plotSTC_OCGOL(animal_data)

%% Import variables



%% Define tuned cells across trials

%spatial information criterion
Atuned = animal_data{1}.Place_cell{1}.Spatial_Info.significant_ROI == 1;
Btuned = animal_data{1}.Place_cell{2}.Spatial_Info.significant_ROI == 1;

AandB_tuned =  Atuned & Btuned;

%% Extract STCs with tuned ROIs

%definition of spatial tuning curve
%Gaussian smoothed onset rate map / spatial bin occupancy time (sec)
%Normalization from (0-1) for each ROI (ROI-by-ROI)

%these contain NaNs
A_STC_both = animal_data{1}.Place_cell{1}.Spatial_tuning_curve(:,AandB_tuned);
B_STC_both = animal_data{1}.Place_cell{2}.Spatial_tuning_curve(:,AandB_tuned);

%% Sort A STCs by maximum rate

%maxBin - spatial bin where activity is greatest for each ROI
[~,maxBin] = max(A_STC_both, [], 1);

%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,sortIdx] = sort(maxBin,'ascend');

%sort activity from both trials
A_STC_both_sorted = A_STC_both(:,sortIdx);
B_STC_both_sorted = B_STC_both(:,sortIdx);

% STC_sorted_nonan=STC_sorted(~any(isnan(STC_sorted),2),:);
% STC_dF_sorted_nonan=STC_dF_sorted(~any(isnan(STC_sorted),2),:);

find(isnan(animal_data{1}.Place_cell{1}.Spatial_tuning_curve) ==1)

%% Plot STCs side-by-side

figure;
subplot(1,2,1)
imagesc(A_STC_both_sorted')
hold on;
ylabel('Neuron #');
title('A trials')
colormap('jet');

subplot(1,2,2)
imagesc(B_STC_both_sorted')
hold on;
title('B trials')
colormap('jet')
colorbar


end

