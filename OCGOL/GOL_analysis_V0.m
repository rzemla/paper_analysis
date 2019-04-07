%% Import variables and define options

%run componenet registration across sessions
options.register = 0;

%input directories to matching function
path_dir = {'G:\GOL\I55_RTLS_RF_GOL_022219\1',...
    'G:\GOL\I55_RTLS_RF_GOL_022219\2',...
    'G:\GOL\I55_RTLS_RF_GOL_022219\3',...
    'G:\GOL\I55_RTLS_RF_GOL_022219\4',...
    'G:\GOL\I55_RTLS_RF_GOL_022219\5',...
    'G:\GOL\I55_RTLS_RF_GOL_022219\6',...
    'G:\GOL\I55_RTLS_RF_GOL_022219\7'};

%cross session directory
crossdir = 'G:\GOL\I55_RTLS_RF_GOL_022219\crossSession';

%load place cell variables for each session
%get mat directories in each output folder
for ii=1:size(path_dir,2)
    %get matfile names for each session
    matfiles{ii} = dir([path_dir{ii},'\output','\*.mat']);
end
%load in place cell variables (and others later)
for ii = 1:size(path_dir,2)
    session_vars{ii} = load(fullfile(matfiles{ii}.folder,matfiles{ii}.name),'Place_cell', 'Behavior');
end

%% Match ROIs from across GOL days

if options.register == 1
    %run cross registration
    [registered] = match_ROIs_V2(path_dir,crossdir);
    
    %save registration variable in crosssession
    save(fullfile(crossdir,'registered.mat'),'registered');
    
elseif options.register == 0
    %load the registered struct
    load(fullfile(crossdir,'registered.mat'),'registered');
    
end

%% Plot smoothed event rate across track (function)



%% Calculate the distance from centroid to reward (function)

%compare to which reward (obsolete)
%if 1 --> GOL Block 1 location; if 2 --> GOL Block 2
options.rewardBlock = 1;

[centroids] = centroid_diff_reward(session_vars, options);


%% Plot mean centroid diff to reward for each session against performance (D vs. P) 

%fraction of licks in reward zone for each session
for ss=1:size(session_vars,2)
    frac_rew(ss) = session_vars{ss}.Behavior.performance.frac_rew;
end

%mean of centroid diff to reward of each session for TS sig neurons
for ss=1:size(session_vars,2)
    mean_Distance_TS_sig(ss) = nanmean(centroids.TS_sig.angleDiffReward{ss});
end

%plot scatterplot
distanceToRew_vs_performance_scatter(mean_Distance_TS_sig,frac_rew)

%% Plot performance across sessions

figure
hold on
title('Performance')
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
ylabel('Fraction of licks in reward zone');
yticks([0 0.3 0.7 1]);
ylim([0 1]);
plot(frac_rew, 'k', 'LineWidth',2)
xticks([1 2 3 4 5 6 7]);
xticklabels({'RF Day 0','GOL 1 Day 1','GOL 1 Day 2','GOL 1 Day 3',...
    'GOL 2 Day 4','GOL 2 Day 5','GOL 2 Day 6'})
xtickangle(45);
%reward zone switch line
plot([4.5 4.5],[0 1], '--b','LineWidth',2);

%% Plot D0 dF/F raster vs. D3 GOL 1 raster

figure;
subplot(1,2,1)
imagesc(session_vars{1, 1}.Place_cell{1, 1}.Spatial_tuning_dF_sorted)
hold on
plot([30 30],[536 1],'k','LineWidth',2)
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 0.5;
title('Day 0 RF')
colormap('jet')
caxis([0 1]) 
ylabel('Neuron #');
xlabel('Position')
xticks([1 100]);
xticklabels({'0', '1'})
hold off


subplot(1,2,2)
imagesc(session_vars{1, 4}.Place_cell{1, 1}.Spatial_tuning_dF_sorted)
hold on
plot([30 30],[406 1],'k','LineWidth',2)
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1;
title('Day 3 GOL 1')
caxis([0 1]) 
xlabel('Position')
xticks([1 100]);
xticklabels({'0', '1'})

%% Bar graph of distance to reward on 1st vs last day
%only tuned cells according to TS tuning criteria

plotBarMeanDistanceToReward(mean_Distance_TS_sig)

%% Compute the width of each cell

%use Dombeck code to detect place fields
[Place_cell_field] =  placeFieldDetection_barthosBeta(session_vars{1, 1}.Place_cell{1, 1} );

%plot dF spatial tuning curve and associated place field

figure;
for ii=1:size(Place_cell_field.Tuned_ROI_mask,2)
    hold on
    plot(Place_cell_field.Spatial_tuning_dF(:,ii),'k');
    %binary field definition
    plot(Place_cell_field.placeField.binary(:,ii),'r');
    hold off
    
    pause;
    clf
end

%Danielson
%(for cells with multiple fields
%the field with the highest in-field transient rate was used to calculate the centroid)
% The place field sensitivity (defined as the fraction of complete forward passes through the
% place field associated with a significant Ca2+ transient) was significantly higher in deep than in superficial (n = 14
% mice, p = 0.001, paired T-Test), but the (ii) place cell specificity (defined as the fraction of running-related transients
% occurring within the place field)



