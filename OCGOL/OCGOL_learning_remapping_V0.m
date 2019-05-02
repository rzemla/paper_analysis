%% Import variables and define options

%run componenet registration across sessions
options.register = 0;

%input directories to matching function
path_dir = {'F:\OCGOL_training\I56_RLTS_041019\5A5B',...
    'F:\OCGOL_training\I56_RLTS_041019\ABrand_no_punish_041619'};

%cross session directory
crossdir = 'F:\OCGOL_training\I56_RLTS_041019\crossSession';

%load place cell variables for each session
%get mat directories in each output folder
for ii=1:size(path_dir,2)
    %get matfile names for each session
    matfiles{ii} = dir([path_dir{ii},'\output','\*.mat']);
end
%load in place cell variables (and others later)
for ii = 1:size(path_dir,2)
    %add event variables
    session_vars{ii} = load(fullfile(matfiles{ii}.folder,matfiles{ii}.name),'Place_cell', 'Behavior');
end

%% Match ROIs from across GOL days

if options.register == 1
    %run cross registration
    disp('Running registration of components');
    [registered] = match_ROIs_V2(path_dir,crossdir);
    
    %save registration variable in crosssession
    disp('Saving registered component matchings');
    save(fullfile(crossdir,'registered.mat'),'registered');
    
elseif options.register == 0
    %load the registered struct
    disp('Loading registered component matchings');
    load(fullfile(crossdir,'registered.mat'),'registered');
    
end

%% Visualize the matching ROIs that were matched above
%number of ROIs (rows) by sessions (cols)
rows = 20;
cols = 2; %take # of sessions as input
ROI_zooms = registered.multi.ROI_zooms;
ROI_outlines = registered.multi.ROI_outlines;

visualize_matches(rows,cols,ROI_zooms,ROI_outlines);


%% Plot smoothed event rate across track (function)

%% Plot dF/F rasters side by side by laps for matching ROIs
%all correct laps
figure;
for ii=1:size(registered.multi.assigned_all,1)
    
    subplot(2,3,1)
    imagesc(session_vars{1}.Place_cell{1, 3}.dF_lap_map_ROI{registered.multi.assigned_all(ii,1)})
    hold on;
    caxis([0 2])
    colormap(gca,'jet');
    hold off;
    
    subplot(2,3,2)
    
    
    subplot(2,3,3)
    imagesc(ROI_zooms{ii,1})
    hold on;
    colormap(gca, 'gray')
    xticks([])
    yticks([])
    b = bwboundaries(ROI_outlines{ii,1},'noholes');
    plot(b{1}(:,2),b{1}(:,1),'r')
    hold off
    
    
    subplot(2,3,4)
    imagesc(session_vars{2}.Place_cell{1, 3}.dF_lap_map_ROI{registered.multi.assigned_all(ii,2)})
    hold on;
    caxis([0 2])
    colormap(gca, 'jet');
    hold off;
    
    subplot(2,3,6)
    imagesc(ROI_zooms{ii,2})
    %imagesc(ROI_zooms{registered.multi.assigned_all(ii,2),2})
    hold on;
    colormap(gca, 'gray')
    xticks([])
    yticks([])
    %b = bwboundaries(ROI_outlines{registered.multi.assigned_all(ii,2),2},'noholes');
    b = bwboundaries(ROI_outlines{ii,2},'noholes');
    plot(b{1}(:,2),b{1}(:,1),'r')
    hold off
    
    pause;
    clf;
    
end


%% Learning pre and post event spiral
%display event spiral
%dF/F across laps between unlearned and learned day
%ROI outlines that show matching between days - need day from match ROI
%script



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

%make scatter plot
figure;
hold on
%ylim([0 0.8]);
%xlim([0.5 2.25]);
ylabel('Fraction of licks in reward zone')
xticks([0.7854, 1.1781,1.5708, 1.9635])
xticklabels({'\pi/4','3\pi/8','\pi/2','5\pi/8'})
xlabel('Mean distance to reward');
scatter(mean_Distance_TS_sig,frac_rew,'b');

%% Bar graph of distance to reward on 1st vs last day
%only tuned cells according to TS tuning criteria

figure;
hold on
ylabel('Distance to reward');
xticks([1 2 3 4 5 6 7]);
xticklabels({'RF Day 0','GOL 1 Day 1','GOL 1 Day 2','GOL 1 Day 3',...
    'GOL 2 Day 4','GOL 2 Day 5','GOL 2 Day 6'})
xtickangle(45)
yticks([0.7854,1.5708])
ylim([0.7854,2])
yticklabels({'\pi/4','\pi/2'})
bar(mean_Distance_TS_sig,'b')
refline(0,1.5708)





