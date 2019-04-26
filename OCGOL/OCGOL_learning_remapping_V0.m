%% Import variables and define options

%run componenet registration across sessions
options.register = 1;

%input directories to matching function
path_dir = {'G:\OCGOL_training\I56_RLTS_041019\5A5B',...
    'G:\OCGOL_training\I56_RLTS_041019\ABrand_no_punish_041619'};

%cross session directory
crossdir = 'G:\OCGOL_training\I56_RLTS_041019\crossSession';

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





