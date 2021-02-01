function [si_ts_score_dist_data] = si_ts_score_distributions(path_dir)

%% Color scheme for figures
color_codes = [65,105,225; 220,20,60; 139, 0, 139; 128 128 128]./255;

%% Load the scores for each animal
for aa=1:size(path_dir,2)
    load_data_path{aa} = fullfile(path_dir{aa},'cumul_analysis','tuning_scores.mat');
    scores{aa} = load(string(load_data_path{aa}));
end

%% Create cumulative distribution and histograms for each class of neurons
%A,B,A&B,Neither

%first column - Alaps; second column - Blaps
for aa=1:size(path_dir,2)
    
    %SI
    %A only neurons
    si.Aonly{aa}(:,1) = scores{aa}.tuning_scores.si.Aonly.Alaps;
    si.Aonly{aa}(:,2) = scores{aa}.tuning_scores.si.Aonly.Blaps;
    
    %B only neurons
    si.Bonly{aa}(:,1) = scores{aa}.tuning_scores.si.Bonly.Alaps;
    si.Bonly{aa}(:,2) = scores{aa}.tuning_scores.si.Bonly.Blaps;
    
    %A&B neurons
    si.AB{aa}(:,1) = scores{aa}.tuning_scores.si.AB.Alaps;
    si.AB{aa}(:,2) = scores{aa}.tuning_scores.si.AB.Blaps;
    %neither neurons
    si.N{aa}(:,1) = scores{aa}.tuning_scores.si.neither.Alaps;
    si.N{aa}(:,2) = scores{aa}.tuning_scores.si.neither.Blaps;
    
    %TS
    %A only neurons
    ts.Aonly{aa}(:,1) = scores{aa}.tuning_scores.ts.Aonly.Alaps;
    ts.Aonly{aa}(:,2) = scores{aa}.tuning_scores.ts.Aonly.Blaps;
    
    %B only neurons
    ts.Bonly{aa}(:,1) = scores{aa}.tuning_scores.ts.Bonly.Alaps;
    ts.Bonly{aa}(:,2) = scores{aa}.tuning_scores.ts.Bonly.Blaps;
    
    %A&B neurons
    ts.AB{aa}(:,1) = scores{aa}.tuning_scores.ts.AB.Alaps;
    ts.AB{aa}(:,2) = scores{aa}.tuning_scores.ts.AB.Blaps;
    %neither neurons
    ts.N{aa}(:,1) = scores{aa}.tuning_scores.ts.neither.Alaps;
    ts.N{aa}(:,2) = scores{aa}.tuning_scores.ts.neither.Blaps;

end

%collapse into single matrix
%SI
si.Aonly_cumul = cell2mat(si.Aonly');
si.Bonly_cumul = cell2mat(si.Bonly');
si.AB_cumul = cell2mat(si.AB');
si.N_cumul = cell2mat(si.N');

%TS
ts.Aonly_cumul = cell2mat(ts.Aonly');
ts.Bonly_cumul = cell2mat(ts.Bonly');
ts.AB_cumul = cell2mat(ts.AB');
ts.N_cumul = cell2mat(ts.N');

%get mean for each animal
for aa=1:size(path_dir,2)
    %si
    mean_si_each.Aonly(aa,:) = nanmean(si.Aonly{aa},1);
    mean_si_each.Bonly(aa,:) = nanmean(si.Bonly{aa},1);
    mean_si_each.AB(aa,:) = nanmean(si.AB{aa},1);
    mean_si_each.N(aa,:) = nanmean(si.N{aa},1);
    
    %ts
    mean_ts_each.Aonly(aa,:) = nanmean(ts.Aonly{aa},1);
    mean_ts_each.Bonly(aa,:) = nanmean(ts.Bonly{aa},1);
    mean_ts_each.AB(aa,:) = nanmean(ts.AB{aa},1);
    mean_ts_each.N(aa,:) = nanmean(ts.N{aa},1);    
end

%% Return the values for plotting

si_ts_score_dist_data.si = si;
si_ts_score_dist_data.ts = ts;
si_ts_score_dist_data.mean_si_each = mean_si_each;
si_ts_score_dist_data.mean_ts_each = mean_ts_each;


%% Get histogram distributions around centerline - S.I.

%for si
options.xlims = [-0.2 0.2];
unity_hist_scatter_spatial_scores(si,options)

%for ts
options.xlims = [-1 1];
unity_hist_scatter_spatial_scores(ts,options)



end

