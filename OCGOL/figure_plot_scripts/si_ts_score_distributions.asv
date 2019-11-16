function [outputArg1,outputArg2] = si_ts_score_distributions(path_dir)


%% Load the scores for each animal
for aa=1:size(path_dir,2)
    load_data_path{aa} = fullfile(path_dir{aa},'cumul_analysis','tuning_scores.mat');
    scores{aa} = load(string(load_data_path{aa}));
end

%% Create cumulative distribution and histograms for each class of neurons
%A,B,A&B,Neither

%first column - Alaps; second column - Blaps
for aa=1:size(path_dir,2)
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

end

%collapse into single matrix
si.Aonly_cumul = cell2mat(si.Aonly');
si.Bonly_cumul = cell2mat(si.Bonly');
si.AB_cumul = cell2mat(si.AB');
si.N_cumul = cell2mat(si.N');

%try scatterplot
figure
hold on
axis square
xlim([0 0.25])
ylim([0 0.24])
scatter(si.Aonly_cumul(:,1),si.Aonly_cumul(:,2),'MarkerFaceColor','b')
scatter(si.Bonly_cumul(:,1),si.Bonly_cumul(:,2),'MarkerFaceColor','r')
scatter(si.AB_cumul(:,1),si.AB_cumul(:,2),'MarkerFaceColor','m')
scatter(si.N_cumul(:,1),si.N_cumul(:,2),'MarkerFaceColor',[1 1 1]*0.5)



end

