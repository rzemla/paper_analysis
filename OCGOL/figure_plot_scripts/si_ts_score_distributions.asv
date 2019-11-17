function [outputArg1,outputArg2] = si_ts_score_distributions(path_dir)

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

%% Scatterplots of mean scores
%% Plot
figure('Position',[2295 322 1080 420]);
subplot(1,2,1)
hold on
axis square
xlim([0 0.25])
ylim([0 0.25])
xticks([0 0.1 0.2])
yticks([0 0.1 0.2])
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
xlabel('A trials')
ylabel('B trials')
title('S.I. score')

%plot connecting lines for same animal points
% for ee=1:size(path_dir,2)
%     plot([mean_AUC.A(ee,1),mean_AUC.B(ee,1)],[mean_AUC.A(ee,2),mean_AUC.B(ee,2)] ,'Color',[1 1 1]*0.5,'LineWidth',0.5)
% end

%plot center line
plot([0 10], [0 10],'Color',[1 1 1]*0,'LineStyle', '--','LineWidth',2)

s1 = scatter(mean_si_each.Aonly(:,1),mean_si_each.Aonly(:,2),'filled','MarkerFaceColor',color_codes(1,:));
s2 = scatter(mean_si_each.Bonly(:,1),mean_si_each.Bonly(:,2),'filled','MarkerFaceColor',color_codes(2,:));
s4 = scatter(mean_si_each.N(:,1),mean_si_each.N(:,2),'filled','MarkerFaceColor',color_codes(4,:));
s3 = scatter(mean_si_each.AB(:,1),mean_si_each.AB(:,2),'filled','MarkerFaceColor',color_codes(3,:));
legend([s1 s2 s3 s4],{'A','B','A&B','Neither'},'location','southeast')

subplot(1,2,2)
hold on
axis square
xlim([0 1])
ylim([0 1])
xticks([0 0.5 1])
yticks([0 0.5 1])
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
xlabel('A trials')
ylabel('B trials')
title('T.S. score')

%plot connecting lines for same animal points
% for ee=1:size(path_dir,2)
%     plot([mean_AUC.A(ee,1),mean_AUC.B(ee,1)],[mean_AUC.A(ee,2),mean_AUC.B(ee,2)] ,'Color',[1 1 1]*0.5,'LineWidth',0.5)
% end

%plot center line
plot([0 10], [0 10],'Color',[1 1 1]*0,'LineStyle', '--','LineWidth',2)

s1 = scatter(mean_ts_each.Aonly(:,1),mean_ts_each.Aonly(:,2),'filled','MarkerFaceColor',color_codes(1,:),'MarkerEdgeColor','none');
s2 = scatter(mean_ts_each.Bonly(:,1),mean_ts_each.Bonly(:,2),'filled','MarkerFaceColor',color_codes(2,:),'MarkerEdgeColor','none');
s4 = scatter(mean_ts_each.N(:,1),mean_ts_each.N(:,2),'filled','MarkerFaceColor',color_codes(4,:),'MarkerEdgeColor','none');
s3 = scatter(mean_ts_each.AB(:,1),mean_ts_each.AB(:,2),'filled','MarkerFaceColor',color_codes(3,:),'MarkerEdgeColor','none');
legend([s1 s2 s3 s4],{'A','B','A&B','Neither'},'location','southeast')

%% Get histogram distributions around centerline - S.I.

%for si
options.xlims = [-0.2 0.2];
unity_hist_scatter_spatial_scores(si,options)

%for ts
options.xlims = [-1 1];
unity_hist_scatter_spatial_scores(ts,options)

%% Plot scatter plots of all neurons for each category
marker_size = 5;

%try scatterplot
figure('Position',[2208 244 1198 512])
%si
subplot(1,2,1)
hold on
title('S.I. score')
axis square
xlim([-0.01 0.25])
ylim([-0.01 0.25])
xticks([0 0.1 0.2])
yticks([0 0.1 0.2])
xlabel('A laps')
ylabel('B laps')
scatter(si.N_cumul(:,1),si.N_cumul(:,2),marker_size,'MarkerFaceColor',color_codes(4,:),'MarkerEdgeColor',color_codes(4,:))
scatter(si.Aonly_cumul(:,1),si.Aonly_cumul(:,2),marker_size,'MarkerFaceColor',color_codes(1,:),'MarkerEdgeColor',color_codes(1,:))
scatter(si.Bonly_cumul(:,1),si.Bonly_cumul(:,2),marker_size,'MarkerFaceColor',color_codes(2,:),'MarkerEdgeColor',color_codes(2,:))
scatter(si.AB_cumul(:,1),si.AB_cumul(:,2),marker_size,'MarkerFaceColor',color_codes(3,:),'MarkerEdgeColor',color_codes(3,:))
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%center line
plot([-0.01 0.25],[-0.01 0.25],'k--')

%ts
subplot(1,2,2)
hold on
title('T.S. score')
axis square
xlim([0 1.1])
ylim([0 1.1])
xticks([0 0.5 1])
yticks([0 0.5 1])
xlabel('A laps')
ylabel('B laps')
scatter(ts.N_cumul(:,1),ts.N_cumul(:,2),marker_size,'MarkerFaceColor',color_codes(4,:),'MarkerEdgeColor',color_codes(4,:))
scatter(ts.Aonly_cumul(:,1),ts.Aonly_cumul(:,2),marker_size,'MarkerFaceColor',color_codes(1,:),'MarkerEdgeColor',color_codes(1,:))
scatter(ts.Bonly_cumul(:,1),ts.Bonly_cumul(:,2),marker_size,'MarkerFaceColor',color_codes(2,:),'MarkerEdgeColor',color_codes(2,:))
scatter(ts.AB_cumul(:,1),ts.AB_cumul(:,2),marker_size,'MarkerFaceColor',color_codes(3,:),'MarkerEdgeColor',color_codes(3,:))
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%center line
plot([0 1.1],[0 1.1],'k--')

end

