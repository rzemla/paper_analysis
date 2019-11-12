%% Define experiment core directories
%screened and good datasets for analysis
%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_120618_1'}; %
%field rate error (after PSEM silencing
path_dir{1} = {'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1'};
path_dir{2} = {'G:\Figure_2_3_selective_remap\I42R_AB_d1_032118_1'};
path_dir{3} = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1'};
path_dir{4} = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_2'};
path_dir{5} = {'G:\Figure_2_3_selective_remap\I53LT_AB_sal_113018_1'}; %place field finder problem - adjust
path_dir{6} = {'G:\Figure_2_3_selective_remap\I56_RTLS_AB_prePost_sal_042419_1'}; %place field finder problem - adjust
path_dir{7} = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_113018_1'};
path_dir{8} = {'G:\Figure_2_3_selective_remap\I57_RTLS_AB_prePost_792_042519_1'};
path_dir{9} = {'G:\Figure_2_3_selective_remap\I45_RT_AB_d1_062018_1'};
path_dir{10} = {'G:\Figure_2_3_selective_remap\I46_AB_d1_062018_1'};
path_dir{11} = {'G:\Figure_2_3_selective_remap\I57_LT_ABrand_no_punish_042119_1'};

%N =10; FOV = 11 (2 FOV from single animal)

%% Fraction of place cell (A sel, B sel, A&B remapping,neither) - bar chart and pie chart
%Figure 2C
fraction_place_cells(path_dir)


%% Remapping centroids
%Figure 3E/F
options.lowPVcorr = [6 7 8];
%organize this 
remapping_centroids(path_dir,options)

%% Centroid difference for A&B tuned neurons and centroid diff as fxn of max bin 
%scatterplot of centroid difference as a function of center between
%centroid of max place field
centroid_difference(path_dir)

%% Figure 3C - fractional distribution of remapping neuron subtypes

frac_remapping_neurons(path_dir)

%% AUC/min scatterplots of A vs B neurons for each animal 
%read relevant data
for ee=1:size(path_dir,2)
    load_data_path_auc{ee} = fullfile(path_dir{ee},'cumul_analysis','auc.mat');
    auc_data{ee} = load(string(load_data_path_auc{ee}));
end

%add colormap to this with cbrewer
%blues (A trials)
cmap_blue = cbrewer('seq','Blues',size(path_dir,2)+2);
%reds (B trials)
cmap_red = cbrewer('seq','Reds',size(path_dir,2)+2);

figure;
hold on
axis square
xlim([0 10])
ylim([0 10])
xticks(0:2:10)
yticks(0:2:10)
set(gca,'FontSize',14)
set(gca,'LineWidth',2)
xlabel('AUC/min - A trials')
ylabel('AUC/min - B trials')
title('Activity rate for task selective place cells')
%A selective
for ee=1:size(path_dir,2)
    scatter(mean(auc_data{ee}.total_AUC_min.A(1,:)),mean(auc_data{ee}.total_AUC_min.A(2,:)),'filled','MarkerFaceColor',cmap_blue(ee+2,:))
end
%B selective
for ee=1:size(path_dir,2)
    scatter(mean(auc_data{ee}.total_AUC_min.B(1,:)),mean(auc_data{ee}.total_AUC_min.B(2,:)),'filled','MarkerFaceColor',cmap_red(ee+2,:))
end
%plot center line
plot([0 10], [0 10],'Color',[0.5 0.5 0.5],'LineStyle', '--','LineWidth',2)

%% Centroid distribution

%read relevant data
for ee=1:size(path_dir,2)
    load_data_path_cent{ee} = fullfile(path_dir{ee},'cumul_analysis','centroid.mat');
    centroid_data{ee} = load(string(load_data_path_cent{ee}));
end

%create matrices with centroid counts for each subclass of tuned neuron
for ee=1:size(path_dir,2)
    centroid_mat_A(ee,:) = centroid_data{ee}.centroid_ct.A;
    centroid_mat_B(ee,:) = centroid_data{ee}.centroid_ct.B;
end

%normalize to itself and get mean for A and B
centroid_mat_A_norm = centroid_mat_A./sum(centroid_mat_A,2);
centroid_mat_B_norm = centroid_mat_B./sum(centroid_mat_B,2);

%take mean in each bin from each animal
centroid_mat_A_norm_mean = mean(centroid_mat_A_norm,1);
centroid_mat_B_norm_mean = mean(centroid_mat_B_norm,1);

%get normalized standard deviation for centroid of task-selective neurons
centroid_mat_A_norm_std = std(centroid_mat_A_norm,0,1);
centroid_mat_B_norm_std = std(centroid_mat_B_norm,0,1);

%plot as fraction of neurons tuned in each bin with shaded std
figure('Position', [2015 120 635 820]);
subplot(2,1,1)
hold on;
set(gca,'FontSize',14)
set(gca,'LineWidth',3)
title('A-selective')
ylabel('Normalized density')
ylim([-0.03 0.13])
plot(centroid_mat_A_norm_mean,'b')
s = shadedErrorBar(1:100,centroid_mat_A_norm,{@mean,@std},'lineprops','-','transparent',true,'patchSaturation',0.20);
set(s.edge,'LineWidth',1.5,'LineStyle','-','Color',[[65,105,225]/255, 0.2]) %last value add transparency value
s.mainLine.LineWidth = 2;
s.mainLine.Color = [65,105,225]/255;
s.patch.FaceColor = [65,105,225]/255;
hold off

subplot(2,1,2)
hold on;
set(gca,'FontSize',14)
set(gca,'LineWidth',3)
title('B-selective')
ylabel('Normalized density')
ylim([-0.03 0.13])
plot(centroid_mat_B_norm_mean,'r')
s = shadedErrorBar(1:100,centroid_mat_B_norm,{@mean,@std},'lineprops','-','transparent',true,'patchSaturation',0.20);
set(s.edge,'LineWidth',1.5,'LineStyle','-','Color',[[220,20,60]/255, 0.2]) %last value add transparency value
s.mainLine.LineWidth = 2;
s.mainLine.Color = [220,20,60]/255;
s.patch.FaceColor = [220,20,60]/255;
xlabel('Spatial bin')

%plot this an empirical distribution curve (histogram)
figure;
subplot(2,1,1)
hold on
ylim([0 0.03])
histogram('BinEdges',0.5:1:100.5,'BinCounts',sum(centroid_mat_A,1),'Normalization', 'probability')
subplot(2,1,2)
hold on
ylim([0 0.03])
histogram('BinEdges',0.5:1:100.5,'BinCounts',sum(centroid_mat_B,1),'Normalization', 'probability')

%plot in polar space (cumulative)
figure;
pax1 = subplot(1,2,1,polaraxes);
hold on
title('A selective')
pax1.ThetaAxisUnits = 'degrees';
pax1.RAxisLocation = 45;
pax1.RLim = [0 0.04];
pax1.RColor = 'r';
pax1.FontSize = 14;
polarhistogram(pax1,'BinEdges',0:(2*pi)/100:2*pi,'BinCounts', sum(centroid_mat_A,1),'Normalization','probability',...
    'FaceColor','blue')

pax2 = subplot(1,2,2,polaraxes);
hold on
title('B selective')
pax2.ThetaAxisUnits = 'degrees';
pax2.RAxisLocation = 45;
pax2.RLim = [0 0.04];
pax2.RColor = 'r';
pax2.FontSize = 14;
polarhistogram(pax2,'BinEdges',0:(2*pi)/100:2*pi,'BinCounts', sum(centroid_mat_B,1),'Normalization', 'probability',...
    'FaceColor','red')

%% Polar plot mean
figure;
pax1 = subplot(1,2,1,polaraxes);
hold on
title('A selective')
pax1.ThetaAxisUnits = 'degrees';
pax1.RAxisLocation = 45;
pax1.RLim = [0 0.06];
pax1.RColor = 'r';
pax1.FontSize = 14;
polarhistogram(pax1,'BinEdges',0:(2*pi)/100:2*pi,'BinCounts', centroid_mat_A_norm_mean,'Normalization','probability',...
    'FaceColor','blue')

pax2 = subplot(1,2,2,polaraxes);
hold on
title('B selective')
pax2.ThetaAxisUnits = 'degrees';
pax2.RAxisLocation = 45;
pax2.RLim = [0 0.06];
pax2.RColor = 'r';
pax2.FontSize = 14;
polarhistogram(pax2,'BinEdges',0:(2*pi)/100:2*pi,'BinCounts', centroid_mat_B_norm_mean,'Normalization', 'probability',...
    'FaceColor','red')

%empirical cdf
if 0
    figure;
    hold on
    ecdf((1:100),'Frequency',sum(centroid_mat_A,1)/sum(sum(centroid_mat_A,1)))
    ecdf((1:100),'Frequency',sum(centroid_mat_B,1)/sum(sum(centroid_mat_B,1)))
end

% figure;
% hold on
% subplot(2,1,1)
% plot(sum(centroid_mat_A,1),'b')
% subplot(2,1,2)
% plot(sum(centroid_mat_B,1),'r')



%% PV/TC correlation
%read relevant data
for ee=1:size(path_dir,2)
    load_data_path_corr{ee} = fullfile(path_dir{ee},'cumul_analysis','corr.mat');
    correlation_data{ee} = load(string(load_data_path_corr{ee}));
end

%get diagonal of each PV correlation matrix
for ee=1:size(path_dir,2)
    diag_PVcorr_mat(ee,:) = diag(correlation_data{ee}.correlation.PVcorr);
end
%get mean of PV correlation
mean_diag_PVcorr = mean(diag_PVcorr_mat,1);
%get std at each spatial bin
std_diag_PVcorr = std(diag_PVcorr_mat,0,1);
%get sem at each spatial bin
sem_diag_PVcorr = std_diag_PVcorr./sqrt(size(diag_PVcorr_mat,1));

%plot the pv correlation (each animal) across track length and assn sem at
%each bin (around mean)
figure('Position',[2050 530 1380 420]);
hold on
ylim([-0.05 1])
yticks(-0.2:0.2:1)
xlim([1,100])
xticks(0:10:100)
ylabel('Correlation coef');
xlabel('Spatial bin');
set(gca,'FontSize',14);
title('Population vector correlation')
%plot each correlation traces from each animal
for ee=1:size(path_dir,2)
    plot(diag_PVcorr_mat(ee,:),'Color',[0.7 0.7 0.7],'LineWidth',1','LineStyle','-');
end
%plot the mean PV correlation from each bin
%plot(mean_diag_PVcorr,'Color','k','LineWidth',1.5)
%plot std around the PV correlation
%use shaded plot for this
%plot upper std line
s = shadedErrorBar(1:100,diag_PVcorr_mat,{@mean,@std},'lineprops','-','transparent',true,'patchSaturation',0.20);
set(s.edge,'LineWidth',1.5,'LineStyle','-','Color',[[0 153 0]/255, 0.2]) %last value add transparency value
s.mainLine.LineWidth = 2;
s.mainLine.Color = [0 153 0]/255;
s.patch.FaceColor = [0 153 0]/255;
%behaviorally relevant lines
plot([30 30],[-0.05 1],'r--','LineWidth',1.5);
text([31 31],[0.9 0.9],'Reward zone B','Color','r','FontSize',14)
plot([70 70],[-0.05 1],'b--','LineWidth',1.5);
text([71 71],[0.9 0.9],'Reward zone A','Color','b','FontSize',14)
plot([10 10],[-0.05 1],'--','Color',[0 153 153]/255,'LineWidth',1.5);
text([11 11],[0.9 0.9],'Odor zone end','Color',[0 153 153]/255,'FontSize',14)

%TC
%place in separate cells (combine diagnonal values)
for ee=1:size(path_dir,2)
    %A tuned by stat sig TS only
    diag_TC.A{ee} = diag(correlation_data{ee}.correlation.TCcorr.Aonly);
    %B tuned by stat sig TS only
    diag_TC.B{ee} = diag(correlation_data{ee}.correlation.TCcorr.Bonly);
    %A selective (with added filters)
    diag_TC.Asel{ee} = diag(correlation_data{ee}.correlation.TCcorr.Aselective);
    %B selective (with added filters)
    diag_TC.Bsel{ee} = diag(correlation_data{ee}.correlation.TCcorr.Bselective);
    %current A&B tuned without filter
    diag_TC.AB{ee} = diag(correlation_data{ee}.correlation.TCcorr.AB);
    %all neurons
    diag_TC.all{ee} = diag(correlation_data{ee}.correlation.TCcorr.all);
end

%combine correlation values into single matrix
comb_TC.A = cell2mat(diag_TC.A');
comb_TC.B = cell2mat(diag_TC.B');
comb_TC.AB = cell2mat(diag_TC.AB');
comb_TC.Asel = cell2mat(diag_TC.Asel');
comb_TC.Bsel = cell2mat(diag_TC.Bsel');

%plot the mean TC for each animal with bar and SD with A-selective and
%B-selective and for (current, not-filtered) AB neurons 
for ee=1:size(path_dir,2)
    mean_TC.Asel(ee) = nanmean(diag_TC.Asel{ee});
    mean_TC.Bsel(ee) = nanmean(diag_TC.Bsel{ee});
    mean_TC.AB(ee) = nanmean(diag_TC.AB{ee}); 
    mean_TC.all(ee) = nanmean(diag_TC.all{ee}); 
end

%plot the mean TC boxplots for A-sel,B-sel and AB

figure('Position',[650 190 500 420])
hold on
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)
ylim([0 0.8])
title('Tuning curve correlation between task laps')
boxplot([mean_TC.Asel;mean_TC.Bsel;mean_TC.AB;mean_TC.all]','Colors',[65,105,225; 220,20,60; 255,0,255; 0 0 0]/255)
xticklabels({'A-selective','B-selective','A&B','All'})
xtickangle(45)
ylabel('Correlation coef.')

%plot combined TC corr values as single plot
figure;
hold on;
c1 = cdfplot(comb_TC.A);
c1.Color = 'b';
c1.LineWidth = 2;
c2 = cdfplot(comb_TC.B);
c2.Color = 'r';
c2.LineWidth = 2;
c3 = cdfplot(comb_TC.AB);
c3.Color = 'm';
c3.LineWidth = 2;
title('Combined tuning curve corr. CDF')
xlabel('Correlation coef.');
ylabel('Cumulative fraction');
grid off
set(gca,'FontSize',14)
xlim([-0.5 1]);
xticks([-0.5 0 0.5 1])
legend([c1 c2 c3],{'A','B','A&B'},'Location','northwest')

%plot this as ecdf curve fraction vs correlation (individual)
figure('Position', [2060 550 1400 420]);
subplot(1,3,1)
hold on
for ee=1:size(path_dir,2)
    %A selective
    c1 = cdfplot(diag(correlation_data{ee}.correlation.TCcorr.Aonly));
    c1.Color = 'b';
    c1.LineWidth = 1.5;
end
title('A-selective tuning curve corr. CDF')
xlabel('Correlation coef.');
ylabel('Cumulative fraction');
grid off
xlim([-0.5 1]);
xticks([-0.5 0 0.5 1])
axis square
set(gca,'FontSize', 12)
set(gca,'LineWidth',2)

subplot(1,3,2)
hold on
for ee=1:size(path_dir,2)
    %B-selective
    c2 = cdfplot(diag(correlation_data{ee}.correlation.TCcorr.Bonly));
    c2.Color = 'r';
    c2.LineWidth = 1.5;
end
title('B-selective tuning curve corr. CDF')
xlabel('Correlation coef.');
ylabel('');
grid off
xlim([-0.5 1]);
xticks([-0.5 0 0.5 1])
axis square
set(gca,'FontSize', 12)
set(gca,'LineWidth',2)

subplot(1,3,3)
hold on
for ee=1:size(path_dir,2)
    %A&B
    c3 = cdfplot(diag(correlation_data{ee}.correlation.TCcorr.AB));
    c3.Color = 'm';
    c3.LineWidth = 1.5;
end
title('A&B tuning curve corr. CDF')
xlabel('Correlation coef.');
ylabel('');
grid off
xlim([-0.5 1]);
xticks([-0.5 0 0.5 1])
axis square
set(gca,'FontSize', 12)
set(gca,'LineWidth',2)
%% Place field analysis (width and number for selective neurons

%read relevant data
for ee=1:size(path_dir,2)
    load_data_path_pf{ee} = fullfile(path_dir{ee},'cumul_analysis','placeField_dist.mat');
    placeField_data{ee} = load(string(load_data_path_pf{ee}));
end

%combine field counts for Asel and Bsel into 1 matrix
%create matrices with centroid counts for each subclass of tuned neuron
for ee=1:size(path_dir,2)
    %pf_count_mat_A(ee,:) = placeField_data{ee}.placeField_dist.field_count_A;
    pf_count_mat_A(ee,:) = placeField_data{ee}.placeField_dist.task_sel.A.field_count;
    % placeField_data{ee}.placeField_dist.all.A.ts
    %pf_count_mat_B(ee,:) = placeField_data{ee}.placeField_dist.field_count_B;
pf_count_mat_B(ee,:) = placeField_data{ee}.placeField_dist.task_sel.B.field_count;
end

%normalize as fraction of neurons for each animal/exp for A-sel/B-sel
pf_count_mat_A_norm = pf_count_mat_A./sum(pf_count_mat_A,2);
pf_count_mat_B_norm = pf_count_mat_B./sum(pf_count_mat_B,2);

%get means for each subclass
mean_pf_norm_A = mean(pf_count_mat_A_norm,1);
mean_pf_norm_B = mean(pf_count_mat_B_norm,1);

%get std and sem for each group
std_pf_norm_A = std(pf_count_mat_A_norm,0,1);
std_pf_norm_B = std(pf_count_mat_B_norm,0,1);
%sem
sem_pf_norm_A = std_pf_norm_A./sqrt(size(pf_count_mat_A,1));
sem_pf_norm_B = std_pf_norm_B./sqrt(size(pf_count_mat_B,1));
%grouped sem
grouped_norm_sem = [sem_pf_norm_A',sem_pf_norm_B'];

%combined means from norm counts
grouped_norm_mean = [mean_pf_norm_A', mean_pf_norm_B'];

%sum A and B
grouped_pf_counts = [sum(pf_count_mat_A,1)',sum(pf_count_mat_B,1)'];
%normalized for each group
grouped_pf_counts_norm = [(sum(pf_count_mat_A,1)./sum(sum(pf_count_mat_A,1)))',...
            (sum(pf_count_mat_B,1)./sum(sum(pf_count_mat_B,1)))'];
        
      
%Place field analysis plotting
%plot bar
figure;
hold on;
title('Place fields per neuron - S.I.');
%bar the mean for each group
b = bar(1:3,grouped_norm_mean,'FaceColor', 'flat');
pause(0.1)
%plot the sem for each mean for each group
for ib = 1:numel(b)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    if ib ==1
        xData = b(ib).XData + b(ib).XOffset;
    elseif ib ==2
        xData = b(ib).XData + b(ib).XOffset;
    end
    errorbar(xData,grouped_norm_mean(:,ib)',grouped_norm_sem(:,ib),'k.')
end

%set A group bars to blue
b(1).CData(1:3,:) =  repmat([0 0 1],3,1);
%set B group bars to red
b(2).CData(1:3,:) =  repmat([1 0 0],3,1);
xticks([1 2 3]);
xticklabels({'1','2','3+'});
ylabel('Fraction of neurons');
legend('A','B')

