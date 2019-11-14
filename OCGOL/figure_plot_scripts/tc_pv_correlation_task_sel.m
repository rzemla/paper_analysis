function [outputArg1,outputArg2] = tc_pv_correlation_task_sel(path_dir)

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

%% Tuning curve (TC) correlation

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


end
