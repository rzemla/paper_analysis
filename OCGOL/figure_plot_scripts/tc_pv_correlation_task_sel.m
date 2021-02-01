function [tc_corr_sel_data] = tc_pv_correlation_task_sel(path_dir)

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

%% plot the pv correlation (each animal) across track length and assn sem at
%each bin (around mean)
if 0
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
end

%% Tuning curve (TC) correlation (export this data and related curve!!!)

%place in separate cells (combine diagnonal values)
for ee=1:size(path_dir,2)
    %A tuned by stat sig TS only
    %diag_TC.A{ee} = diag(correlation_data{ee}.correlation.TCcorr.Aonly);
    %B tuned by stat sig TS only
    %diag_TC.B{ee} = diag(correlation_data{ee}.correlation.TCcorr.Bonly);
    %A selective (with added filters)
    diag_TC.Asel{ee} = diag(correlation_data{ee}.correlation.TCcorr.Aselective);
    %B selective (with added filters)
    diag_TC.Bsel{ee} = diag(correlation_data{ee}.correlation.TCcorr.Bselective);
    %A&B tuned - filtered; tuned to both criteria by either SI or TS
    %criterion
    diag_TC.AB{ee} = diag(correlation_data{ee}.correlation.TCcorr.si_ts.AB);
    %all neurons
    diag_TC.all{ee} = diag(correlation_data{ee}.correlation.TCcorr.all);
end

%combine correlation values into single matrix
%comb_TC.A = cell2mat(diag_TC.A');
%comb_TC.B = cell2mat(diag_TC.B');

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
%cell array - row - each subgroup; columns - each category
merge_new = cell(1,3);
merge_new{1} =  mean_TC.Asel;
merge_new{2} = mean_TC.Bsel;
merge_new{3} = mean_TC.AB;

%% Export data for TC correlation boxplots

tc_corr_sel_data.mean_TC.Asel = mean_TC.Asel;
tc_corr_sel_data.mean_TC.Bsel = mean_TC.Bsel;
tc_corr_sel_data.mean_TC.AB = mean_TC.AB;
%add individual neurons
tc_corr_sel_data.all_neurons.Asel = diag_TC.Asel; 
tc_corr_sel_data.all_neurons.Bsel = diag_TC.Bsel;
tc_corr_sel_data.all_neurons.AB = diag_TC.AB; 

if 0
%% Separate plot that includes the global remapping correlation scores (for grant)

%load the global r scores
load(fullfile('G:\Figure_2_3_selective_remap\cumulative_data_output','r_global.mat'),'r_global');

mean_TC.global = cellfun(@mean, r_global);

%plot the mean TC boxplots for A-sel,B-sel and AB
%cell array - row - each subgroup; columns - each category
merge_new = cell(1,4);
merge_new{1} = mean_TC.Asel;
merge_new{2} = mean_TC.Bsel;
merge_new{3} = mean_TC.AB;
merge_new{4} = mean_TC.global;


xlab={'1','2','3','4'} 
col=[220,20,60, 255; 65,105,225, 255; 139, 0, 139,255];  
col=col/255;

f=figure('Position',[2263 287 377 460])
hold on
title('Change colors manually in illustrator')
%plot zone separator lines
[f,x,group,positions,labelpos] =  multiple_boxplot(merge_new',xlab,{'Global','A&B','B','A'},col([1],:)') 
%overlay boxplot to add median line 
z= boxplot(x,group, 'positions', positions);
lines = findobj(z, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
set(lines, 'LineWidth',2)
xticks([1.25, 1.75, 2.25, 2.75])
xlim([1 3])
yticks([-0.3 -0.1 0 0.1 0.3 0.5 0.7])
ylim([-0.35 0.75])
xticklabels({'A sel.', 'B sel.', 'A&B','Global'})
ylabel('Correlation coef.')
%plot([1.7500 ,1.7500],[-10 110],'k--','LineWidth',1.5)
%plot([2.500 ,2.500],[-10 110],'k--','LineWidth',1.5)
set(gca,'FontSize',16)

%plot 0 line
plot([1 3],[0 0],'k-')

%% Statistics - TC correlation - paired Wilcoxon rank sum (with Dunn sidak correction for multi comp in Prism)

%place each set of mean in columns of a matrix
mean_matrix = [mean_TC.Asel', mean_TC.Bsel', mean_TC.AB'];


%si paired Wilcoxon tests (each type against other type)
for ii=1:3
    for jj=1:3
        [p_sr.si(ii,jj),~,stats_sr.si(ii,jj)] = signrank(mean_matrix(:,ii),mean_matrix(:,jj));
    end
end

%% Plot combined CDF of TC corr values

figure;
hold on;
c1 = cdfplot(comb_TC.Asel);
c1.Color = 'b';
c1.LineWidth = 2;
c2 = cdfplot(comb_TC.Bsel);
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

%% Plot this as ecdf curve fraction vs correlation (individual)
figure('Position', [2060 550 1400 420]);
subplot(1,3,1)
hold on
for ee=1:size(path_dir,2)
    %A selective
    c1 = cdfplot(diag(correlation_data{ee}.correlation.TCcorr.Aselective));
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
    c2 = cdfplot(diag(correlation_data{ee}.correlation.TCcorr.Bselective));
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
    c3 = cdfplot(diag(correlation_data{ee}.correlation.TCcorr.si_ts.AB));
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

end

