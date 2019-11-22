function [outputArg1,outputArg2] = event_speed_scatter(path_dir)


%% Load in lap speed data

%read relevant data
for ee=1:size(path_dir,2)
    load_data_path{ee} = fullfile(path_dir{ee},'cumul_analysis','lap_and_event_speed.mat');
    event_speed_data{ee} = load(string(load_data_path{ee}),'mean_event_speed');
end


%save(fullfile(path_dir{1},'cumul_analysis','lap_and_event_speed.mat'),'mean_bin_speed', 'lap_bin_split','mean_event_speed');

%% Cumulate all points for scatterplot

%extract speed data for selective neurons from each animal
for ii=1:size(path_dir,2)
    %A speed, then B speed order
    %A selective
    cum_speed.Asel{ii,1} = event_speed_data{ii}.mean_event_speed.Asel';
    cum_speed.Asel{ii,2} = event_speed_data{ii}.mean_event_speed.Asel_B_laps';
    %B selective
    cum_speed.Bsel{ii,1} = event_speed_data{ii}.mean_event_speed.Bsel_A_laps';
    cum_speed.Bsel{ii,2} = event_speed_data{ii}.mean_event_speed.Bsel';
end

%combine into single matrix
Asel_speed_cum =cell2mat(cum_speed.Asel);
Bsel_speed_cum =cell2mat(cum_speed.Bsel);

%% Mean difference for each neuron (cumulative)

%cumulative speed difference between A and B
Asel_speed_diff = Asel_speed_cum(:,1) - Asel_speed_cum(:,2);
Bsel_speed_diff = Bsel_speed_cum(:,1) - Bsel_speed_cum(:,2);

nb_Asel = sum(~isnan(Asel_speed_diff));
nb_Bsel = sum(~isnan(Bsel_speed_diff));

sem_Asel = nanstd(Asel_speed_diff)./sqrt(nb_Asel);
sem_Bsel = nanstd(Bsel_speed_diff)./sqrt(nb_Bsel);

mean_Asel = nanmean(Asel_speed_diff);
mean_Bsel = nanmean(Bsel_speed_diff);

%95% probability intervals of t-distribution
ci95_Asel = tinv([0.025 0.975], nb_Asel-1);  
%33% probability intervals of t-disttribution
ci33_Asel = tinv([0.165 0.835], nb_Asel-1); 

%95% probability intervals of t-distribution
ci95_Bsel = tinv([0.025 0.975], nb_Bsel-1);  
%33% probability intervals of t-disttribution
ci33_Bsel = tinv([0.165 0.835], nb_Bsel-1);  

%calculate 95% confidence interval
Asel_95ci = ci95_Asel.*sem_Asel;
Asel_33ci = ci33_Asel.*sem_Asel;

Bsel_95ci = ci95_Bsel.*sem_Bsel;

%% Find fraction of A sel neurons and B cell neurons that are within +5/-5 speed range diff

ROI_wn_5cms.A = find(Asel_speed_diff <=5 & Asel_speed_diff >=-5);
ROI_wn_5cms.B = find(Bsel_speed_diff <=5 & Bsel_speed_diff >=-5);

%fraction of A cells within -5/+5 cm/s
fracA = size(ROI_wn_5cms.A,1)/nb_Asel;
fracB = size(ROI_wn_5cms.B,1)/nb_Bsel;

Asel_1SD_range = [mean_Asel - nanstd(Asel_speed_diff), mean_Asel + nanstd(Asel_speed_diff)];
Bsel_1SD_range = [mean_Bsel - nanstd(Bsel_speed_diff), mean_Bsel + nanstd(Bsel_speed_diff)];

skewA = skewness(Asel_speed_diff);
skewB = skewness(Bsel_speed_diff);

kurtA = kurtosis(Asel_speed_diff);
kurtB = kurtosis(Bsel_speed_diff);

%% Statistics - test distribution for normslity
%Shapiro-Wilk test (p> 0.05 - data is normally distrubuted)
[h_sw, p_sw, sw_stat] = swtest(Asel_speed_diff)

[h_sw, p_sw, sw_stat] = swtest(Bsel_speed_diff)

%no normal
h = lillietest(Asel_speed_diff)

%% Colors
color_mat=[ 65,105,225; 220,20,60; 139, 0, 139]./255;

%plot histogram
figure('Position',[2384 433 920 420])
subplot(1,2,1)
hold on
title('A selective')
axis square
yticks([0 0.1 0.2])
ylabel('Fraction of neurons')
xlim([-20 20])
hA = histogram(Asel_speed_diff,[-20:2:20],'Normalization','probability');
hA.FaceColor = color_mat(1,:);
xlabel('Speed difference (A-B) [cm/s]')
%5 cm/s cutoffs
plot([-5 -5],[0 0.25],'k','LineWidth',1)
plot([5 5],[0 0.25],'k','LineWidth',1)

set(gca,'FontSize',12)
set(gca,'LineWidth',1.5)

subplot(1,2,2)
hold on
title('B selective')
axis square
yticks([0 0.1 0.2])
ylabel('Fraction of neurons')
xlim([-20 20])
hA = histogram(Bsel_speed_diff,[-20:2:20],'Normalization','probability');
hA.FaceColor = color_mat(2,:);
xlabel('Speed difference (A-B) [cm/s]')
%5 cm/s cutoffs
plot([-5 -5],[0 0.25],'k','LineWidth',1)
plot([5 5],[0 0.25],'k','LineWidth',1)

set(gca,'FontSize',12)
set(gca,'LineWidth',1.5)


%% Plot scatters 
figure('Position',[2247 400 1060 420])
subplot(1,2,1)
hold on
title('A selective')
axis square
xlim([0 30])
ylim([0 30])
xlabel('A lap speed [cm/s]')
ylabel('B lap speed [cm/s]')
scatter(Asel_speed_cum(:,1),Asel_speed_cum(:,2),6,'filled','MarkerFaceColor',color_mat(1,:))
plot([0 30],[0 30],'k--')

set(gca,'FontSize',12)
set(gca,'LineWidth',1.5)

subplot(1,2,2)
hold on
title('B selective')
axis square
xlim([0 30])
ylim([0 30])
xlabel('A lap speed [cm/s]')
ylabel('B lap speed [cm/s]')
scatter(Bsel_speed_cum(:,1),Bsel_speed_cum(:,2),6,'filled','MarkerFaceColor',color_mat(2,:))
plot([0 30],[0 30],'k--')

set(gca,'FontSize',12)
set(gca,'LineWidth',1.5)



end

