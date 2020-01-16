function [outputArg1,outputArg2] = centroid_dist_selective(path_dir,options)

%read relevant data
for ee=1:size(path_dir,2)
    load_data_path_cent{ee} = fullfile(path_dir{ee},'cumul_analysis','centroid.mat');
    centroid_data{ee} = load(string(load_data_path_cent{ee}));
end

%create matrices with centroid counts for each subclass of tuned neuron
for aa=1:size(path_dir,2)
    %bin count accross 100 bins
    centroid_mat_A(aa,:) = centroid_data{aa}.centroid_ct.A;
    centroid_mat_B(aa,:) = centroid_data{aa}.centroid_ct.B;
    
    %bin assignment for each neuron
    bin_assign.A{aa} = centroid_data{aa}.centroid_bins.A; 
    bin_assign.B{aa} = centroid_data{aa}.centroid_bins.B;
end

%normalize to itself and get mean for A and B (100 bins)
centroid_mat_A_norm = centroid_mat_A./sum(centroid_mat_A,2);
centroid_mat_B_norm = centroid_mat_B./sum(centroid_mat_B,2);

%take mean in each bin from each animal
centroid_mat_A_norm_mean = mean(centroid_mat_A_norm,1);
centroid_mat_B_norm_mean = mean(centroid_mat_B_norm,1);

%get normalized standard deviation for centroid of task-selective neurons
centroid_mat_A_norm_std = std(centroid_mat_A_norm,0,1);
centroid_mat_B_norm_std = std(centroid_mat_B_norm,0,1);

%% Downsample to respective bin nb (50 25 20) from 100 bin input
%1 - 20 (every 5) relative to 100 bins
%2 - 25 (every 4) relative to 100 bins
%3 - 50 (every 2) relative to 100 bins
bin_choose = options.bin_choose;

bin_skip = [5, 4 ,2];
bin_nb = [20, 25, 50];

%generate count for fewer bin sizes
for aa=1:size(path_dir,2)
    for bb=1:3
        %A
        centroid_count_ds.A{bb}(aa,:) = histcounts(bin_assign.A{aa}, 0.5:bin_skip(bb):100.5);
        %B
        centroid_count_ds.B{bb}(aa,:) = histcounts(bin_assign.B{aa}, 0.5:bin_skip(bb):100.5);
    end
end

%normalized counts for each bin range
for bb=1:3
    %A
    centroid_count_ds_norm.A{bb} = centroid_count_ds.A{bb}./sum(centroid_count_ds.A{bb},2);
    %B
    centroid_count_ds_norm.B{bb} = centroid_count_ds.B{bb}./sum(centroid_count_ds.B{bb},2);
end

%% Generate plot of mean and sem at each bin (25), plot and do stats

%
figure
hold on
ylim([0 0.2])
plot(mean(centroid_count_ds_norm.A{2},1),'b')
plot(mean(centroid_count_ds_norm.B{2},1),'r')

for bb=1:25
    p(bb) = signrank(centroid_count_ds_norm.A{2}(:,bb), centroid_count_ds_norm.B{2}(:,bb));
end

%% Resume here to generate linearize plot and get Rayleigh test
figure
hold on
ylim([0 0.2])
yticks([0 0.05 0.1 0.15 0.2])
h1 = histogram(common_bins_combined(1,:),0:10:100,'Normalization','probability');
xticks([0 100])
xticklabels({'0','1'})
xlabel('Normalized position')
ylabel('Normalized density')
h1.FaceColor = [139 0 139]./255;
h1.FaceAlpha = 1;

plot([B_zone_end B_zone_end], [0 0.2],'LineStyle','--', 'LineWidth',2,'Color',[220,20,60]./255)
plot([A_zone_end A_zone_end], [0 0.2],'LineStyle','--', 'LineWidth',2,'Color',[65,105,225]./255)

%% Do Rayleigh test of uniformity on dataset
%convert bins to radians
common_rad_angles = circ_ang2rad(common_bins_combined(1,:)./100.*360);
[pval, z] = circ_rtest(common_rad_angles)

%% Generate shuffle data
%get number of ROI for each animal for A and B laps
nb_ROI_A = sum(centroid_mat_A,2);
nb_ROI_B = sum(centroid_mat_B,2);

nb_bins = bin_nb(bin_choose);
nb_shuffle = 1000;

%randomize location of max centroid
for aa=1:size(path_dir,2)
    shuffled_bin.A{aa} = randi([1 nb_bins],nb_ROI_A(aa),nb_shuffle);
    shuffled_bin.B{aa} = randi([1 nb_bins],nb_ROI_B(aa),nb_shuffle);
end

%count number of events in each bin
for aa=1:size(path_dir,2)
    %generate a histcounts
    for ss=1:nb_shuffle
        [bin_count_shuffle.A{aa}(:,ss),~] = histcounts(shuffled_bin.A{aa}(:,ss), 0.5:1:(nb_bins+0.5));
        [bin_count_shuffle.B{aa}(:,ss),~] = histcounts(shuffled_bin.B{aa}(:,ss), 0.5:1:(nb_bins+0.5));
    end
end

%normalize each bin count
for aa=1:size(path_dir,2)
    bin_count_shuffle_norm.A{aa} = bin_count_shuffle.A{aa}./sum(bin_count_shuffle.A{aa});
    bin_count_shuffle_norm.B{aa} = bin_count_shuffle.B{aa}./sum(bin_count_shuffle.B{aa});
end

%take mean across each bin for each animal
for aa=1:size(path_dir,2)
    bin_shuffle_mean.A(:,aa) = mean(bin_count_shuffle_norm.A{aa},2);
    bin_shuffle_mean.B(:,aa) = mean(bin_count_shuffle_norm.B{aa},2);
end


%get 3 std thresholds
thres_A = mean(bin_shuffle_mean.A(:)) + 3*std(bin_shuffle_mean.A(:));
thres_B = mean(bin_shuffle_mean.B(:)) + 3*std(bin_shuffle_mean.B(:));

%% Range of bins above threshold for A or B
%for 25 bins
%at least 2 neighboring bins above threshold to be considered in sig range

%define the cm range for the sig bins

%bin width of eahc bin
binSize_cm = 196/25;
%
bin_distance = [0:binSize_cm:196];
%cm centers of the each bin
bin_center_cm = bin_distance + (binSize_cm/2);
bin_center_cm = bin_center_cm(1:end-1);

%find bin range with at least 2 bins neighboring
sig_bins.A = mean(centroid_count_ds_norm.A{bin_choose},1)>thres_A
sig_bins.B = mean(centroid_count_ds_norm.B{bin_choose},1)>thres_B

bin_center_cm(sig_bins.A)
bin_center_cm(sig_bins.B)

%% Define colors
%blue, red, light blue, light red      
color_mat = [65,105,225; 220,20,60; 135,206,250; 240,128,128]./255;

%% Run ks test between the all neurons and generate cdf plot 

%all neurons
[~,p_all_neurons,ks2stat_all_neurons] = kstest2(cell2mat(bin_assign.A),cell2mat(bin_assign.B));

%all means


%every 4 cm from 18-100
edges = [1:4:100];
edge_center = edges(1:end-1)+1;
for aa=1:size(path_dir,2)
    [N_Asel(aa,:),~] = histcounts(bin_assign.A{aa} ,edges,'Normalization','probability');
    [N_Bsel(aa,:),~] = histcounts(bin_assign.B{aa},edges,'Normalization','probability');
end

%get cumulative probability
cum_prob_A_sel = cumsum(N_Asel,2);
cum_prob_B_sel = cumsum(N_Bsel,2);

%get sem at each fractional point
sem_cum_prob_A = nanstd(cum_prob_A_sel,0,1)./sqrt(size(path_dir,2));
sem_cum_prob_B = nanstd(cum_prob_B_sel,0,1)./sqrt(size(path_dir,2));

lw = 2;
figure
hold on
axis square
xlabel('Normalized position')
ylabel('Fraction of place cells')
yticks([0 0.2 0.4 0.6 0.8 1])

%e3 = errorbar(edge_center,mean(cum_prob_AB_A,1),sem_cum_prob_AB_A,'LineWidth',lw,'Color',color_mat(3,:));
%e4 = errorbar(edge_center,mean(cum_prob_AB_B,1),sem_cum_prob_AB_B,'LineWidth',lw,'Color',color_mat(4,:));
e1 = errorbar(edge_center,mean(cum_prob_A_sel,1),sem_cum_prob_A,'LineWidth',lw,'Color',color_mat(1,:));
e2 = errorbar(edge_center,mean(cum_prob_B_sel,1),sem_cum_prob_B,'LineWidth',lw,'Color',color_mat(2,:));

legend([e1 e2],{'A','B'},'Location','southeast')
set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)


% figure;
% hold on;
% axis square
% c1 = cdfplot(cell2mat(bin_assign.A));
% c1.Color = [65,105,225]./255;
% c1.LineWidth = 2;
% c2 = cdfplot(cell2mat(bin_assign.B));
% c2.Color = [220,20,60]./255;
% c2.LineWidth = 2;
% 
% xlabel('Normalized position');
% ylabel('Cumulative fraction');
% grid off
% set(gca,'FontSize',14)
% set(gca,'LineWidth',1.5)
% %xlim([-0.5 1]);
% xticks([1 50 100])
% xticklabels({'0', '0.5', '100'})
% legend([c1 c2],{'A','B'},'Location','northwest')

%% Get distribution histograms for A and B trials and 3*std line
figure
subplot(1,2,1)
hold on
ylim([0 0.5])
histogram(bin_shuffle_mean.A,'Normalization','probability')
%plot upper threshold
plot([thres_A thres_A], [0 0.2],'r')
xlabel('Bin density')
ylabel('Probability')

subplot(1,2,2)
hold on
ylim([0 0.5])
histogram(bin_shuffle_mean.B,'Normalization','probability')
%plot upper threshold
plot([thres_B thres_B], [0 0.2],'r')
xlabel('Bin density')
ylabel('Probability')

%define marker vectors
reward_A_vector = exp(1i*2*pi*0.70);
reward_A_vector_end = exp(1i*2*pi*0.75);

reward_B_vector = exp(1i*2*pi*0.30);
reward_B_vector_end = exp(1i*2*pi*0.35);

lap_start_vector = exp(1i*2*pi*0);
odor_zone_end_vector = exp(1i*2*pi*0.1);

%% Polar plot mean
figure('Position',[2220 146 1327 707]);
pax1 = subplot(1,2,1,polaraxes);
hold on

%plot the relevant lap markers
%lap start
polarplot([0+0i,0.15*lap_start_vector],'Color',[1 1 1]*0.7,'LineWidth',1.5)
%odor zone end
polarplot([0+0i,0.15*odor_zone_end_vector],'Color',[0 146 69]./255,'LineWidth',1.5)
%reward zone A start
polarplot([0+0i,0.15*reward_A_vector],'Color',[65,105,225]./255,'LineWidth',1.5)
%reward zone A end
polarplot([0+0i,0.15*reward_A_vector_end],'Color',[65,105,225]./255,'LineWidth',1.5)
%reward zone B start
polarplot([0+0i,0.15*reward_B_vector],'Color',[220,20,60]./255,'LineWidth',1.5)
%reward zone B start
polarplot([0+0i,0.15*reward_B_vector_end],'Color',[220,20,60]./255,'LineWidth',1.5)

title('A selective')
pax1.ThetaAxisUnits = 'degrees';
pax1.ThetaTick = [0,rad2deg(angle(odor_zone_end_vector)),rad2deg(angle(reward_B_vector)),rad2deg(angle(reward_A_vector))+360,];
pax1.ThetaTickLabel = {'Lap start','Odor end','B reward','A reward'};
pax1.RAxisLocation = 330;
pax1.RLim = [0 0.15];
pax1.RColor = 'k';
pax1.ThetaColorMode = 'manual';
pax1.ThetaColor = 'k';
pax1.GridColorMode = 'manual';
pax1.GridColor = 'k';
pax1.GridAlpha = 0.5;
pax1.FontSize = 18;
polarhistogram(pax1,'BinEdges',0:(2*pi)/bin_nb(bin_choose):2*pi,'BinCounts', mean(centroid_count_ds_norm.A{bin_choose},1),...
    'FaceColor',[65,105,225]./255,'FaceAlpha',.9)

%circle random dist
th = linspace(0,2*pi,50);
polarplot(th,thres_A+zeros(size(th)),'--','LineWidth',1.5,'Color',[255,69,0]./255)


pax2 = subplot(1,2,2,polaraxes);
hold on

%plot the relevant lap markers
%lap start
polarplot([0+0i,0.15*lap_start_vector],'Color',[1 1 1]*0.7,'LineWidth',1.5)
%odor zone end
polarplot([0+0i,0.15*odor_zone_end_vector],'Color',[0 146 69]./255,'LineWidth',1.5)
%reward zone A start
polarplot([0+0i,0.15*reward_A_vector],'Color',[65,105,225]./255,'LineWidth',1.5)
%reward zone A end
polarplot([0+0i,0.15*reward_A_vector_end],'Color',[65,105,225]./255,'LineWidth',1.5)
%reward zone B start
polarplot([0+0i,0.15*reward_B_vector],'Color',[220,20,60]./255,'LineWidth',1.5)
%reward zone B start
polarplot([0+0i,0.15*reward_B_vector_end],'Color',[220,20,60]./255,'LineWidth',1.5)

title('B selective')
pax2.ThetaAxisUnits = 'degrees';
pax2.ThetaTick = [0,rad2deg(angle(odor_zone_end_vector)),rad2deg(angle(reward_B_vector)),rad2deg(angle(reward_A_vector))+360,];
pax2.ThetaTickLabel = {'Lap start','Odor end','B reward','A reward'};
pax2.RAxisLocation = 330;
pax2.RLim = [0 0.15];
pax2.RColor = 'k';
pax2.ThetaColorMode = 'manual';
pax2.ThetaColor = 'k';
pax2.GridColorMode = 'manual';
pax2.GridColor = 'k';
pax2.GridAlpha = 0.5;
pax2.FontSize = 18;
polarhistogram(pax2,'BinEdges',0:(2*pi)/bin_nb(bin_choose):2*pi,'BinCounts', mean(centroid_count_ds_norm.B{bin_choose},1),...
    'FaceColor',[220,20,60]./255,'Facealpha',0.9)

%circle random dist
th = linspace(0,2*pi,50);
polarplot(th,thres_B+zeros(size(th)),'--','LineWidth',2,'Color',[255,69,0]./255)

%empirical cdf
if 0
    figure;
    hold on
    ecdf((1:100),'Frequency',sum(centroid_mat_A,1)/sum(sum(centroid_mat_A,1)))
    ecdf((1:100),'Frequency',sum(centroid_mat_B,1)/sum(sum(centroid_mat_B,1)))
end




%% Plot as fraction of neurons tuned in each bin with shaded std
if 0
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
    'FaceColor',[65,105,225]./255)

pax2 = subplot(1,2,2,polaraxes);
hold on
title('B selective')
pax2.ThetaAxisUnits = 'degrees';
pax2.RAxisLocation = 45;
pax2.RLim = [0 0.04];
pax2.RColor = 'r';
pax2.FontSize = 14;
polarhistogram(pax2,'BinEdges',0:(2*pi)/100:2*pi,'BinCounts', sum(centroid_mat_B,1),'Normalization', 'probability',...
    'FaceColor',[220,20,60]./255)


% figure;
% hold on
% subplot(2,1,1)
% plot(sum(centroid_mat_A,1),'b')
% subplot(2,1,2)
% plot(sum(centroid_mat_B,1),'r')
end


end

