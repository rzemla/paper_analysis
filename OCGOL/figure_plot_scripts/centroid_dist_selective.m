function [outputArg1,outputArg2] = centroid_dist_selective(path_dir)

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

%% Generate shuffle data
%get number of ROI for each animal for A and B laps
nb_ROI_A = sum(centroid_mat_A,2);
nb_ROI_B = sum(centroid_mat_B,2);

nb_shuffle = 1000;

%randomize location of max centroid
for aa=1:size(path_dir,2)
    shuffled_bin.A{aa} = randi([1 100],nb_ROI_A(aa),nb_shuffle);
    shuffled_bin.B{aa} = randi([1 100],nb_ROI_B(aa),nb_shuffle);
end

%count number of events in each bin
for aa=1:size(path_dir,2)
    %generate a histcounts
    for ss=1:nb_shuffle
        [bin_count_shuffle.A{aa}(:,ss),~] = histcounts(shuffled_bin.A{aa}(:,ss), 0.5:1:100.5);
        [bin_count_shuffle.B{aa}(:,ss),~] = histcounts(shuffled_bin.B{aa}(:,ss), 0.5:1:100.5);
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

%get distribution histograms for A and B trials and 3*std line
figure
subplot(1,2,1)
hold on
ylim([0 0.2])
histogram(bin_shuffle_mean.A,'Normalization','probability')
%plot upper threshold
plot([thres_A thres_A], [0 0.2],'r')
xlabel('Bin density')
ylabel('Probability')

subplot(1,2,2)
hold on
ylim([0 0.2])
histogram(bin_shuffle_mean.B,'Normalization','probability')
%plot upper threshold
plot([thres_B thres_B], [0 0.2],'r')
xlabel('Bin density')
ylabel('Probability')



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

th = linspace(0,2*pi,50);
r_thres_A = 0.01;
polarplot(th,r_thres_A+zeros(size(th)),'--','LineWidth',1.5,'Color','k')

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




%% Plot as fraction of neurons tuned in each bin with shaded std
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


% figure;
% hold on
% subplot(2,1,1)
% plot(sum(centroid_mat_A,1),'b')
% subplot(2,1,2)
% plot(sum(centroid_mat_B,1),'r')



end

