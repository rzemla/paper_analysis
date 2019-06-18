%% Define experiment core directories
%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_120618_1'}; % field rate error
path_dir{1} = {'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1'};
path_dir{2} = {'G:\Figure_2_3_selective_remap\I42R_AB_d1_032118_1'};
path_dir{3} = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1'};
path_dir{4} = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_2'};
%path_dir = {'G:\Figure_2_3_selective_remap\I53LT_AB_sal_113018_1'}; %place field finder problem - adjust
%path_dir = {'G:\Figure_2_3_selective_remap\I56_RTLS_AB_prePost_sal_042419_1'}; %place field finder problem - adjust
path_dir{5} = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_113018_1'};
path_dir{6} = {'G:\Figure_2_3_selective_remap\I57_RTLS_AB_prePost_792_042519_1'};
path_dir{7} = {'G:\Figure_2_3_selective_remap\I45_RT_AB_d1_062018_1'};
path_dir{8} = {'G:\Figure_2_3_selective_remap\I46_AB_d1_062018_1'};
path_dir{9} = {'G:\Figure_2_3_selective_remap\I57_LT_ABrand_no_punish_042119_1'};


%% Fraction of A-selective and B-selective neurons by SI and TS criteria

%read relevant data
for ee=1:size(path_dir,2)
    load_data_path{ee} = fullfile(path_dir{ee},'cumul_analysis','frac_tuned.mat');
    fractional_data{ee} = load(string(load_data_path{ee}));
end

%create 1 matrix with all data - row = animal; column = respecitive counts
%A B A&B, neither
for ee=1:size(path_dir,2)
    %SI
    frac_mat_si(ee,:) = fractional_data{ee}.tuned_fractions.tuned_count{1};
    %TS
    frac_mat_ts(ee,:) = fractional_data{ee}.tuned_fractions.tuned_count{2};
end

%total counts by SI/TS
tuned_counts_si = sum(frac_mat_si,1);
tuned_counts_ts = sum(frac_mat_ts,1);
%total neurons by animal
neuron_counts_si = sum(frac_mat_si,2);
neuron_counts_ts = sum(frac_mat_ts,2);
%total neurons (match)
total_neurons = sum(neuron_counts_si);
total_neurons_ts = sum(neuron_counts_ts);

%fractions of total for each animal
frac_tuned_each.si = frac_mat_si./sum(frac_mat_si,2);
%mean
frac_tuned_each_mean.si = mean(frac_tuned_each.si,1);
%std
frac_tuned_each_std.si = std(frac_tuned_each.si,0,1);
%sem
frac_tuned_each_sem.si = frac_tuned_each_std.si./sqrt(size(frac_tuned_each.si,1));

%fractions of total for each animal
frac_tuned_each.ts = frac_mat_ts./sum(frac_mat_ts,2);
%mean
frac_tuned_each_mean.ts = mean(frac_tuned_each.ts,1);
%std
frac_tuned_each_std.ts = std(frac_tuned_each.ts,0,1);
%sem
frac_tuned_each_sem.ts = frac_tuned_each_std.ts./sqrt(size(frac_tuned_each.ts,1));

%fraction of all neurons tuned in each subgroup (cumulative)
frac_all_si = tuned_counts_si/total_neurons;
frac_all_ts = tuned_counts_ts/total_neurons;

%% Plot bar chart of fractions

%plot bar
figure('Position',[2010 380 870 420]);
subplot(1,2,1)
hold on;
title('Fraction tuned - S.I.');
%bar the mean for each group
b = bar(1:4,frac_tuned_each_mean.si,'FaceColor', 'flat');
pause(0.1)
%plot the sem for each mean for each group
for ib = 1:numel(b)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = b(ib).XData + b(ib).XOffset;
    
    errorbar(xData,frac_tuned_each_mean.si',frac_tuned_each_sem.si,'k.')
end

%set each bar to group color
b(1).CData(1:4,:) =  [0 0 1; 1 0 0; 1 0 1; 0.5 0.5 0.5];
xticks([1 2 3 4]);
xticklabels({'A','B','A&B', 'Neither'});
ylabel('Fraction of neurons');
ylim([0 0.5])
yticks([0:0.1:0.5])

subplot(1,2,2)
hold on;
title('Fraction tuned - T.S.');
%bar the mean for each group
b = bar(1:4,frac_tuned_each_mean.ts,'FaceColor', 'flat');
pause(0.1)
%plot the sem for each mean for each group
for ib = 1:numel(b)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = b(ib).XData + b(ib).XOffset;
    
    errorbar(xData,frac_tuned_each_mean.ts',frac_tuned_each_sem.ts,'k.')
end

%set each bar to group color
b(1).CData(1:4,:) =  [0 0 1; 1 0 0; 1 0 1; 0.5 0.5 0.5];
xticks([1 2 3 4]);
xticklabels({'A','B','A&B', 'Neither'});
ylabel('Fraction of neurons');
ylim([0 0.5])
yticks([0:0.1:0.5])




%% Plot pie chart for each type of tuning criterion

figure('Position',[2050 520 1140 450])
subplot(1,2,1)
p = pie(frac_all_si,{['A ', num2str(round(100*frac_all_si(1))), '%'],...
                        ['B ', num2str(round(100*frac_all_si(2))), '%'],...
                        ['A&B ', num2str(round(100*frac_all_si(3))), '%'],...
                        ['     Neither ', num2str(round(100*frac_all_si(4))), '%']});
hold on
title('Percentage of active neurons tuned (S.I.)');
colormap([0 0 1; 1 0 0; 1 0 1; 0.5 0.5 0.5])

subplot(1,2,2)
p = pie(frac_all_ts,{['A ', num2str(round(100*frac_all_ts(1))), '%'],...
                        ['B ', num2str(round(100*frac_all_ts(2))), '%'],...
                        ['A&B ', num2str(round(100*frac_all_ts(3))), '%'],...
                        ['     Neither ', num2str(round(100*frac_all_ts(4))), '%']});
hold on
title('Percentage of active neurons tuned (T.S.)');
colormap([0 0 1; 1 0 0; 1 0 1; 0.5 0.5 0.5])



%% AUC/min scatterplots of A vs B neurons for each animal 




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

%plot this an empirical distribution curve
figure;
subplot(2,1,1)
hold on
ylim([0 0.03])
histogram('BinEdges',0.5:1:100.5,'BinCounts',sum(centroid_mat_A,1),'Normalization', 'probability')
subplot(2,1,2)
hold on
ylim([0 0.03])
histogram('BinEdges',0.5:1:100.5,'BinCounts',sum(centroid_mat_B,1),'Normalization', 'probability')

%plot in polar space
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


%empirical cdf
figure;
hold on
ecdf((1:100),'Frequency',sum(centroid_mat_A,1)/sum(sum(centroid_mat_A,1)))
ecdf((1:100),'Frequency',sum(centroid_mat_B,1)/sum(sum(centroid_mat_B,1)))

% figure;
% hold on
% subplot(2,1,1)
% plot(sum(centroid_mat_A,1),'b')
% subplot(2,1,2)
% plot(sum(centroid_mat_B,1),'r')


%% PV/TC correlation

%% Place field analysis (width and number for selective neurons

%read relevant data
for ee=1:size(path_dir,2)
    load_data_path_pf{ee} = fullfile(path_dir{ee},'cumul_analysis','placeField_dist.mat');
    placeField_data{ee} = load(string(load_data_path_pf{ee}));
end

%combine field counts for Asel and Bsel into 1 matrix
%create matrices with centroid counts for each subclass of tuned neuron
for ee=1:size(path_dir,2)
    pf_count_mat_A(ee,:) = placeField_data{ee}.placeField_dist.field_count_A;
    pf_count_mat_B(ee,:) = placeField_data{ee}.placeField_dist.field_count_B;
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

