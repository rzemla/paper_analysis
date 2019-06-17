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

%fraction of all neurons tuned in each subgroup
frac_all_si = tuned_counts_si/total_neurons;
frac_all_ts = tuned_counts_ts/total_neurons;

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



