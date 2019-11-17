%% Define experiment core directories
%screened and good datasets for analysis
%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_120618_1'}; %
%field rate error (after PSEM silencing
path_dir{1} = {'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1'};
path_dir{2} = {'G:\Figure_2_3_selective_remap\I42R_AB_d1_032118_1'};

path_dir{3} = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1'};
path_dir{4} = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_2'};

path_dir{5} = {'G:\Figure_2_3_selective_remap\I53LT_AB_sal_113018_1'}; 
path_dir{6} = {'G:\Figure_2_3_selective_remap\I56_RTLS_AB_prePost_sal_042419_1'};

path_dir{7} = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_113018_1'};
path_dir{8} = {'G:\Figure_2_3_selective_remap\I57_RTLS_AB_prePost_792_042519_1'};

path_dir{9} = {'G:\Figure_2_3_selective_remap\I45_RT_AB_d1_062018_1'};
path_dir{10} = {'G:\Figure_2_3_selective_remap\I46_AB_d1_062018_1'};
path_dir{11} = {'G:\Figure_2_3_selective_remap\I57_LT_ABrand_no_punish_042119_1'};

%path_dir = {'G:\Figure_2_3_selective_remap\I47_LP_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I42R_AB_d1_032118_1'};

%path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I42L_AB_d1_032118_2'};
 
%path_dir = {'G:\Figure_2_3_selective_remap\I53LT_AB_sal_113018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I56_RTLS_AB_prePost_sal_042419_1'};

%path_dir = {'G:\Figure_2_3_selective_remap\I52RT_AB_sal_113018_1'};  
%path_dir = {'G:\Figure_2_3_selective_remap\I57_RTLS_AB_prePost_792_042519_1'};

%path_dir = {'G:\Figure_2_3_selective_remap\I45_RT_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I46_AB_d1_062018_1'};
%path_dir = {'G:\Figure_2_3_selective_remap\I57_LT_ABrand_no_punish_042119_1'};

%N =10; FOV = 11 (2 FOV from single animal)

%% Fraction of place cell (A sel, B sel, A&B remapping,neither) - bar chart and pie chart
%Figure 2C
fraction_place_cells(path_dir)

%% Cumulative task-selective STCs

cumulative_task_sel_STC(path_dir)

%% A/B/A&B/N cell SI/TS score distributions plots (Figure 2C/D)

si_ts_score_distributions(path_dir)

%% Centroid difference for A&B tuned neurons and centroid diff as fxn of max bin 
%scatterplot of centroid difference as a function of center between
%centroid of max place field - not used
%centroid_difference(path_dir)

%% Centroid distribution for A or B selective neurons - 
%how many bins to display
%1 - 20 (every 5) relative to 100 bins
%2 - 25 (every 4) relative to 100 bins
%3 - 50 (every 2) relative to 100 bins
options.bin_choose = 2;

centroid_dist_selective(path_dir,options)

%% Fractional distribution of remapping neuron subtypes - Figure 3C

frac_remapping_neurons(path_dir)

%% AUC/min scatterplots of A vs B neurons for each animal 

auc_scatterplots(path_dir)

%% PV/TC correlation

tc_pv_correlation_task_sel(path_dir)

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

%% Remapping centroids
%Figure 3E/F ( all Figure 3 code is here - organize this)
options.lowPVcorr = [6 7 8];
%organize this 
remapping_centroids(path_dir,options)

