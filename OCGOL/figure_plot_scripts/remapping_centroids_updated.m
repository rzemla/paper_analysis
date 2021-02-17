function [reward_zones_all_animal,partial_idx_by_animal_zone,remap_prop_figs,global_dist_scatter] = remapping_centroids_updated(path_dir)

%% Load in data from each session directory
for ee=1:size(path_dir,2)
    load_data_path{ee} = fullfile(path_dir{ee},'cumul_analysis','place_field_centers_remap.mat');
    binCenter_data{ee} = load(string(load_data_path{ee}));
end


%% Load in reward positions from each animal
for ee=1:size(path_dir,2)
    load_reward_path{ee} = fullfile(path_dir{ee},'cumul_analysis','reward_positions.mat');
    reward_pos_data{ee} = load(string(load_reward_path{ee}));
end

%% Load in the indices of the neurons for each animal to identify partial remappers for visualization in supp.

for ee=1:size(path_dir,2)
    load_neuron_idx_path{ee} = fullfile(path_dir{ee},'cumul_analysis','remap_corr_idx.mat');
    remapping_idxs{ee} = load(string(load_neuron_idx_path{ee}));
end


%% Mean normalized reward 
for ee=1:size(path_dir,2)
    rew_A_pos(ee) = reward_pos_data{ee}.reward_positions.mean_rew_A_norm;
    rew_B_pos(ee) = reward_pos_data{ee}.reward_positions.mean_rew_B_norm;
    track_end_cm(ee) = reward_pos_data{ee}.reward_positions.mean_track_end_cm;
end

%mean track end position
mean_track_end_cm = mean(track_end_cm);

%mean reward zone start position
mean_rew_A_pos = mean(rew_A_pos);
mean_rew_B_pos = mean(rew_B_pos);

%10 cm zone
reward_10cm_dist = 10/mean_track_end_cm;

%zone ends in bins
A_zone_start = 100*round(mean_rew_A_pos,2);
B_zone_start = 100*round(mean_rew_B_pos,2);
%end zone
A_zone_end = 100*round(mean_rew_A_pos + reward_10cm_dist,2);
B_zone_end = 100*round(mean_rew_B_pos + reward_10cm_dist,2);



%% Bundle reward bins for export

reward_zones_all_animal.A_zone_start = A_zone_start;
reward_zones_all_animal.B_zone_start = B_zone_start;
reward_zones_all_animal.A_zone_end = A_zone_end;
reward_zones_all_animal.B_zone_end = B_zone_end;


%% Simple scatter plot for common/global remapper (global remapping field scatter for supplement figure)

figure
hold on
title('Global')
xlabel('Normalized location of place field A')
ylabel('Normalized location of place field B')
%generate axis for place fields
xticks(0:10:100)
yticks(0:10:100)
xlabels_norm = 0:0.1:1;
x_labels = {};

for xt=1:11
    x_labels{xt} = num2str(xlabels_norm(xt));
end

xticklabels(x_labels)
yticklabels(x_labels)

for ee=1:size(path_dir,2)
    scatter(binCenter_data{ee}.bin_center.final.global(1,:),binCenter_data{ee}.bin_center.final.global(2,:),12,'filled',...
            'MarkerEdgeColor',[139,0,139]./255,'MarkerFaceColor',[139,0,139]./255)
end

plot([0 100],[0 100],'k','LineWidth',2)

plot([B_zone_end B_zone_end], [0 100],'k--', 'LineWidth',2)
plot([A_zone_end A_zone_end], [0 100],'k--', 'LineWidth',2)

%% Export scatter data

global_dist_scatter.binCenter_data = binCenter_data;


% figure
% hold on
% title('Common')
% xlabel('Center location of Place field A')
% ylabel('Center location of Place field B')
% for ee=1:size(path_dir,2)
%     scatter(binCenter_data{ee}.bin_center.final.common(1,:),binCenter_data{ee}.bin_center.final.common(2,:),'filled','m')
% end
% 
% plot([0 100],[0 100],'k--')
% 
% plot([B_zone_end B_zone_end], [0 B_zone_end],'k--')
% plot([A_zone_end A_zone_end], [0 A_zone_end],'k--')


%% Look for difference and scatter
%combine global place fields into 1 matrix
for ee=1:size(path_dir,2)
    A_global_bincenter{ee} = binCenter_data{ee}.bin_center.final.global(1,:);
    B_global_bincenter{ee} = binCenter_data{ee}.bin_center.final.global(2,:);
    
    A_common_bincenter{ee} = binCenter_data{ee}.bin_center.final.common(1,:);
    B_common_bincenter{ee} = binCenter_data{ee}.bin_center.final.common(2,:);
    
end

%combine bins for common and global
global_bins_combined = [cell2mat(A_global_bincenter); cell2mat(B_global_bincenter)];
common_bins_combined = [cell2mat(A_common_bincenter); cell2mat(B_common_bincenter)]; 

%get difference
A_minus_B_global = global_bins_combined(1,:) - global_bins_combined(2,:);
A_minus_B_common = common_bins_combined(1,:) - common_bins_combined(2,:);

%scatter of A-B
% figure
% hold on
% ylabel('A-B')
% boxplot(A_minus_B_global)
%scatter(ones(size(A_minus_B_global)),A_minus_B_global)

[p,h,~] = signrank(A_minus_B_common)

[p,h,~] = signrank(A_minus_B_global)

% zone_global_idx(1).A = find(global_far.A > 0 & global_far.A <= B_zone_end);
% zone_global_idx(2).A = find(global_far.A > B_zone_end & global_far.A <= A_zone_end);
% zone_global_idx(3).A = find(global_far.A > A_zone_end & global_far.A <= 100);
%define the zones of the track
mean(A_minus_B_common)

%% Split global remappers between A & B into distinct zones (combined)
%QC checked
% for ee=1:size(path_dir,2)
%     scatter(binCenter_data{ee}.bin_center.final.common(1,:),binCenter_data{ee}.bin_center.final.common(2,:),'filled','m')
% end

%preallocate
ct.zoneI = 0;
ct.zoneII = 0;
ct.zoneIII = 0;

ct.zoneI_II = 0;
ct.zoneI_III = 0;

ct.zoneII_I = 0;
ct.zoneII_III = 0;

ct.zoneIII_II = 0;
ct.zoneIII_I = 0;

%for each ROI
for rr=1:size(global_bins_combined,2)
    %if both in zone I
    if (global_bins_combined(1,rr) > 0 && global_bins_combined(1,rr) <= B_zone_end) &&...
            (global_bins_combined(2,rr) > 0 && global_bins_combined(2,rr) <= B_zone_end)
        ct.zoneI = ct.zoneI + 1;
    end

    %if both in zone II
    if (global_bins_combined(1,rr) > B_zone_end && global_bins_combined(1,rr) <= A_zone_end) &&...
            (global_bins_combined(2,rr) > B_zone_end && global_bins_combined(2,rr) <= A_zone_end)
        ct.zoneII = ct.zoneII + 1;
    end    

    %if both in zone III
    if (global_bins_combined(1,rr) > A_zone_end && global_bins_combined(1,rr) <= 100) &&...
            (global_bins_combined(2,rr) > A_zone_end && global_bins_combined(2,rr) <= 100)
        ct.zoneIII = ct.zoneIII + 1;
    end       
    
    %if both in zone I --> II
    if (global_bins_combined(1,rr) > 0 && global_bins_combined(1,rr) <= B_zone_end) &&...
            (global_bins_combined(2,rr) > B_zone_end && global_bins_combined(2,rr) <= A_zone_end)
        ct.zoneI_II = ct.zoneI_II + 1;
    end 

    %if both in zone I --> III
    if (global_bins_combined(1,rr) > 0 && global_bins_combined(1,rr) <= B_zone_end) &&...
            (global_bins_combined(2,rr) > A_zone_end && global_bins_combined(2,rr) <= 100)
        ct.zoneI_III = ct.zoneI_III + 1;
    end     
    
    %if both in zone II --> I
    if (global_bins_combined(1,rr) > B_zone_end && global_bins_combined(1,rr) <= A_zone_end) &&...
            (global_bins_combined(2,rr) > 0 && global_bins_combined(2,rr) <= B_zone_end)
        ct.zoneII_I = ct.zoneII_I + 1;
    end       
    
    %if both in zone II --> III
    if (global_bins_combined(1,rr) > B_zone_end && global_bins_combined(1,rr) <= A_zone_end) &&...
            (global_bins_combined(2,rr) > A_zone_end && global_bins_combined(2,rr) <= 100)
        ct.zoneII_III = ct.zoneII_III + 1;
    end   
    
    %if both in zone III --> II
    if (global_bins_combined(1,rr) > A_zone_end && global_bins_combined(1,rr) <= 100) &&...
            (global_bins_combined(2,rr) > B_zone_end && global_bins_combined(2,rr) <= A_zone_end)
        ct.zoneIII_II = ct.zoneIII_II + 1;
    end   

    %if both in zone III --> I
    if (global_bins_combined(1,rr) > A_zone_end && global_bins_combined(1,rr) <= 100) &&...
            (global_bins_combined(2,rr) > 0 && global_bins_combined(2,rr) <= B_zone_end)
        ct.zoneIII_I = ct.zoneIII_I + 1;
    end     
    
end

all_zone_counts = [ct.zoneI, ct.zoneII, ct.zoneIII, ct.zoneI_II,ct.zoneI_III,...
    ct.zoneII_I, ct.zoneII_III, ct.zoneIII_II,ct.zoneIII_I];

%plot bar chart
% figure
% hold on
% bar(all_zone_counts)
% xticks(1:9)
% xticklabels({'I-I','II-II','III-III','I-II','I-III','II-I','II-III','III-II','III-I'})

% zone_global_idx(1).A = find(global_far.A > 0 & global_far.A <= B_zone_end);
% zone_global_idx(2).A = find(global_far.A > B_zone_end & global_far.A <= A_zone_end);
% zone_global_idx(3).A = find(global_far.A > A_zone_end & global_far.A <= 100);

%% Split global remappers between A & B into distinct zones (combined) (Fig. 3g)
%QC checked

session_nb = size(A_global_bincenter,2);

%preallocate
ct_split.zoneI = zeros(session_nb,1);
ct_split.zoneII = zeros(session_nb,1);
ct_split.zoneIII = zeros(session_nb,1);

ct_split.zoneI_II = zeros(session_nb,1);
ct_split.zoneI_III = zeros(session_nb,1);

ct_split.zoneII_I = zeros(session_nb,1);
ct_split.zoneII_III = zeros(session_nb,1);

ct_split.zoneIII_II = zeros(session_nb,1);
ct_split.zoneIII_I = zeros(session_nb,1);

%merge together with by-animal split of global remapping neurons
for aa=1:session_nb
    global_bincenter_split{aa} = [A_global_bincenter{aa}; B_global_bincenter{aa}];
end


%for each session
for ss=1:session_nb
    %for each ROI
    for rr=1:size(global_bincenter_split{ss},2)
        %if both in zone I
        if (global_bincenter_split{ss}(1,rr) > 0 && global_bincenter_split{ss}(1,rr) <= B_zone_end) &&...
                (global_bincenter_split{ss}(2,rr) > 0 && global_bincenter_split{ss}(2,rr) <= B_zone_end)
            ct_split.zoneI(ss) = ct_split.zoneI(ss) + 1;
        end
        
        %if both in zone II
        if (global_bincenter_split{ss}(1,rr) > B_zone_end && global_bincenter_split{ss}(1,rr) <= A_zone_end) &&...
                (global_bincenter_split{ss}(2,rr) > B_zone_end && global_bincenter_split{ss}(2,rr) <= A_zone_end)
            ct_split.zoneII(ss) = ct_split.zoneII(ss) + 1;
        end
        
        %if both in zone III
        if (global_bincenter_split{ss}(1,rr) > A_zone_end && global_bincenter_split{ss}(1,rr) <= 100) &&...
                (global_bincenter_split{ss}(2,rr) > A_zone_end && global_bincenter_split{ss}(2,rr) <= 100)
            ct_split.zoneIII(ss) = ct_split.zoneIII(ss) + 1;
        end
        
        %if both in zone I --> II
        if (global_bincenter_split{ss}(1,rr) > 0 && global_bincenter_split{ss}(1,rr) <= B_zone_end) &&...
                (global_bincenter_split{ss}(2,rr) > B_zone_end && global_bincenter_split{ss}(2,rr) <= A_zone_end)
            ct_split.zoneI_II(ss) = ct_split.zoneI_II(ss) + 1;
        end
        
        %if both in zone I --> III
        if (global_bincenter_split{ss}(1,rr) > 0 && global_bincenter_split{ss}(1,rr) <= B_zone_end) &&...
                (global_bincenter_split{ss}(2,rr) > A_zone_end && global_bincenter_split{ss}(2,rr) <= 100)
            ct_split.zoneI_III(ss) = ct_split.zoneI_III(ss) + 1;
        end
        
        %if both in zone II --> I
        if (global_bincenter_split{ss}(1,rr) > B_zone_end && global_bincenter_split{ss}(1,rr) <= A_zone_end) &&...
                (global_bincenter_split{ss}(2,rr) > 0 && global_bincenter_split{ss}(2,rr) <= B_zone_end)
            ct_split.zoneII_I(ss) = ct_split.zoneII_I(ss) + 1;
        end
        
        %if both in zone II --> III
        if (global_bincenter_split{ss}(1,rr) > B_zone_end && global_bincenter_split{ss}(1,rr) <= A_zone_end) &&...
                (global_bincenter_split{ss}(2,rr) > A_zone_end && global_bincenter_split{ss}(2,rr) <= 100)
            ct_split.zoneII_III(ss) = ct_split.zoneII_III(ss) + 1;
        end
        
        %if both in zone III --> II
        if (global_bincenter_split{ss}(1,rr) > A_zone_end && global_bincenter_split{ss}(1,rr) <= 100) &&...
                (global_bincenter_split{ss}(2,rr) > B_zone_end && global_bincenter_split{ss}(2,rr) <= A_zone_end)
            ct_split.zoneIII_II(ss) = ct_split.zoneIII_II(ss) + 1;
        end
        
        %if both in zone III --> I
        if (global_bincenter_split{ss}(1,rr) > A_zone_end && global_bincenter_split{ss}(1,rr) <= 100) &&...
                (global_bincenter_split{ss}(2,rr) > 0 && global_bincenter_split{ss}(2,rr) <= B_zone_end)
            ct_split.zoneIII_I(ss) = ct_split.zoneIII_I(ss) + 1;
        end
        
        
    end
end

%turn in to fractional count
%organize by zone with first digit being A zone
%animal x zone_switch
all_zone_counts_split = [ct_split.zoneI, ct_split.zoneII, ct_split.zoneIII,...
    ct_split.zoneI_II,ct_split.zoneII_I,...
    ct_split.zoneII_III,ct_split.zoneIII_II,...
    ct_split.zoneI_III ,ct_split.zoneIII_I];

%get sum for each animal and get fractional count
total_global_neurons_by_animal = sum(all_zone_counts_split,2);

%fractional count
frac_zones = all_zone_counts_split./total_global_neurons_by_animal;

%QC - check - adds up to 1
%sum(frac_zones,2)

%get mean and sem
mean_fraction_zone_split = mean(frac_zones,1);
sem_fraction_zone_split = std(frac_zones,0,1)./sqrt(session_nb);


%plot bar chart
% figure
% hold on
% bar(all_zone_counts)
% xticks(1:9)
% xticklabels({'I-I','II-II','III-III','I-II','I-III','II-I','II-III','III-II','III-I'})

%plot mean bar chart
figure
hold on
bar([1 2 3 4.5 5.5 7 8 9.5 10.5], mean(frac_zones,1),'FaceColor',[139, 0, 139]/255)
xticks([1 2 3 4.5 5.5 7 8 9.5 10.5])
errorbar([1 2 3 4.5 5.5 7 8 9.5 10.5],mean_fraction_zone_split,sem_fraction_zone_split,'k.')
xticklabels({'AI <--> BI','AII <--> BII','AIII <--> BIII','AI --> BII','BI --> AII',...
    'AII --> BIII','BII --> AIII','AI --> BIII','BI --> AIII'})
xtickangle(45)
%zone I,II,III switches
%sigstar({[1,2], [1,3], [2,3]})
%I vs. II switch, II vs. III switch, I vs. III
%sigstar({[4.5 5.5],[7 8],[9.5 10.5]})
yticks([0 0.1 0.2 0.3])
ylabel('Fraction of neurons');

%% Export Fig 3g data (global remapping between reward zones)

remap_prop_figs.global_zones.frac_zones = frac_zones;
remap_prop_figs.global_zones.mean_fraction_zone_split = mean_fraction_zone_split;
remap_prop_figs.global_zones.sem_fraction_zone_split = sem_fraction_zone_split;


%% Plot common histogram (Fig 3e)
figure;
hold on;
scatter(binCenter_data{1}.bin_center.final.common(1,:),binCenter_data{1}.bin_center.final.common(2,:))
%combine common centers into one matrix
common_centers = [];
for ee=1:size(path_dir,2)
    common_centers = [common_centers, binCenter_data{ee}.bin_center.final.common];
end

%use A trial centers and not mean
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

%% Common place field distribution (relative to A place field center bins - 3e)

remap_prop_figs.common_dist.common_bins_combined = common_bins_combined;
remap_prop_figs.common_dist.B_zone_end = B_zone_end;
remap_prop_figs.common_dist.A_zone_end = A_zone_end;

%% Split into animal by animal plot - display histogram for each animal
if 0
figure('Position', [2663 91 301 887])
for ee=1:size(path_dir,2)
    subplot(6,2,ee)
    hold on
    title(num2str(ee))
    ylim([0 0.5])
    if (rem(ee,2) ~= 0)
        ylabel('Normalized density');
    end
    xticks([0 100])
    xticklabels({'0','1'})
    %add number of neurons for each animal
    txt = ['n=' num2str(size(binCenter_data{ee}.bin_center.final.common(1,:),2))];
    text(60,0.4,txt)
    
    h1 = histogram(binCenter_data{ee}.bin_center.final.common(1,:),0:10:100,'Normalization','probability');
    h1.FaceColor = [139 0 139]./255;
    h1.FaceAlpha = 1;
    
    %add reward zones
    plot([B_zone_end B_zone_end], [0 0.3],'LineStyle','--', 'LineWidth',2,'Color',[220,20,60]./255)
    plot([A_zone_end A_zone_end], [0 0.3],'LineStyle','--', 'LineWidth',2,'Color',[65,105,225]./255)
    
    if ee ==10 || ee==11
        xlabel('Normalized position')
    end
end
end

%% Generate mean and sem for common neurons for supplement (not used)

for ee=1:size(path_dir,2)
    %get normalization counts
    [N(ee,:),edges] = histcounts(binCenter_data{ee}.bin_center.final.common(1,:),0:10:100,'Normalization','probability');
end

%get standard error of mean
sem_counts = std(N,0,1)./sqrt(11);

if 0
figure
hold on
yticks([0 0.05 0.1 0.15 0.2])
xticks([0 100])
xticklabels({'0','1'})

h1 = bar(5:10:100,mean(N),'BarWidth',1);
errorbar(5:10:100, mean(N),sem_counts,'LineStyle','none','Color','k' )
xlabel('Normalized position')
ylabel('Normalized density')
h1.FaceColor = [139 0 139]./255;
h1.FaceAlpha = 1;

plot([B_zone_end B_zone_end], [0 0.25],'LineStyle','--', 'LineWidth',2,'Color',[220,20,60]./255)
plot([A_zone_end A_zone_end], [0 0.25],'LineStyle','--', 'LineWidth',2,'Color',[65,105,225]./255)
end

%% Do Rayleigh test of uniformity on dataset
%convert bins to radians
common_rad_angles = circ_ang2rad(common_bins_combined(1,:)./100.*360);
[pval, z] = circ_rtest(common_rad_angles)

%% Binomial test for shift (edges included) - B-A - relative to A

%zone ends in bins
A_zone_end = 100*round(mean_rew_A_pos + reward_10cm_dist,2);
B_zone_end = 100*round(mean_rew_B_pos + reward_10cm_dist,2);

% global_bins_combined
% A_minus_B_global = global_bins_combined(1,:) - global_bins_combined(2,:);

%ecdf this 
zone1_idx = find(global_bins_combined(1,:) > 0 & global_bins_combined(1,:) <= B_zone_end);
zone2_idx = find(global_bins_combined(1,:) > B_zone_end & global_bins_combined(1,:) <= A_zone_end);
zone3_idx = find(global_bins_combined(1,:) > A_zone_end & global_bins_combined(1,:) <= 100);

%field difference B-A
%if positive - B in front of A
%if negative - B behind A
field_pos_diff = global_bins_combined(2,:) - global_bins_combined(1,:);

zone1_shifts = field_pos_diff(zone1_idx);
zone2_shifts = field_pos_diff(zone2_idx);
zone3_shifts = field_pos_diff(zone3_idx);

%fraction of neurons where B is in front of A ( 
ratio_pos_neg.I = length(find(zone1_shifts > 0))/length(zone1_idx);
ratio_pos_neg.II = length(find(zone2_shifts > 0))/length(zone2_idx);
ratio_pos_neg.III = length(find(zone3_shifts > 0))/length(zone3_idx);

%trapezoid area formula
trap_area = @(a,b,h) (a+b)*h*0.5;  

%calculate area within the trapezoid
%total zone I area
zoneI_area = 1*B_zone_end/100;

%trapezoid params
zoneI.h_upper = B_zone_end/100;
zoneI.a_upper = 1-zoneI.h_upper;
zoneI.b_upper = 1;

%fraction in each area
zoneI.upper_area_frac = trap_area(zoneI.a_upper,zoneI.b_upper,zoneI.h_upper)/zoneI_area; 
zoneI.lower_area_frac = 1-zoneI.upper_area_frac;

%total zone II area
zoneII_area = 1*(A_zone_end - B_zone_end)/100;

%trapezoid params
zoneII.h_lower = (A_zone_end - B_zone_end)/100;
zoneII.a_lower = B_zone_end/100;
zoneII.b_lower = A_zone_end/100;

%fraction in each area
zoneII.lower_area_frac = trap_area(zoneII.a_lower,zoneII.b_lower,zoneII.h_lower)/zoneII_area; 

%total zone III area
zoneIII_area = 1*(1 - A_zone_end/100);

%trapezoid params
zoneIII.h_lower = (1 - A_zone_end/100);
zoneIII.a_lower = A_zone_end/100;
zoneIII.b_lower = 1;

%fraction in each area
zoneIII.lower_area_frac =trap_area(zoneIII.a_lower,zoneIII.b_lower,zoneIII.h_lower)./zoneIII_area;

%number of successful outcomes in lower part , number of all events in that
%area (lower part B < A)

zoneI.frac_lower_actual = length(find(zone1_shifts < 0))/length(zone1_idx);
zoneI.lower_area_frac
%actual: 0.2705
%predicted (area): 0.16
%conclusion: more B before A neurons here
p_zone(1) = myBinomTest(length(find(zone1_shifts < 0)),length(zone1_idx),zoneI.lower_area_frac);

zoneII.frac_lower_actual = length(find(zone2_shifts < 0))/length(zone2_idx);
zoneII.lower_area_frac
%actual: 0.7892
%predicted (area): 0.53
%conclusion: more B before A neurons here
p_zone(2) = myBinomTest(length(find(zone2_shifts < 0)),length(zone2_idx),zoneII.lower_area_frac);

zoneIII.frac_lower_actual = length(find(zone3_shifts < 0))/length(zone3_idx);
zoneIII.lower_area_frac
%actual: 0.907
%predicted (area): 0.87
%conclusion: more B before A neurons here
p_zone(3) = myBinomTest(length(find(zone3_shifts < 0)),length(zone3_idx),zoneIII.lower_area_frac);

% pOut = myBinomTest(length(find(zone1_shifts > 0)),length(zone1_idx),0.5)
% pOut = myBinomTest(length(find(zone2_shifts > 0)),length(zone2_idx),0.5)
% pOut = myBinomTest(length(find(zone3_shifts > 0)),length(zone3_idx),0.5)

%% Split global remappers by animal
for ii=1:size(A_global_bincenter,2)
    global_bincenter_split{ii} = [A_global_bincenter{ii}; B_global_bincenter{ii}];
end

for ii=1:size(A_global_bincenter,2)
    %relative to A
    zone1_idx_split{ii} = find(global_bincenter_split{ii}(1,:) > 0 & global_bincenter_split{ii}(1,:) <= B_zone_end);
    zone2_idx_split{ii} = find(global_bincenter_split{ii}(1,:) > B_zone_end & global_bincenter_split{ii}(1,:) <= A_zone_end);
    zone3_idx_split{ii} = find(global_bincenter_split{ii}(1,:) > A_zone_end & global_bincenter_split{ii}(1,:) <= 100);
end

%bin difference
for ii=1:size(A_global_bincenter,2)
    field_pos_diff_split{ii} = global_bincenter_split{ii}(2,:) - global_bincenter_split{ii}(1,:);
end

for ii=1:size(A_global_bincenter,2)
    zone1_shifts_split{ii} = field_pos_diff_split{ii}(zone1_idx_split{ii});
    zone2_shifts_split{ii} = field_pos_diff_split{ii}(zone2_idx_split{ii});
    zone3_shifts_split{ii} = field_pos_diff_split{ii}(zone3_idx_split{ii});
end

%fraction of neurons where B is behind A (
for ii=1:size(A_global_bincenter,2)
    ratio_B_beforeA.I(ii) = length(find(zone1_shifts_split{ii} < 0))/length(zone1_idx_split{ii});
    ratio_B_beforeA.II(ii) = length(find(zone2_shifts_split{ii} < 0))/length(zone2_idx_split{ii});
    ratio_B_beforeA.III(ii) = length(find(zone3_shifts_split{ii} < 0))/length(zone3_idx_split{ii});
end

%mean(ratio_B_beforeA.I) vs. %0.16
%mean(ratio_B_beforeA.II) vs. %0.53
%mean(ratio_B_beforeA.III) vs. %0.87

%subtract expected value and do Mann Whitney U test
mean(ratio_B_beforeA.I - zoneI.lower_area_frac)

%test against median distribution of 0
p(1) = signrank(ratio_B_beforeA.I - zoneI.lower_area_frac);
p(2) = signrank(ratio_B_beforeA.II - zoneII.lower_area_frac);
p(3) = signrank(ratio_B_beforeA.III - zoneIII.lower_area_frac);

%% Plot scatter for each zone

zoneI_diff = ratio_B_beforeA.I - zoneI.lower_area_frac;
zoneII_diff = ratio_B_beforeA.II - zoneII.lower_area_frac;
zoneIII_diff = ratio_B_beforeA.III - zoneIII.lower_area_frac;

%get the means
zone_med(1) = mean(zoneI_diff);
zone_med(2) = mean(zoneII_diff);
zone_med(3) = mean(zoneIII_diff);
%get the sems
zone_sem(1) = std(zoneI_diff)./sqrt(size(zoneI_diff,2));
zone_sem(2) = std(zoneII_diff)./sqrt(size(zoneII_diff,2));
zone_sem(3) = std(zoneIII_diff)./sqrt(size(zoneIII_diff,2));


%% Difference plot from expected fraction of B before A (Fig. 3f)
figure('renderer','painters','Position',[2186 348 341 420])
hold on
xlim([0 4])
ylim([-0.5 0.5])
xticks([1 2 3])
yticks([-0.4, -0.2, 0, 0.2, 0.4])
xticklabels({'Zone I','Zone II','Zone III'})
dot_size = 14;
scatter(ones(1,size(zoneI_diff,2)),zoneI_diff,dot_size,'filled','MarkerFaceColor',[0.5 0.5 0.5])
scatter(2*ones(1,size(zoneII_diff,2)),zoneII_diff,dot_size,'filled','MarkerFaceColor',[0.5 0.5 0.5])
scatter(3*ones(1,size(zoneIII_diff,2)),zoneIII_diff,dot_size,'filled','MarkerFaceColor',[0.5 0.5 0.5])
bar([1 2 3],zone_med,'FaceColor',[139, 0, 139]/255,'EdgeColor',[0 0 0]/255,'FaceAlpha',0.3)
errorbar([1 2 3],zone_med,zone_sem,'LineStyle','none','Color', [0 0 0])
%0 dash line
%plot([0 4],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',1)

%% Export data (4f - global remapping field shift A relative to B)

remap_prop_figs.global_remap_shift.zoneI_diff = zoneI_diff;
remap_prop_figs.global_remap_shift.zoneII_diff = zoneII_diff;
remap_prop_figs.global_remap_shift.zoneIII_diff = zoneIII_diff;
remap_prop_figs.global_remap_shift.zone_med = zone_med;
remap_prop_figs.global_remap_shift.zone_sem = zone_sem;

%% Generate table with values of number of neurons in each category

zone_fraction_ct = [length(find(zone1_shifts < 0)),0,length(zone1_idx);
    length(find(zone2_shifts < 0)),0,length(zone2_idx);
    length(find(zone3_shifts < 0)),0,length(zone3_idx)];

for ii=1:3
    zone_fraction_ct(ii,2) = zone_fraction_ct(ii,3)- zone_fraction_ct(ii,1);
end


%% Plot stacked columns (not used)
if 0
figure
hold on
% ba = bar([ratio_pos_neg.I, 1-ratio_pos_neg.I; ratio_pos_neg.II, 1-ratio_pos_neg.II; ratio_pos_neg.III, 1-ratio_pos_neg.III],...
%         'stacked','FaceColor','flat')
% xticks([1 2 3])
% xticklabels({'Zone I','Zone II','Zone III'})
% legend([ba(1,1),ba(1,2)],{'A behind B','A after B'},'AutoUpdate','off')
% %plot 50% chance line
% plot([-0.2 4.5],[0.5 0.5],'k--','LineWidth',1)

%try to invert (A after B is below)    
ba = bar([1-ratio_pos_neg.I,ratio_pos_neg.I; 1-ratio_pos_neg.II,ratio_pos_neg.II; 1-ratio_pos_neg.III, ratio_pos_neg.III],...
        'stacked','FaceColor','flat');
xticks([1 2 3])
xticklabels({'Zone I','Zone II','Zone III'})
legend([ba(1,1),ba(1,2)],{'B before A','B after A'},'AutoUpdate','off')
%plot 50% chance line
%plot([-0.2 4.5],[0.5 0.5],'k--','LineWidth',1)

ylim([0 1.2])
xlim([-0.2 4.2])
ylabel('Fraction of neurons');
yticks([0 0.2 0.4 0.6 0.8 1])
%set bar colors
ba(1).CData = [0,100,0]/255;
ba(2).CData = [1,1,1]*0.8;
set(gca,'FontSize',18)
set(gca,'LineWidth',2)
end

%% Partial remapping neurons distribution (updated)

%rearrange partially remapping neurons
%combine common and partial far
for ii=1:size(binCenter_data,2)
    partial_com{ii} = binCenter_data{ii}.bin_center.final.partial_com;
    partial_far{ii} = binCenter_data{ii}.bin_center.final.partial_far;
end

%% Split by animal - combine common and partial fields

%combine the common and partial fields by animal
for ii=1:size(binCenter_data,2)
    partial_com_far_split{ii} = [partial_com{ii};partial_far{ii}];
end

%get A and B indexes for each animal
for ii=1:size(binCenter_data,2)
    partial_A_idx_split{ii} = find(~isnan(partial_com_far_split{ii}(3,:)) == 1);
    partial_B_idx_split{ii} = find(~isnan(partial_com_far_split{ii}(4,:)) == 1);
end

%get common and far fields for A or B type partial neurons
for ii=1:size(binCenter_data,2)
    partial_A_common_far_split{ii} = partial_com_far_split{ii}(1:3,partial_A_idx_split{ii});
    partial_B_common_far_split{ii} = partial_com_far_split{ii}([1,2,4],partial_B_idx_split{ii});
end

%isolate common field (A field) for each animal and partial and combined in
%row
for ii=1:size(binCenter_data,2)
    %common field (A)
    partial_A_common_far_split_input{ii}(1,:) = partial_A_common_far_split{ii}(1,:);
    %partial field (A)
    partial_A_common_far_split_input{ii}(2,:) = partial_A_common_far_split{ii}(3,:);
    
    %common field (A)
    partial_B_common_far_split_input{ii}(1,:) = partial_B_common_far_split{ii}(1,:);
    %partial field (B)
    partial_B_common_far_split_input{ii}(2,:) = partial_B_common_far_split{ii}(3,:);
end

%% Combined analysis

%first 2 rows - common field; last 2 rows - partial field; A, then B
partial_com_far_combined = [cell2mat(partial_com);cell2mat(partial_far)];

%parse into A and B partial remapping neurons (indices)
partial_A_idx = find(~isnan(partial_com_far_combined(3,:)) == 1);
partial_B_idx = find(~isnan(partial_com_far_combined(4,:)) == 1);

%get common for both fields and selective field
%first 2 rows indicate the common field
%third row is the partial field for A trials
%fourth row is the partial field for B trials
partial_A_common_far = partial_com_far_combined(1:3,partial_A_idx);
partial_B_common_far = partial_com_far_combined([1,2,4],partial_B_idx);

%prepare inputs for function below to generate histograms
%first column is mean of common position and latter column is position of
%partial field

%use the A field as common field (OK b/c partial field now collapsed to
%third row)
partial_A_common_far_input(:,1) = partial_A_common_far([1],:);
partial_A_common_far_input(:,2) = partial_A_common_far(3,:);

%use the A field as common field (OK b/c partial field now collapsed to
%third row)
partial_B_common_far_input(:,1) = partial_B_common_far([1],:);
partial_B_common_far_input(:,2) = partial_B_common_far(3,:);

%% Bin as a function of A center and get boxplots for distribution
%place in bin for common A field
%get bin indices
%10 even spaced bins
%uncomment if problem
%partial_bin_idx.A = discretize(partial_A_common_far_input(:,1),[0:10:100]);
%partial_bin_idx.B = discretize(partial_B_common_far_input(:,1),[0:10:100]);

%3 behaviorally spaced bins (1-35, 35-75, 75-100)
%this will assign values to one of three bins defined in parameters using
%mean common position
partial_bin_idx_3.A = discretize(partial_A_common_far_input(:,1),[0,B_zone_end,A_zone_end,100]);
partial_bin_idx_3.B = discretize(partial_B_common_far_input(:,1),[0,B_zone_end,A_zone_end,100]);

%for each bin (3 behavior)
for bb=1:3
    count_3.A{bb} = partial_A_common_far_input(find(partial_bin_idx_3.A == bb),2);
    count_3.B{bb} = partial_B_common_far_input(find(partial_bin_idx_3.B == bb),2);
end

%merge both cell in alternate manner (input into bar plotter)
merge_hist_counts_3(1:2:5) = count_3.A;
merge_hist_counts_3(2:2:6) = count_3.B;

%% Extract partial neuron idx's that are examples of each zone for supplement

%extract all from each zone for each animal and export
for ii=1:size(binCenter_data,2)
    partial_idx_by_animal{ii}.A = remapping_idxs{ii}.remapping_corr_idx.final.partial(partial_A_idx_split{ii});
    partial_idx_by_animal{ii}.B = remapping_idxs{ii}.remapping_corr_idx.final.partial(partial_B_idx_split{ii});
end

%% Bin as a function of A center and get boxplots for distribution by animal/session
%get bin indices
%10 even spaced bins

%split by common field
% for ii=1:size(binCenter_data,2)
%     partial_bin_idx_split{ii}.A = discretize(partial_A_common_far_split_input{ii}(1,:)',[0:10:100]);
%     partial_bin_idx_split{ii}.B = discretize(partial_B_common_far_split_input{ii}(1,:)',[0:10:100]);
% end

%3 behaviorally spaced bins (1-35, 35-75, 75-100)
%this will assign values to one of three bins defined in parameters using
%A common position
for ii=1:size(binCenter_data,2)
    partial_bin_idx_3_split{ii}.A = discretize(partial_A_common_far_split_input{ii}(1,:)',[0,B_zone_end,A_zone_end,100]);
    partial_bin_idx_3_split{ii}.B = discretize(partial_B_common_far_split_input{ii}(1,:)',[0,B_zone_end,A_zone_end,100]);
end

%for each bin (3 behavior)
for ii=1:size(binCenter_data,2)
    for bb=1:3
        count_3_split{ii}.A{bb} = partial_A_common_far_split_input{ii}(2,find(partial_bin_idx_3_split{ii}.A == bb))';
        count_3_split{ii}.B{bb} = partial_B_common_far_split_input{ii}(2,find(partial_bin_idx_3_split{ii}.B == bb))';
        
        %place idx extraction here using same approach (export these for
        %animal plots)
        partial_idx_by_animal_zone{ii}.A.zone{bb} = partial_idx_by_animal{ii}.A(find(partial_bin_idx_3_split{ii}.A == bb));
        partial_idx_by_animal_zone{ii}.B.zone{bb} = partial_idx_by_animal{ii}.B(find(partial_bin_idx_3_split{ii}.B == bb));
        
    end
end

%merge both cell in alternate manner (input into bar plotter) for each
%animal
for ii=1:size(binCenter_data,2)
    merge_hist_counts_3_split{ii}(1:2:5) = count_3_split{ii}.A;
    merge_hist_counts_3_split{ii}(2:2:6) = count_3_split{ii}.B;
end

%QC check if split adds to total number of partial neurons (checks out)
% for ii=1:size(binCenter_data,2)
%     element_partial{ii} = cell2mat(merge_hist_counts_3_split{ii}');
% end
% 
% nb_partial_remappers_check = length(cell2mat(element_partial'));
% isequal(nb_partial_remappers_check, size(partial_com_far_combined,2))

%% Statistics
[h,p] = ranksum(count_3.A{1},count_3.B{1})
[h,p] = ranksum(count_3.A{2},count_3.B{2})
[h,p] = ranksum(count_3.A{3},count_3.B{3})


%% Boxplot of distributions relative to mean of centroid in each zone

merge_new = [merge_hist_counts_3([1,3,5]);merge_hist_counts_3([2,4,6])];

xlab={'1','2','3'} ;
col=[220,20,60, 255;
65,105,225, 255;
]; 
%0, 0, 255, 200]; 
col=col/255;

if 0
f=figure 
hold on
%plot zone separator lines
%[f2,x,group,positions,labelpos] =  multiple_boxplot(merge_new',xlab,{'A','B'},col');
%overlay boxplot to add median line 
%z= boxplot(x,group, 'positions', positions);
lines = findobj(z, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
set(lines, 'LineWidth',2)
ylim([0 100])
xticks([1.3750, 2.1250, 2.875])
xticklabels({'Zone I', 'Zone II', 'Zone III'})
ylabel('Place field position')
%yticks([0 20 40 60 80 100])
%yticklabels({'0', '0.2', '0.4','0.6','0.8','1'})
plot([1.7500 ,1.7500],[-10 110],'k--','LineWidth',1.5)
plot([2.500 ,2.500],[-10 110],'k--','LineWidth',1.5)
set(gca,'FontSize',16)

end

%% Make plot of of scatter and bars to show bimodal distribution of A lap neurons

%get colormap
paper_cmap = return_paper_colormap();

figure('Position', [2335 351 915 597])
hold on
xlim([0 10])
marker_size = 18;
ylabel('Partial place field position along track')
yticks([0 20 40 60 80 100])
yticklabels({'0','0.2','0.4','0.6','0.8','1'})

xticks([1.5 4.5 7.5])
xticklabels({'Zone I','Zone II','Zone III'})

%AI
bar(1, median(merge_new{1,1}),'FaceColor','none','EdgeColor',paper_cmap(1,:),'LineWidth',1.5)
scatter(ones(1,size(merge_new{1,1},1)),merge_new{1,1},marker_size,'filled','MarkerFaceColor',paper_cmap(1,:))
%BI
bar(2, median(merge_new{2,1}),'FaceColor','none','EdgeColor',paper_cmap(2,:),'LineWidth',1.5)
scatter(2*ones(1,size(merge_new{2,1},1)),merge_new{2,1},marker_size,'filled','MarkerFaceColor',paper_cmap(2,:))

%AII
bar(4, median(merge_new{1,2}),'FaceColor','none','EdgeColor',paper_cmap(1,:),'LineWidth',1.5)
scatter(4*ones(1,size(merge_new{1,2},1)),merge_new{1,2},marker_size,'filled','MarkerFaceColor',paper_cmap(1,:))
%BII
bar(5, median(merge_new{2,2}),'FaceColor','none','EdgeColor',paper_cmap(2,:),'LineWidth',1.5)
scatter(5*ones(1,size(merge_new{2,2},1)),merge_new{2,2},marker_size,'filled','MarkerFaceColor',paper_cmap(2,:))

%AIII
bar(7, median(merge_new{1,3}),'FaceColor','none','EdgeColor',paper_cmap(1,:),'LineWidth',1.5)
scatter(7*ones(1,size(merge_new{1,3},1)),merge_new{1,3},marker_size,'filled','MarkerFaceColor',paper_cmap(1,:))
%BIII
bar(8, median(merge_new{2,3}),'FaceColor','none','EdgeColor',paper_cmap(2,:),'LineWidth',1.5)
scatter(8*ones(1,size(merge_new{2,3},1)),merge_new{2,3},marker_size,'filled','MarkerFaceColor',paper_cmap(2,:))

%add significance lines
% sigstar([1,2])
% sigstar([4,5])
% sigstar([7,8])

set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)
%set(gca,'XTick',[])

%% Boxplot of partial remapping distributions as subplot for each animal

%merge for each animal
for ii=1:size(binCenter_data,2)
    merge_new_split{ii} = [merge_hist_counts_3_split{ii}([1,3,5]);merge_hist_counts_3_split{ii}([2,4,6])];
end

%count each type cell type in each subgroup for each animal
for ii=1:size(binCenter_data,2)
    %place each count in to array
    partial_ct_by_animal(:,:,ii) = cellfun(@(x) size(x,1),merge_new_split{ii},'UniformOutput',true);
end

%replace empty cells with nans
for ii=1:size(binCenter_data,2)
    %find index of empty cells
    empty_cells{ii} = find(cellfun(@isempty,merge_new_split{ii}) ==1);
    %if index found
    if ~isempty(empty_cells{ii})
        merge_new_split{ii}(empty_cells{ii}) = {nan};
    end
end

%% For each animal, boxplots of A vs. B counts
if 1
    f=figure
    for ii=1:size(binCenter_data,2)
        subplot(6,2,ii)
        hold on
        %plot zone separator lines
        [f2,x,group,positions,labelpos] =  multiple_boxplot_subplot(merge_new_split{ii}',xlab,{'A','B'},col');
        
        title(num2str(ii))
        %overlay boxplot to add median line
        z= boxplot(x,group, 'positions', positions);
        lines = findobj(z, 'type', 'line', 'Tag', 'Median');
        set(lines, 'Color', 'k');
        set(lines, 'LineWidth',2)
        ylim([0 100])
        xticks([1.3750, 2.1250, 2.875])
        xticklabels({'Zone I', 'Zone II', 'Zone III'})
        ylabel('Place field position')
        %yticks([0 20 40 60 80 100])
        %yticklabels({'0', '0.2', '0.4','0.6','0.8','1'})
        plot([1.7500 ,1.7500],[-10 110],'k--','LineWidth',1.5)
        plot([2.500 ,2.500],[-10 110],'k--','LineWidth',1.5)
        %set(gca,'FontSize',16)
        
        %add number count to display
        txtA_1 = [num2str(partial_ct_by_animal(1,1,ii))];
        text(60,0.8,txtA_1)
        
    end
end


%% Distribution of partially remapping fields for A and B plot and KS test for significance (Figure 3H)

figure('Position',[2200 317 508 433])
hold on
title('Distribution of partially remapping fields')
[fA,xA] = ecdf(partial_A_common_far_input(:,2));
[fB,xB] = ecdf(partial_B_common_far_input(:,2));
ap = stairs(xA,fA,'LineWidth',2,'Color',[65,105,225]/255);
bp = stairs(xB,fB,'LineWidth',2,'Color',[220,20,60]/255);
xlabel(gca,'Normalized positon')
ylabel('Cumulative probability')

xticks([0 100])
xticklabels({'0','1'})
set(gca,'LineWidth',1.5)
set(gca,'FontSize',16)
legend([ap bp],{'A','B'},'Location','northwest','AutoUpdate','off')

plot([B_zone_end B_zone_end], [0 1],'--','Color',[220,20,60]/255, 'LineWidth',2)
plot([A_zone_end A_zone_end], [0 1],'--','Color',[65,105,225]/255, 'LineWidth',2)

%% Export Fig. 3h data partial field remap

remap_prop_figs.partial_remap.partial_A_common_far_input = partial_A_common_far_input;
remap_prop_figs.partial_remap.partial_B_common_far_input = partial_B_common_far_input;
remap_prop_figs.partial_remap.B_zone_end = B_zone_end;
remap_prop_figs.partial_remap.A_zone_end = A_zone_end;

%% Plot cumulative density function (CDF) for each animal (not used)
if 0
figure
title('Distribution of partially remapping fields for each animal')
for ii=1:11
    subplot(6,2,ii)
    hold on
    title(num2str(ii))
    [fA,xA] = ecdf(partial_A_common_far_split_input{ii}(2,:));
    [fB,xB] = ecdf(partial_B_common_far_split_input{ii}(2,:));
    xlim([0 100])
    ap = stairs(xA,fA,'LineWidth',2,'Color',[65,105,225]/255);
    bp = stairs(xB,fB,'LineWidth',2,'Color',[220,20,60]/255);
    %xlabel(gca,'Normalized position')
    %ylabel('Cumulative probability')
    
    xticks([0 100])
    xticklabels({'0','1'})
    
    %add label with neuron count
    txtA = ['n(A)=' num2str(size(partial_A_common_far_split_input{ii}(2,:),2))];
    txtB = ['n(B)=' num2str(size(partial_B_common_far_split_input{ii}(2,:),2))];
    text(68,0.4,txtA)
    text(68,0.2,txtB)
    
    if (rem(ii,2) ~= 0)
        ylabel('Cumulative probability');
    end
    
    %add A vs. B legend on last animal
    if ii==11
    legend([ap bp], {'A','B'},'location','northwest')
    end
    
    if ii ==10 || ii==11
        xlabel('Normalized position')
    end
end
end

%% Make histograms for partial fields (skip for now - come back if necessary)

figure('Position', [2663 91 301 887])
for ee=1:size(path_dir,2)
    subplot(6,2,ee)
    hold on
    title(num2str(ee))
    ylim([0 0.5])
    if (rem(ee,2) ~= 0)
        ylabel('Normalized density');
    end
    xticks([0 100])
    xticklabels({'0','1'})
    %add number of neurons for each animal
    txt = ['n=' num2str(size(binCenter_data{ee}.bin_center.final.common(1,:),2))];
    text(60,0.4,txt)
    
    h1 = histogram(binCenter_data{ee}.bin_center.final.common(1,:),0:10:100,'Normalization','probability');
    h1.FaceColor = [139 0 139]./255;
    h1.FaceAlpha = 1;
    
    %add reward zones
    plot([B_zone_end B_zone_end], [0 0.3],'LineStyle','--', 'LineWidth',2,'Color',[220,20,60]./255)
    plot([A_zone_end A_zone_end], [0 0.3],'LineStyle','--', 'LineWidth',2,'Color',[65,105,225]./255)
    
    if ee ==10 || ee==11
        xlabel('Normalized position')
    end
end


%% Plot scatter for each animal
%common bin on x axis
%partial bin on y axis


%% Statistics (combined dataset)

[h,p,ksstat] = kstest2(partial_A_common_far_input(:,2), partial_B_common_far_input(:,2));

%number of neurons used in kstest
nb_partial_ks_test = [size(partial_A_common_far_input(:,2),1), size(partial_B_common_far_input(:,2),1)];

%% Minaturized normalized histogram for each set of partial remapping neurons (Fig 3h insets)
%cdf
figure('Position',[1199 400 208 420])
subplot(2,1,1)
hold on;
ylim([0 0.2])
yticks([0 0.1 0.2])
xticks([0 50 100])
xticklabels({'0','0.5','1'})
ylabel('Norm. density')
set(gca,'FontSize',16)
histogram(partial_A_common_far_input(:,2),[0:10:100],'Normalization','probability','FaceColor', [65,105,225]/255,'FaceAlpha',1)

plot([B_zone_end B_zone_end], [0 0.2],'--','Color',[220,20,60]/255, 'LineWidth',2)
plot([A_zone_end A_zone_end], [0 0.2],'--','Color',[65,105,225]/255, 'LineWidth',2)

subplot(2,1,2)
hold on
ylim([0 0.2])
yticks([0 0.1 0.2])
xticks([0 50 100])
xticklabels({'0','0.5','1'})
ylabel('Norm. density')
xlabel('Normalized ')
set(gca,'FontSize',16)
histogram(partial_B_common_far_input(:,2),[0:10:100],'Normalization','probability','FaceColor',[220,20,60]/255,'FaceAlpha',1)

plot([B_zone_end B_zone_end], [0 0.2],'--','Color',[220,20,60]/255, 'LineWidth',2)
plot([A_zone_end A_zone_end], [0 0.2],'--','Color',[65,105,225]/255, 'LineWidth',2)

%% Export Fig. 3h data partial field remap (insets)

remap_prop_figs.partial_remap.insets.partial_A_common_far_input = partial_A_common_far_input;
remap_prop_figs.partial_remap.insets.partial_B_common_far_input = partial_B_common_far_input;
remap_prop_figs.partial_remap.insets.B_zone_end = B_zone_end;
remap_prop_figs.partial_remap.insets.A_zone_end = A_zone_end;

%% OLD

%boxplot 2 wrapper(plot distributions of diff sizes
% col=@(x)reshape(x,numel(x),1);
% boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});


% 
% merge_far_binCenters = [];
% figure
% hold on
% title('near + far ')
% xlabel('Center of Place field A')
% ylabel('Center of Place field B')
% for ee=1:size(path_dir,2)
%     scatter(binCenter_data{ee}.bin_center.global_near(1,:),binCenter_data{ee}.bin_center.global_near(2,:),'g')
%     %scatter(binCenter_data{ee}.bin_center.global_far(1,:),binCenter_data{ee}.bin_center.global_far(2,:),'r')
% %merge far centroids
% merge_far_binCenters = [merge_far_binCenters,binCenter_data{ee}.bin_center.global_far,binCenter_data{ee}.bin_center.global_near];
% end
% 
% plot([0 100],[0 100],'k--')

%tile histogram
% figure;
% hold on
% nbins = 15;
% xedges = 0:5:100;
% xedges = [0,33,66, 100];
% yedges = xedges;
% xlim([0 100]);
% ylim([0 100])
% h = histogram2(merge_far_binCenters(1,:),merge_far_binCenters(2,:),xedges,yedges,'DisplayStyle','tile','ShowEmptyBins','on');
% %h = histogram2(merge_far_binCenters(1,:),merge_far_binCenters(2,:),nbins,'DisplayStyle','tile','ShowEmptyBins','on');
% colormap('hot')
% colorbar
% 
% %find all y's below x values
% less_than_unity = merge_far_binCenters(1,:)>merge_far_binCenters(2,:);
% figure
% hold on
% histogram2(merge_far_binCenters(1,less_than_unity),merge_far_binCenters(2,less_than_unity),xedges,yedges,'DisplayStyle','tile','ShowEmptyBins','on')
% 
% %only low track PV correlated animals
% figure
% hold on
% title('Low PV correlated animal - far global')
% xlabel('Center of Place field A')
% ylabel('Center of Place field B')
% for ee=options.lowPVcorr
%     scatter(binCenter_data{ee}.bin_center.global_far(1,:),binCenter_data{ee}.bin_center.global_far(2,:),'m')
% end
% plot([0 100],[0 100],'k--')
% 
% %plot partial remappers
% figure
% hold on
% for ee=1:size(path_dir,2)
% for rr=1:size(binCenter_data{ee}.bin_center.partial_com,2)
%     if isnan(binCenter_data{ee}.bin_center.partial_far(1,rr))
%         scatter(mean(binCenter_data{ee}.bin_center.partial_com(:,rr)),binCenter_data{ee}.bin_center.partial_far(2,rr),'r')
%     else
%         scatter(mean(binCenter_data{ee}.bin_center.partial_com(:,rr)),binCenter_data{ee}.bin_center.partial_far(1,rr),'b')
%     end
% end
% end

%
%for each bin (10 even)
% for bb=1:10
%     count.A{bb} = partial.com_parA(find(partial_bin_idx.A == bb),2)
%     count.B{bb} = partial.com_parB(find(partial_bin_idx.B == bb),2)
% end


%merge both cell in alternate manner
% merge_hist_counts(1:2:20) = count.A;
% merge_hist_counts(2:2:21) = count.B;


% figure;
% %Distribution of partial indexes  
% hold on
% boxplot2(merge_hist_counts_3)
% 
% %use matlab boxplot to generate groups
% grouping_box{1} = repmat(1,size(count_3.A{1},1),1)
% grouping_box{2} = repmat(2,size(count_3.B{1},1),1)
% grouping_box{3} = repmat(3,size(count_3.A{2},1),1)
% grouping_box{4} = repmat(4,size(count_3.B{2},1),1)
% grouping_box{5} = repmat(5,size(count_3.A{3},1),1)
% grouping_box{6} = repmat(6,size(count_3.B{3},1),1)
% %catx=repmat({'I' 'II','III'},1,100);
% c=cell2mat(grouping_box')';
% c(find(c==1 | c==3 | c==5)) = 1;
% c(find(c==2 |c==4 | c==6)) = 2

% g=gramm('x',cell2mat(grouping_box'),'y',cell2mat(merge_hist_counts_3'),'color',c);
% g.set_color_options('map',[65,105,225; 220,20,60; 255,0,255; 0 0 0]/255);

%%%
% Plot raw data as points
% g.geom_point();
% g.stat_boxplot('width',1,'dodge',2);
% g.geom_vline('xintercept',[2.5,4.5])
% % Set appropriate names for legends
% g.set_names('column','Origin','x','Zone','y','Place field center','color','Trial type');
% %%%
% % Set figure title
% %g.set_title('Fuel economy of new cars between 1970 and 1982');
% %%%
% % Do the actual drawing
% figure('Position',[100 100 800 400]);
% g.draw();

%{

%% 2 subplots - store common vs. partial center for each A or B find
storeIdx =0;
figure
subplot(1,2,1)
hold on
for ee=1:size(path_dir,2)
    for rr=1:size(binCenter_data{ee}.bin_center.partial_com,2)
        if isnan(binCenter_data{ee}.bin_center.partial_far(1,rr))
            storeIdx = storeIdx + 1;
            scatter(mean(binCenter_data{ee}.bin_center.partial_com(:,rr)),binCenter_data{ee}.bin_center.partial_far(2,rr),'r')
            %store B mean common vs. B partial in matrix
%all animals
            partial.com_parB(storeIdx,:) = [mean(binCenter_data{ee}.bin_center.partial_com(:,rr)), binCenter_data{ee}.bin_center.partial_far(2,rr)];
        else
            %scatter(mean(binCenter_data{ee}.bin_center.partial_com(:,rr)),binCenter_data{ee}.bin_center.partial_far(1,rr),'b')
        end
    end
end
storeIdx =0;

subplot(1,2,2)
hold on
for ee=1:size(path_dir,2)
for rr=1:size(binCenter_data{ee}.bin_center.partial_com,2)
    if isnan(binCenter_data{ee}.bin_center.partial_far(1,rr))
        %scatter(mean(binCenter_data{ee}.bin_center.partial_com(:,rr)),binCenter_data{ee}.bin_center.partial_far(2,rr),'r')
    else
storeIdx = storeIdx + 1;
        scatter(mean(binCenter_data{ee}.bin_center.partial_com(:,rr)),binCenter_data{ee}.bin_center.partial_far(1,rr),'b')
            %store A mean common vs. A partial in matrix
%all animals            
partial.com_parA(storeIdx,:) = [mean(binCenter_data{ee}.bin_center.partial_com(:,rr)), binCenter_data{ee}.bin_center.partial_far(1,rr)];
    end
end
end

%% Convert partial mappings to negative of difference and plot histogram
diff_par.A = (-1)*(partial.com_parA(:,1) -  partial.com_parA(:,2));
diff_par.B = (-1)*(partial.com_parB(:,1) -  partial.com_parB(:,2));

A_fwd = find(diff_par.A > 0);
A_rew = find(diff_par.A < 0);
[h,p] = kstest2((-1)*diff_par.A(A_fwd), diff_par.A(A_rew))

B_fwd = find(diff_par.B > 0);
B_rew = find(diff_par.B < 0);
[h,p] = kstest2((-1)*diff_par.B(B_fwd), diff_par.B(B_rew))


%test difference of number of neurons ahead of centroid
kstest2(diff_par.A,diff_par.A);

%Ks test of difference between the 2 distribtuions
[h,p] = kstest2(diff_par.A,diff_par.B);

%ks test on distribution of partial field centers
[h,p] =kstest2(partial.com_parA(:,2),partial.com_parB(:,2))



%% Distance relative to PF center position
figure;
hold on
title('')
xlabel('Distance of field from mean of place field center')
xlim([-100 100])
histogram(diff_par.A,-100:10:100)
histogram(diff_par.B,-100:10:100)
plot([0 0],[0 35],'k--')

figure;
hold on
title('')
xlabel('')
xlim([0 100])
ylim([0 40])
histogram(partial.com_parA(:,2),0:5:100)
histogram(partial.com_parB(:,2),0:5:100)
%histogram(diff_par.B,-100:10:100)
%plot([0 0],[0 35],'k--')

%%
%% For partial, scatter plot of centroids - sort my mean of common place field center

figure
hold on
for ss=1:11
    
    [~,Isort_com] = sort(mean(binCenter_data{ss}.bin_center.partial_com),'descend');
    %sort columns by common mean center
    bin_center.partial_com_sort = binCenter_data{ss}.bin_center.partial_com(:,Isort_com);
    bin_center.partial_far_sort = binCenter_data{ss}.bin_center.partial_far(:,Isort_com);
    
    
    %plot common
    hold on
    for rr=1:size(binCenter_data{ss}.bin_center.partial_com,2)
        scatter(bin_center.partial_com_sort(1,rr),rr,'g')
        scatter(bin_center.partial_com_sort(2,rr),rr,'g')
    end
    
    %plot partial
    for rr=1:size(binCenter_data{ss}.bin_center.partial_com,2)
        if ~isnan(bin_center.partial_far_sort(1,rr))
            scatter(bin_center.partial_far_sort(1,rr),rr,'b')
        else
            scatter(bin_center.partial_far_sort(2,rr),rr,'r')
        end
    end
end

figure;
for ss=1:11
%sort global far by B trials
[~,Isort_glo_far_B] = sort(binCenter_data{1}.bin_center.global_far(2,:),'descend');
[~,Isort_glo_far_A] = sort(binCenter_data{1}.bin_center.global_far(1,:),'descend');
%sort columns by common mean center
bin_center.global_far_sortB = binCenter_data{1}.bin_center.global_far(:,Isort_glo_far_B);
bin_center.global_far_sortA = binCenter_data{1}.bin_center.global_far(:,Isort_glo_far_A);
%bin_center.partial_far_sort =bin_center.partial_far(:,Isort_com);

%Plot global far

subplot(1,2,1)
hold on;
for rr=1:size(binCenter_data{1}.bin_center.global_far,2)
    scatter(bin_center.global_far_sortA(1,rr),rr,'b')
    scatter(bin_center.global_far_sortA(2,rr),rr,'r')
end
subplot(1,2,2)
hold on;
for rr=1:size(binCenter_data{1}.bin_center.global_far,2)
    scatter(bin_center.global_far_sortB(1,rr),rr,'b')
    scatter(bin_center.global_far_sortB(2,rr),rr,'r')
end
end

%% Make one common, global far, global near partial matrix
global_far.A = [];
global_far.B = [];
global_near.A = [];
global_near.B = [];
common = [];
partial_common = [];
rate = [];

for ee=1:11
    global_far.A = [global_far.A, binCenter_data{ee}.bin_center.global_far(1,:)];
    global_far.B = [global_far.B, binCenter_data{ee}.bin_center.global_far(2,:)];
    common = [common, binCenter_data{ee}.bin_center.common];
    partial_common = [partial_common,mean(binCenter_data{ee}.bin_center.partial_com,1) ];
    rate = [rate, binCenter_data{ee}.bin_center.rate];
    global_near.A = [global_near.A, binCenter_data{ee}.bin_center.global_near(1,:)]
    global_near.B = [global_near.B, binCenter_data{ee}.bin_center.global_near(2,:)]
end

%split by animal into separate cell
for ee=1:11
    global_near_each{ee} = binCenter_data{ee}.bin_center.global_near;
end


%combine global fars
global_far_comb = [global_far.A; global_far.B];

bin_space = 10;
ylim_def = 0.3;

figure
hold on
plot(global_far_comb)

figure
hold on
%scatter(global_near.A,global_near.B)
%unity line
axis square
plot([0 100], [0 100],'k-');
less_than_unity = global_near.A >= global_near.B;
greater_than_unity = global_near.A < global_near.B;
scatter(global_near.A(less_than_unity),global_near.B(less_than_unity),14,'MarkerFaceColor','b','MarkerEdgeColor','b') 
scatter(global_near.A(greater_than_unity),global_near.B(greater_than_unity),14,'MarkerFaceColor','r','MarkerEdgeColor','r')

%fit line
x = global_near.A';                                        % Create Data
y = global_near.B';               % Create Data
mdl = fitlm(x, y);                                  % Fit Data
B = mdl.Coefficients.Estimate;                      % Coefficients
CI = coefCI(mdl);                                   % Coefficient Confidence Intervals
[Ypred,YCI] = predict(mdl, x);                      % Fitted Regression Line & Confidence Intervals

Xnew = linspace(min(x), max(x), 1000);
[~,YCI] = predict(mdl, Xnew'); 

[B,BINT,R,RINT,STATS] = regress( global_near.B',[global_near.A',ones(size(global_near.B,2),1)] )

plot(Xnew, YCI, '--k')
plot(x, Ypred,'-r')
hold off
%grid



%% 6 zones - not more informative than 3 zones
if 0
    %ecdf this - break into 6 zones
    zone1_idx = find(global_near.A > 5 & global_near.A <=20);
    zone2_idx = find(global_near.A > 20 & global_near.A <=35);
    zone3_idx = find(global_near.A > 35 & global_near.A <= 60);
    zone4_idx = find(global_near.A > 60 & global_near.A <= 75);
    zone5_idx = find(global_near.A > 75 & global_near.A <= 85);
    zone6_idx = find(global_near.A > 85 & global_near.A <= 95);
    
    %field difference B-A
    field_pos_diff = global_near.B - global_near.A;
    
    zone1_shifts = field_pos_diff(zone1_idx);
    zone2_shifts = field_pos_diff(zone2_idx);
    zone3_shifts = field_pos_diff(zone3_idx);
    zone4_shifts = field_pos_diff(zone4_idx);
    zone5_shifts = field_pos_diff(zone5_idx);
    zone6_shifts = field_pos_diff(zone6_idx);
    
    ratio_pos_neg.I = length(find(zone1_shifts > 0))/length(zone1_idx)
    ratio_pos_neg.II = length(find(zone2_shifts > 0))/length(zone2_idx)
    ratio_pos_neg.III = length(find(zone3_shifts > 0))/length(zone3_idx)
    ratio_pos_neg.IV = length(find(zone4_shifts > 0))/length(zone4_idx)
    ratio_pos_neg.V = length(find(zone5_shifts > 0))/length(zone5_idx)
    ratio_pos_neg.VI = length(find(zone6_shifts > 0))/length(zone6_idx)
    
    pOut = myBinomTest(length(find(zone1_shifts > 0)),length(zone1_idx),0.5)
    pOut = myBinomTest(length(find(zone2_shifts > 0)),length(zone2_idx),0.5)
    pOut = myBinomTest(length(find(zone3_shifts > 0)),length(zone3_idx),0.5)
    pOut = myBinomTest(length(find(zone4_shifts > 0)),length(zone4_idx),0.5)
    pOut = myBinomTest(length(find(zone5_shifts > 0)),length(zone5_idx),0.5)
    pOut = myBinomTest(length(find(zone6_shifts > 0)),length(zone6_idx),0.5)
end

%% Classify global remappers into respective zones
%I-II
%II-III
%I-III
%ecdf this
zone_global_idx(1).A = find(global_far.A > 0 & global_far.A <= 35);
zone_global_idx(2).A = find(global_far.A > 35 & global_far.A <= 75);
zone_global_idx(3).A = find(global_far.A > 75 & global_far.A <= 100);

zone_global_idx(1).B = find(global_far.B > 0 & global_far.B <= 35);
zone_global_idx(2).B = find(global_far.B > 35 & global_far.B <= 75);
zone_global_idx(3).B = find(global_far.B > 75 & global_far.B <= 100);

%merge into one matrix and assign zones
merge_zone_matrix = zeros(2,size(global_far.A,2));

merge_zone_matrix(1,zone_global_idx(1).A) = 1;
merge_zone_matrix(1,zone_global_idx(2).A) = 2;
merge_zone_matrix(1,zone_global_idx(3).A) = 3;

merge_zone_matrix(2,zone_global_idx(1).B) = 1;
merge_zone_matrix(2,zone_global_idx(2).B) = 2;
merge_zone_matrix(2,zone_global_idx(3).B) = 3;

%count each type of transition
%I-II, II-I
one_two = (merge_zone_matrix(1,:) == 1) & (merge_zone_matrix(2,:) == 2);
two_one = (merge_zone_matrix(2,:) == 1) & (merge_zone_matrix(1,:) == 2);
one_two_all = one_two | two_one;

%1-III, III-I
one_three = (merge_zone_matrix(1,:) == 1) & (merge_zone_matrix(2,:) == 3);
three_one = (merge_zone_matrix(2,:) == 1) & (merge_zone_matrix(1,:) == 3);
one_three_all = one_three | three_one;

%II-III,III-II
two_three = (merge_zone_matrix(1,:) == 2) & (merge_zone_matrix(2,:) == 3);
three_two = (merge_zone_matrix(2,:) == 2) & (merge_zone_matrix(1,:) == 3);
two_three_all = two_three | three_two;

%self I-I, II-II, III-III
self_one = (merge_zone_matrix(1,:) == 1) & (merge_zone_matrix(2,:) == 1);
self_two = (merge_zone_matrix(1,:) == 2) & (merge_zone_matrix(2,:) == 2);
self_three = (merge_zone_matrix(1,:) == 3) & (merge_zone_matrix(2,:) == 3);

self_all = ((self_one | self_two) | self_three);

zone_count = [length(find(self_all == 1)), length(find(one_two_all == 1)), length(find(two_three_all == 1)),length(find(one_three_all == 1))]

%set up count for each transition (I-II, II-III, I-III)
transition_stacked = [length(find(one_two ==1)),length(find(two_one ==1));...
                    length(find(two_three ==1)),length(find(three_two ==1));...
                    length(find(one_three ==1)),length(find(three_one ==1))];

group_count = sum(transition_stacked,2);


figure;
hold on
ba2 = bar(transition_stacked./sum(transition_stacked,2),'stacked','FaceColor','flat');
xticks([1 2 3])
xticklabels({'Early <-> Mid','Mid <-> Late','Early <-> Late'})
ylim([0 1.2])
xlim([-0.2 4.2])
ylabel('Fraction of neurons');
yticks([0 0.2 0.4 0.6 0.8 1])
xtickangle(45);
%set bar colors
ba2(1).CData = [65,105,225]/255;
ba2(2).CData = [220,20,60]/255;
set(gca,'FontSize',18)
set(gca,'LineWidth',2)

%binomial test for transitions
pOut = myBinomTest(transition_stacked(1,1),group_count(1),0.5)
pOut = myBinomTest(transition_stacked(2,1),group_count(2),0.5)
pOut = myBinomTest(transition_stacked(3,1),group_count(3),0.5)

%% Common distribution 

%account for edge when calculating the mean of the common centers 
%find neurons with difference greater than 
split_common_idxs = find(diff(common,1) > 20 | diff(common,1) < -20);

%find min idx
%add 100 to min
%take mean
%subtract 100
for ii=1:size(split_common_idxs,2)
    [val_temp, idx_temp] = min(common(:,split_common_idxs(ii)));
    if idx_temp == 2
        upd_common_val(ii) = ceil(((val_temp+100) + common(1,split_common_idxs(ii)))/2 - 100);
    elseif idx_temp == 1
        upd_common_val(ii) = ceil(((val_temp+100) + common(2,split_common_idxs(ii)))/2 - 100);
    end
end
%change 0 outputs to 100 bin
upd_common_val((upd_common_val == 0)) = 100;

%take the mean and update values
mean_common = mean(common,1);
mean_common(split_common_idxs) = upd_common_val;


figure
hold on
ylabel('Normalized density')
xlabel('Track position')
histogram(mean_common,[0:10:100],'Normalization','probability','FaceColor',[139, 0, 139]/255)
set(gca,'FontSize',18)
set(gca,'LineWidth',2)

%test against uniform distribition (count of common neurons in each bin)
[h,p]= chi2gof(mean_common, 'Edges', [0:10:100],'Expected',ones(1,10)*(size(mean_common,2)/10))

%% Plot scatter as a function of place field center difference - account for edge effects

%field difference B-A
field_pos_diff = global_near.B - global_near.A;

%global_near.B(track_fwd_idx(1)) - global_near.A(track_fwd_idx(1))

%account for edge effect on difference calculation
%shifting behind track edge
track_rev_idx = find(field_pos_diff >50);
%shifting forward of track edge
track_fwd_idx = find(field_pos_diff <-50);

%account for edge effect
field_pos_diff(track_rev_idx) = field_pos_diff(track_rev_idx)-100;
field_pos_diff(track_fwd_idx) = field_pos_diff(track_fwd_idx)+100;

%function of mean position 
%mean_pos = mean([global_near.A; global_near.B],1); 

%convert to normalized space, assign as input variables
x = (global_near.A./100)';
y = (field_pos_diff./100)';

figure;
hold on
scatter(x,  y',14,'MarkerFaceColor',[139,0,139]/255,'MarkerEdgeColor',[139,0,139]/255)
axis square

%fit line
%x = global_near.A';                                        % Create Data
%y = field_pos_diff';               % Create Data
mdl = fitlm(x, y);                                  % Fit Data
B = mdl.Coefficients.Estimate;                      % Coefficients
CI = coefCI(mdl);                                   % Coefficient Confidence Intervals
[Ypred,YCI] = predict(mdl, x);                      % Fitted Regression Line & Confidence Intervals

Xnew = linspace(min(x), max(x), 1000);
[~,YCI] = predict(mdl, Xnew'); 

%[B,BINT,R,RINT,STATS] = regress(y,[x,ones(size(x,1),1)])

%plot no correlation 0 line
plot([0 1],[0 0],'LineStyle',':','Color',[34, 139, 34]/255,'LineWidth',2)
plot(Xnew, YCI,'LineStyle', '--','Color',[0.5 0.5 0.5],'LineWidth',1)

plot(x, Ypred,'-k','LineWidth',2)
xlabel('Normalized position relative to A')
ylabel('B - A field position difference')
set(gca,'FontSize',18)
set(gca,'LineWidth',2)
hold off

%pearson R correlation
[r,p] = corrcoef(x,  y)


%% Scatter - Try this for global remapping neurons


%field difference B-A
field_pos_diff = global_far.B - global_far.A;

%merge global positions
global_far_merge = [global_far.A; global_far.B];

%global_near.B(track_fwd_idx(1)) - global_near.A(track_fwd_idx(1))

%account for edge effect on difference calculation - no account 
%shifting behind track edge
track_rev_idx = find(field_pos_diff >50);
%shifting forward of track edge
track_fwd_idx = find(field_pos_diff <-50);

%function of mean position 
mean_pos = mean([global_far.A; global_far.B],1); 

%account for edge effect
field_pos_diff(track_rev_idx) = field_pos_diff(track_rev_idx)-100;
field_pos_diff(track_fwd_idx) = field_pos_diff(track_fwd_idx)+100;


%convert to normalized space, assign as input variables
x = (global_far.A./100)';
y = (field_pos_diff./100)';

figure;
hold on
scatter(x,  y',14,'MarkerFaceColor',[139,0,139]/255,'MarkerEdgeColor',[139,0,139]/255)
axis square

%fit line
%x = global_near.A';                                        % Create Data
%y = field_pos_diff';               % Create Data
mdl = fitlm(x, y);                                  % Fit Data
B = mdl.Coefficients.Estimate;                      % Coefficients
CI = coefCI(mdl);                                   % Coefficient Confidence Intervals
[Ypred,YCI] = predict(mdl, x);                      % Fitted Regression Line & Confidence Intervals

Xnew = linspace(min(x), max(x), 1000);
[~,YCI] = predict(mdl, Xnew'); 

%[B,BINT,R,RINT,STATS] = regress(y,[x,ones(size(x,1),1)])

%plot no correlation 0 line
plot([0 1],[0 0],'LineStyle',':','Color',[34, 139, 34]/255,'LineWidth',2)
plot(Xnew, YCI,'LineStyle', '--','Color',[0.5 0.5 0.5],'LineWidth',1)

plot(x, Ypred,'-k','LineWidth',2)
%max position difference line for global
plot([0 1],[0.15 0.15],'r--')
plot([0 1],[-0.15 -0.15],'r--')
xlabel('Normalized position relative to A')
ylabel('B - A field position difference')
set(gca,'FontSize',18)
set(gca,'LineWidth',2)
hold off

%pearson R correlation
[r,p] = corrcoef(x,  y)

%% Combine all 3 classes

%merge global positions
all_far_merge = [[global_far.A; global_far.B], [global_near.A; global_near.B], common, rate];

%field difference B-A
field_pos_diff = all_far_merge(2,:) -  all_far_merge(1,:);

%global_near.B(track_fwd_idx(1)) - global_near.A(track_fwd_idx(1))

%account for edge effect on difference calculation - no account 
%shifting behind track edge
if 0 
track_rev_idx = find(field_pos_diff >50);
%shifting forward of track edge
track_fwd_idx = find(field_pos_diff <-50);

%account for edge effect
field_pos_diff(track_rev_idx) = field_pos_diff(track_rev_idx)-100;
field_pos_diff(track_fwd_idx) = field_pos_diff(track_fwd_idx)+100;
end

%convert to normalized space, assign as input variables
x = (all_far_merge(1,:)./100)';
y = (field_pos_diff./100)';



figure;
hold on
scatter(x,  y',9,'MarkerFaceColor',[139,0,139]/255,'MarkerEdgeColor',[139,0,139]/255)
axis square

%fit line
%x = global_near.A';                                        % Create Data
%y = field_pos_diff';               % Create Data
mdl = fitlm(x, y);                                  % Fit Data
B = mdl.Coefficients.Estimate;                      % Coefficients
CI = coefCI(mdl);                                   % Coefficient Confidence Intervals
[Ypred,YCI] = predict(mdl, x);                      % Fitted Regression Line & Confidence Intervals

Xnew = linspace(min(x), max(x), 1000);
[~,YCI] = predict(mdl, Xnew'); 

%[B,BINT,R,RINT,STATS] = regress(y,[x,ones(size(x,1),1)])

%plot no correlation 0 line
plot([0 1],[0 0],'LineStyle',':','Color',[34, 139, 34]/255,'LineWidth',2)
plot(Xnew, YCI,'LineStyle', '--','Color',[0.5 0.5 0.5],'LineWidth',1)

plot(x, Ypred,'-k','LineWidth',2)
%max position difference line for global
%plot([0 1],[0.15 0.15],'r--')
%plot([0 1],[-0.15 -0.15],'r--')
xlabel('Normalized position relative to A')
ylabel('B - A field position difference')
set(gca,'FontSize',18)
set(gca,'LineWidth',2)
hold off

%pearson R correlation
[r,p] = corrcoef(x,  y)




%% plot histograms
figure;
hold on
histogram(global_near.A(less_than_unity),[0:20:100])
histogram(global_near.A(greater_than_unity),[0:20:100])



[h,p] =kstest2(global_near.A(less_than_unity),global_near.B(less_than_unity))

scatter(global_far_comb(1,:),global_far_comb(2,:))

figure
subplot(5,1,1)
hold on
ylim([0 ylim_def])
histogram(mean(common,1),[0:bin_space:100],'Normalization','probability')
subplot(5,1,2)
hold on
ylim([0 ylim_def])
histogram(global_far_comb(1,:),[0:bin_space:100],'Normalization','probability')
subplot(5,1,3)
hold on
ylim([0 ylim_def])
histogram(global_far_comb(2,:),[0:bin_space:100],'Normalization','probability')
subplot(5,1,4)
hold on
ylim([0 ylim_def])
histogram(partial_common,[0:bin_space:100],'Normalization','probability')
subplot(5,1,5)
hold on
ylim([0 ylim_def])
histogram(partial.com_parA(:,2),[0:bin_space:100],'Normalization','probability')
histogram(partial.com_parB(:,2),[0:bin_space:100],'Normalization','probability')

%% Ecdf of values above
figure;
hold on
ecdf(mean(common,1))
ecdf(mean(rate,1))
ecdf(global_far_comb(1,:))
ecdf(global_far_comb(2,:))
% ecdf(partial_common)
% ecdf(partial.com_parA(:,2))
% ecdf(partial.com_parB(:,2))

figure;
hold on
histogram(binCenter_data{1}.bin_center.global_far(1,:),0:10:100)

figure;
hold on
histogram(binCenter_data{1}.bin_center.global_far(2,:),0:10:100)

figure
hold on
title('Far global remappers - center of A vs center of B')
xlabel('A field center')
ylabel('B field center')
scatter(binCenter_data{1}.bin_center.global_far_sortA(1,:),binCenter_data{1}.bin_center.global_far_sortA(2,:))

figure
hold on
title('Near global remappers - center of A vs center of B')
xlabel('A field center')
ylabel('B field center')
scatter(binCenter_data{1}.bin_center.global_near(1,:),binCenter_data{1}.bin_center.global_near(2,:))
plot([0 100],[0 100],'k--')

figure
hold on
title('Common - center of A vs center of B')
xlabel('A field center')
ylabel('B field center')
scatter(binCenter_data{1}.bin_center.common(1,:),binCenter_data{1}.bin_center.common(2,:))
plot([0 100],[0 100],'k--')

%% For partial, scatter plot of centroids - sort my mean of common place field center
%play with later
%{
[~,Isort_com] = sort(mean(bin_center.partial_com),'descend');
%sort columns by common mean center
bin_center.partial_com_sort = bin_center.partial_com(:,Isort_com);
bin_center.partial_far_sort =bin_center.partial_far(:,Isort_com);

figure
%plot common
hold on
for rr=1:size(bin_center.partial_com,2)
    scatter(bin_center.partial_com_sort(1,rr),rr,'g')
    scatter(bin_center.partial_com_sort(2,rr),rr,'g')
end
%plot partial
for rr=1:size(bin_center.partial_com,2)
if ~isnan(bin_center.partial_far_sort(1,rr))
    scatter(bin_center.partial_far_sort(1,rr),rr,'b')
else
    scatter(bin_center.partial_far_sort(2,rr),rr,'r')
end
end

%sort global far by B trials
[~,Isort_glo_far_B] = sort(bin_center.global_far(2,:),'descend');
[~,Isort_glo_far_A] = sort(bin_center.global_far(1,:),'descend');
%sort columns by common mean center
bin_center.global_far_sortB = bin_center.global_far(:,Isort_glo_far_B);
bin_center.global_far_sortA = bin_center.global_far(:,Isort_glo_far_A);
%bin_center.partial_far_sort =bin_center.partial_far(:,Isort_com);

%Plot global far
figure;
subplot(1,2,1)
hold on;
for rr=1:size(bin_center.global_far,2)
    scatter(bin_center.global_far_sortA(1,rr),rr,'b')
    scatter(bin_center.global_far_sortA(2,rr),rr,'r')
end
subplot(1,2,2)
hold on;
for rr=1:size(bin_center.global_far,2)
    scatter(bin_center.global_far_sortB(1,rr),rr,'b')
    scatter(bin_center.global_far_sortB(2,rr),rr,'r')
end

figure;
hold on
histogram(bin_center.global_far(1,:),0:10:100)

figure;
hold on
histogram(bin_center.global_far(2,:),0:10:100)

figure
hold on
title('Far global remappers - center of A vs center of B')
xlabel('A field center')
ylabel('B field center')
scatter(bin_center.global_far_sortA(1,:),bin_center.global_far_sortA(2,:))

figure
hold on
title('Near global remappers - center of A vs center of B')
xlabel('A field center')
ylabel('B field center')
scatter(bin_center.global_near(1,:),bin_center.global_near(2,:))
plot([0 100],[0 100],'k--')

figure
hold on
title('Common - center of A vs center of B')
xlabel('A field center')
ylabel('B field center')
scatter(bin_center.common(1,:),bin_center.common(2,:))
plot([0 100],[0 100],'k--')

%}
%}
end

