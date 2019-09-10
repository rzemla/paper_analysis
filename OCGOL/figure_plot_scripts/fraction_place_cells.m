function [outputArg1,outputArg2] = fraction_place_cells(path_dir)

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
    frac_mat_si(ee,:) = fractional_data{ee}.tuned_fractions.tuned_count_filt{1};
    %TS
    frac_mat_ts(ee,:) = fractional_data{ee}.tuned_fractions.tuned_count_filt{2};
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

end

