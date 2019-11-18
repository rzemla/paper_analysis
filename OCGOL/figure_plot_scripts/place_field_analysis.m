function [outputArg1,outputArg2] = place_field_analysis(path_dir)

%% 

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


end

