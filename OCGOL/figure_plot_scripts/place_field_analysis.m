function [outputArg1,outputArg2] = place_field_analysis(path_dir)

%% Load relevant place field data

%read relevant data
for ee=1:size(path_dir,2)
    load_data_path_pf{ee} = fullfile(path_dir{ee},'cumul_analysis','placeField_dist.mat');
    placeField_data{ee} = load(string(load_data_path_pf{ee}));
end

%combine field counts for Asel and Bsel into 1 matrix
%create matrices with centroid counts for each subclass of tuned neuron
for ee=1:size(path_dir,2)
    %counts of each field
    %A selective
    pf_count.Asel(ee,:) = placeField_data{ee}.placeField_dist.task_sel.A.field_count_total;
    %B selective
    pf_count.Bsel(ee,:) = placeField_data{ee}.placeField_dist.task_sel.B.field_count_total;
    
    %A and B selective using either criterion
    %A trials
    pf_count.AB.A(ee,:) = placeField_data{ee}.placeField_dist.other_classes.si_ts.AB.A.field_count_total;
    %B trials
    pf_count.AB.B(ee,:) = placeField_data{ee}.placeField_dist.other_classes.si_ts.AB.B.field_count_total;

%width of each field
    pf_width.Asel{ee} = placeField_data{ee}.placeField_dist.task_sel.A.width_cm;
    %B selective
    pf_width.Bsel{ee} = placeField_data{ee}.placeField_dist.task_sel.B.width_cm;
    
    %A and B selective using either criterion
    %A trials
    pf_width.AB.A{ee} = placeField_data{ee}.placeField_dist.other_classes.si_ts.AB.A.width_cm;
    %B trials
    pf_width.AB.B{ee} = placeField_data{ee}.placeField_dist.other_classes.si_ts.AB.B.width_cm;

end

%normalize as fraction of neurons for each animal/exp for A-sel/B-sel
pf_count.Asel_norm = pf_count.Asel./sum(pf_count.Asel,2);
pf_count.Bsel_norm = pf_count.Bsel./sum(pf_count.Bsel,2);

pf_count.AB.A_norm = pf_count.AB.A./sum(pf_count.AB.A,2);
pf_count.AB.B_norm = pf_count.AB.B./sum(pf_count.AB.B,2);

%get means for each subclass (norm data)
mean_pf_count.Asel = nanmean(pf_count.Asel_norm,1);
mean_pf_count.Bsel = nanmean(pf_count.Bsel_norm,1);

mean_pf_count.AB.A = nanmean(pf_count.AB.A_norm,1);
mean_pf_count.AB.B = nanmean(pf_count.AB.B_norm,1);

%get std and sem for each group
std_pf_count.Asel = nanstd(pf_count.Asel_norm,0,1);
std_pf_count.Bsel = nanstd(pf_count.Bsel_norm,0,1);

std_pf_count.AB.A = nanstd(pf_count.AB.A_norm,0,1);
std_pf_count.AB.B = nanstd(pf_count.AB.B_norm,0,1);

%sem
sem_pf_count.Asel = std_pf_count.Asel./sqrt(size(pf_count.Asel_norm,1));
sem_pf_count.Bsel = std_pf_count.Bsel./sqrt(size(pf_count.Bsel_norm,1));

sem_pf_count.AB.A = std_pf_count.AB.A./sqrt(size(pf_count.AB.A_norm,1));
sem_pf_count.AB.B = std_pf_count.AB.B./sqrt(size(pf_count.AB.B_norm,1));

%grouped sem
grouped_norm_sem = [sem_pf_count.Asel',sem_pf_count.Bsel',sem_pf_count.AB.A',sem_pf_count.AB.B'];

%combined means from norm counts
grouped_norm_mean = [mean_pf_count.Asel', mean_pf_count.Bsel',mean_pf_count.AB.A',mean_pf_count.AB.B'];

%% Statistics on place field numbers (A vs. B; A vs A-sel; B vs B-sel)

%single
single_pf_frac = [pf_count.Asel_norm(:,1),pf_count.Bsel_norm(:,1),pf_count.AB.A_norm(:,1),pf_count.AB.B_norm(:,1)];
double_pf_frac = [pf_count.Asel_norm(:,2),pf_count.Bsel_norm(:,2),pf_count.AB.A_norm(:,2),pf_count.AB.B_norm(:,2)];
triple_pf_frac = [pf_count.Asel_norm(:,3),pf_count.Bsel_norm(:,3),pf_count.AB.A_norm(:,3),pf_count.AB.B_norm(:,3)];

p_single(1) = signrank(single_pf_frac(:,1),single_pf_frac(:,2))
p_single(2) = signrank(single_pf_frac(:,1),single_pf_frac(:,3))
p_single(3) = signrank(single_pf_frac(:,2),single_pf_frac(:,4))

p_double(1) = signrank(double_pf_frac(:,1),double_pf_frac(:,2))
p_double(2) = signrank(double_pf_frac(:,1),double_pf_frac(:,3))
p_double(3) = signrank(double_pf_frac(:,2),double_pf_frac(:,4))

p_triple(1) = signrank(triple_pf_frac(:,1),triple_pf_frac(:,2))
p_triple(2) = signrank(triple_pf_frac(:,1),triple_pf_frac(:,3))
p_triple(3) = signrank(triple_pf_frac(:,2),triple_pf_frac(:,4))

%% Define colors
%blue, red, light blue, light red      
color_mat = [65,105,225; 220,20,60; 135,206,250; 240,128,128]./255;

%% Width mean and sem 

%get mean place field width for each animal
for aa=1:size(path_dir,2)
    mean_width.Asel(aa) = nanmean(pf_width.Asel{aa});
    mean_width.Bsel(aa) = nanmean(pf_width.Bsel{aa});
    mean_width.AB.A(aa) = nanmean(pf_width.AB.A{aa});
    mean_width.AB.B(aa) = nanmean(pf_width.AB.B{aa});
end

mean(mean_width.Asel)
mean(mean_width.Bsel)
mean(mean_width.AB.A)
mean(mean_width.AB.B)

%every 4 cm from 18-100
edges = [10:2:70];
edge_center = edges(1:end-1)+1;
for aa=1:size(path_dir,2)
    [N_Asel(aa,:),~] = histcounts(pf_width.Asel{aa},edges,'Normalization','probability');
    [N_Bsel(aa,:),~] = histcounts(pf_width.Bsel{aa},edges,'Normalization','probability');
    [N_AB_A(aa,:),~] = histcounts(pf_width.AB.A{aa},edges,'Normalization','probability');
    [N_AB_B(aa,:),~] = histcounts(pf_width.AB.B{aa},edges,'Normalization','probability');
end

%get cumulative probability
cum_prob_A_sel = cumsum(N_Asel,2);
cum_prob_B_sel = cumsum(N_Bsel,2);
cum_prob_AB_A = cumsum(N_AB_A,2);
cum_prob_AB_B = cumsum(N_AB_B,2);

%get sem at each fractional point
sem_cum_prob_A = nanstd(cum_prob_A_sel,0,1)./sqrt(size(path_dir,2));
sem_cum_prob_B = nanstd(cum_prob_B_sel,0,1)./sqrt(size(path_dir,2));
sem_cum_prob_AB_A = nanstd(cum_prob_AB_A,0,1)./sqrt(size(path_dir,2));
sem_cum_prob_AB_B = nanstd(cum_prob_AB_B,0,1)./sqrt(size(path_dir,2));

%linewidth
lw = 2;
figure
hold on
xlabel('Place field width [cm]')
ylabel('Fraction of place fields')
yticks([0 0.2 0.4 0.6 0.8 1])
% plot(edge_center,mean(cum_prob_A_sel,1),'LineWidth',lw,'Color',color_mat(1,:))
% plot(edge_center,mean(cum_prob_B_sel,1),'LineWidth',lw,'Color',color_mat(2,:))
% plot(edge_center,mean(cum_prob_AB_A,1),'LineWidth',lw,'Color',color_mat(3,:))
% plot(edge_center,mean(cum_prob_AB_B,1),'LineWidth',lw,'Color',color_mat(4,:))

e3 = errorbar(edge_center,mean(cum_prob_AB_A,1),sem_cum_prob_AB_A,'LineWidth',lw,'Color',color_mat(3,:));
e4 = errorbar(edge_center,mean(cum_prob_AB_B,1),sem_cum_prob_AB_B,'LineWidth',lw,'Color',color_mat(4,:));
e1 = errorbar(edge_center,mean(cum_prob_A_sel,1),sem_cum_prob_A,'LineWidth',lw,'Color',color_mat(1,:));
e2 = errorbar(edge_center,mean(cum_prob_B_sel,1),sem_cum_prob_B,'LineWidth',lw,'Color',color_mat(2,:));

legend([e1 e2 e3 e4],{'A selective','B selective','A&B - A', 'A&B - B'},'location','southeast')

set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)

%% Pool place field width data together for export for ks test analysis

%move this data to prism
pf_width_pool.Asel = cell2mat(pf_width.Asel)';
pf_width_pool.Bsel = cell2mat(pf_width.Bsel)';
pf_width_pool.AB.A = cell2mat(pf_width.AB.A)';
pf_width_pool.AB.B = cell2mat(pf_width.AB.B)';

[h, p] = kstest2(pf_width_pool.Asel,pf_width_pool.Bsel)

figure
hold on
xlim([10 70])
cdfplot(pf_width_pool.Asel)
cdfplot(pf_width_pool.Bsel)
cdfplot(pf_width_pool.AB.A)
cdfplot(pf_width_pool.AB.B)



%% Place field analysis plotting

%plot bar
figure('Position',[2314 282 427 420]);
hold on;
axis square
title('Place fields per neuron - selective and A&B (both crit)');
%bar the mean for each group
b = bar(1:3,grouped_norm_mean,'FaceColor', 'flat');
pause(0.1)
xlim([0.5 3.5])
ylim([0 1])
%plot the sem for each mean for each group
for ib = 1:numel(b)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    if ib ==1
        xData(1,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==2
        xData(2,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==3
        xData(3,:) = b(ib).XData + b(ib).XOffset;
    elseif ib ==4
        xData(4,:) = b(ib).XData + b(ib).XOffset;
    end
    errorbar(xData(ib,:),grouped_norm_mean(:,ib)',grouped_norm_sem(:,ib),'k.','LineWidth',1)
end

%set A group bars to blue
b(1).CData(1:3,:) =  repmat(color_mat(1,:),3,1);
%set B group bars to red
b(2).CData(1:3,:) =  repmat(color_mat(2,:),3,1);
%set B group bars to red
b(3).CData(1:3,:) =  repmat(color_mat(3,:),3,1);
%set B group bars to red
b(4).CData(1:3,:) =  repmat(color_mat(4,:),3,1);

xticks([1 2 3]);
xticklabels({'1','2','3+'});
ylabel('Fraction of neurons');
legend('A selective','B selective','A&B - A', 'A&B - B')

set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)

%sigstar({xData(1,1),xData(3,1)},[])

%sum A and B
% grouped_pf_counts = [sum(pf_count_mat_A,1)',sum(pf_count_mat_B,1)'];
% %normalized for each group
% grouped_pf_counts_norm = [(sum(pf_count_mat_A,1)./sum(sum(pf_count_mat_A,1)))',...
%             (sum(pf_count_mat_B,1)./sum(sum(pf_count_mat_B,1)))'];
%         


end

