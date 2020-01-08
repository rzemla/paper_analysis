function [outputArg1,outputArg2] = rate_remap_traces(path_dir)

%% Load in data (includes common neurons as well)
for aa=1:size(path_dir,2)
    load_data_path{aa} = fullfile(path_dir{aa},'cumul_analysis','remap_traces.mat');
    remap_traces_idx{aa} = load(string(load_data_path{aa}));
    com_traces_idx{aa} = remap_traces_idx{aa}.com_idx_traces;
end

%% Load in AUC data
for aa=1:size(path_dir,2)
    load_data_path{aa} = fullfile(path_dir{aa},'cumul_analysis','AUC_remapping.mat');
    AUC_data{aa} = load(string(load_data_path{aa}));
    
    AUC.common{aa} = AUC_data{aa}.AUC_remappers.common;
    AUC.remap{aa} = AUC_data{aa}.AUC_remappers.remap;
end


%% Assemble into more compact representation (activity remapping neurons)

%combine into 1 big matrix
mean_dff_matrix.remap.A = zeros(64,300);
mean_dff_matrix.remap.B = zeros(64,300);

sp_idx =1;
for aa=1:11
    for rr=1:size(remap_traces_idx{aa}.remap_idx_traces,2)
        
        %find max number of frames for each mean trace
        max_frame(sp_idx,:) = [length(mean(remap_traces_idx{aa}.remap_idx_traces{rr}.A,1)),length(mean(remap_traces_idx{aa}.remap_idx_traces{rr}.B,1))];
        
        mean_dff_matrix.remap.A(sp_idx,1:max_frame(sp_idx,1)) = mean(remap_traces_idx{aa}.remap_idx_traces{rr}.A,1);
        mean_dff_matrix.remap.B(sp_idx,1:max_frame(sp_idx,2)) = mean(remap_traces_idx{aa}.remap_idx_traces{rr}.B,1);
        
        sp_idx = sp_idx + 1;
    end
end

% remap_traces_idx{1, 1}.remap_idx_traces{1, 1}.A  
% remap_traces_idx{1, 1}.remap_idx_traces{1, 1}.A 

%% Assemble into more compact representation (common  neurons)

%extract values for the matrix size here
sp_idx =1;
for aa=1:11
    for rr=1:size(remap_traces_idx{aa}.com_idx_traces,2)
        
        %find max number of frames for each mean trace
        max_frame_com(sp_idx,:) = [length(mean(remap_traces_idx{aa}.com_idx_traces{rr}.A,1)),length(mean(remap_traces_idx{aa}.com_idx_traces{rr}.B,1))];
        sp_idx = sp_idx + 1;
    end
end
        
%combine into 1 big matrix 
mean_dff_matrix.com.A = zeros(size(max_frame_com,1),max(max(max_frame_com)));
mean_dff_matrix.com.B = zeros(size(max_frame_com,1),max(max(max_frame_com)));

sp_idx =1;
for aa=1:11
    for rr=1:size(remap_traces_idx{aa}.com_idx_traces,2)
        
        %find max number of frames for each mean trace
        max_frame_com(sp_idx,:) = [length(mean(remap_traces_idx{aa}.com_idx_traces{rr}.A,1)),length(mean(remap_traces_idx{aa}.com_idx_traces{rr}.B,1))];
        
        mean_dff_matrix.com.A(sp_idx,1:max_frame_com(sp_idx,1)) = mean(remap_traces_idx{aa}.com_idx_traces{rr}.A,1);
        mean_dff_matrix.com.B(sp_idx,1:max_frame_com(sp_idx,2)) = mean(remap_traces_idx{aa}.com_idx_traces{rr}.B,1);
        
        sp_idx = sp_idx + 1;
    end
end



%% Sort and plot
max_A_value = max(mean_dff_matrix.remap.A,[],2);
max_B_value = max(mean_dff_matrix.remap.B,[],2);

%get duration of each event - later if necessary - purely cosmetic
%mean_dff_matrix.A)

[~,I_sort] = sort(max_A_value - max_B_value,'descend');

figure
subplot(1,2,1)
imagesc(mean_dff_matrix.remap.A(I_sort,:))
hold on
colormap('jet')
caxis([0 3])

subplot(1,2,2)
imagesc(mean_dff_matrix.remap.B(I_sort,:))
hold on
colormap('jet')
caxis([0 3])

max(max(mean_dff_matrix.remap.B))

%% Calculate index based on sum of AUC values (flawed b/c not taking into account number of events)

%calculate AUC sum for A or B trials
for aa=1:size(AUC_data,2)
    %for common neurons
    %for A laps
    AUC_sum_com{aa}(:,1) = cellfun(@sum,AUC.common{aa}(:,1));
    %for B laps
    AUC_sum_com{aa}(:,2) = cellfun(@sum,AUC.common{aa}(:,2));
    
    %for remapping neurons
    %for A laps
    AUC_sum_remap{aa}(:,1) = cellfun(@sum,AUC.remap{aa}(:,1));
    %for B laps
    AUC_sum_remap{aa}(:,2) = cellfun(@sum,AUC.remap{aa}(:,2));
    
end

%combine common sums into 1 matrix
AUC_sum_com_combined = cell2mat(AUC_sum_com');

AUC_sum_remap_combined = cell2mat(AUC_sum_remap');


%calculate index
com_idx_AUC  = abs(AUC_sum_com_combined(:,1) - AUC_sum_com_combined(:,2))./(AUC_sum_com_combined(:,1) + AUC_sum_com_combined(:,2));
remap_idx_AUC  = abs(AUC_sum_remap_combined (:,1) - AUC_sum_remap_combined (:,2))./(AUC_sum_remap_combined (:,1) + AUC_sum_remap_combined (:,2));

figure
hold on
histogram(com_idx_AUC,10,'Normalization','probability')
histogram(remap_idx_AUC,10,'Normalization','probability')

%% Calculate index for common and remapping neurons and plot (this is the index you want)

max_pk_A_value_remap = max(mean_dff_matrix.remap.A,[],2);
max_pk_B_value_remap = max(mean_dff_matrix.remap.B,[],2);

max_pk_A_value_com = max(mean_dff_matrix.com.A,[],2);
max_pk_B_value_com = max(mean_dff_matrix.com.B,[],2);

%index for max peak value
remap_idx_values = abs(max_pk_A_value_remap - max_pk_B_value_remap)./(max_pk_A_value_remap + max_pk_B_value_remap);

common_idx_values = abs(max_pk_A_value_com - max_pk_B_value_com)./(max_pk_A_value_com + max_pk_B_value_com);

%plot both distributions on 1 histogram
figure
hold on
histogram(common_idx_values,20,'Normalization','probability','DisplayStyle','stairs')
histogram(remap_idx_values,20,'Normalization','probability','DisplayStyle','stairs')


[h,p] = kstest2(remap_idx_values,common_idx_values)

%% Plot the cdf of the peak indices for common vs. remapping neurons
figure
hold on
e1 = cdfplot(common_idx_values);
e2 = cdfplot(remap_idx_values);

e1.LineWidth = 1.5;
e1.Color = [0 0 0];

e2.LineWidth = 1.5;
e2.Color = [255,140,0]./255;

grid off
title(' ')
xlim([ 0 1])
xlabel('Peak index')
ylabel('Cumulative fraction')
legend([e1 e2],'Common','Activity remapping','location','southeast')
set(gca,'FontSize',14)
set(gca,'LineWidth',1)

%% Calculate area under mean of the curve and calculate index
% area_mean_pk_A_value_remap = trapz(mean_dff_matrix.remap.A');
% area_mean_pk_B_value_remap = trapz(mean_dff_matrix.remap.B');
% 
% area_mean_pk_A_value_com = trapz(mean_dff_matrix.com.A');
% area_mean_pk_B_value_com = trapz(mean_dff_matrix.com.B');
% 
% remap_idx_area_values = abs(area_mean_pk_A_value_remap - area_mean_pk_B_value_remap)./...
%         (area_mean_pk_A_value_remap + area_mean_pk_B_value_remap);
% 
% com_idx_area_values = abs(area_mean_pk_A_value_com - area_mean_pk_B_value_com)./...
%         (area_mean_pk_A_value_com + area_mean_pk_B_value_com);    
    

%% Subtract the maps from one another and get max value and calculate index

diff_mean_traces_com = mean_dff_matrix.com.A - mean_dff_matrix.com.B;

diff_mean_traces_remap = mean_dff_matrix.remap.A - mean_dff_matrix.remap.B;

area_com = abs(max(diff_mean_traces_com,[],2) - min(diff_mean_traces_com,[],2))./(max(diff_mean_traces_com,[],2) + min(diff_mean_traces_com,[],2))


area_remap = abs(max(diff_mean_traces_remap,[],2) - min(diff_mean_traces_remap,[],2))./(max(diff_mean_traces_remap,[],2) + min(diff_mean_traces_remap,[],2))

%plot both distributions on 1 histogram
figure
hold on
histogram(common_idx_values,20,'Normalization','probability','DisplayStyle','stairs')
histogram(remap_idx_values,20,'Normalization','probability','DisplayStyle','stairs')


%% Plot columns in same order as below except use same jet colormap

ROI_type_order = [1 2; 3 4; 5 6];
subplot_order = [1 3 5; 2 4 6]';
f = figure('Position', [2475 75 372 898]); % event based STC;
for cc =1
    
    %subplot(3,2,subplot_order(cc))
    subaxis(3,2,subplot_order(cc,1),'SpacingHorizontal', 0.02,...
        'SpacingVertical',0.01,...
        'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
    
    imagesc(mean_dff_matrix.remap.A(I_sort,1:200))
    hold on
    ax1 = gca;
    
    colormap(ax1,jet)
    caxis(ax1,[0 3])
    
    %ax1.YTick = [1,100];
    ax1.XTick = [1,30,90,150];
    ax1.XTickLabel = {'0','1','3','5'};
    ax1.YAxis.TickLength = [0 0];
    
    
end

for cc =1
    
    subaxis(3,2,subplot_order(cc,2),'SpacingHorizontal', 0.015,...
        'SpacingVertical',0.01,...
        'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
    imagesc(mean_dff_matrix.remap.B(I_sort,1:200))
    %title('5A5B')
    hold on
    ax2 = gca;
    %cbar= colorbar;
    %cbar.Label.String = 'Normalized activity';
    %cbar.Ticks = [0 1];
    colormap(ax2,jet)
    caxis(ax2,[0 3])

        ax2.YTick = [];%[1,100];
        ax2.XTickLabel = {'0','1'};
        ax2.YAxis.TickLength = [0 0];

            %ax1.YTick = [1,100];
    ax2.XTick = [1,30,90,150];
    ax2.XTickLabel = {'0','1','3','5'};
    ax2.YAxis.TickLength = [0 0];

end

%plot associated colorbar
f = figure('Position', [2475 75 372 898]); % event based STC;
hold on
subaxis(3,2,subplot_order(cc,1),'SpacingHorizontal', 0.02,...
    'SpacingVertical',0.01,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
colormap(jet)
caxis([0 3])
colorbar


%% first column plot

cmap_blue=cbrewer('seq', 'Blues', 32);
cmap_red=cbrewer('seq', 'Reds', 32);

%set bottom value to red to bottom value of blue
cmap_red(1,:) = [1 1 1];
cmap_blue(1,:) = [1 1 1];

ROI_type_order = [1 2; 3 4; 5 6];
subplot_order = [1 3 5; 2 4 6]';
f = figure('Position', [2475 75 372 898]); % event based STC;
for cc =1
    
    %subplot(3,2,subplot_order(cc))
    subaxis(3,2,subplot_order(cc,1),'SpacingHorizontal', 0.02,...
        'SpacingVertical',0.01,...
        'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
    
    imagesc(mean_dff_matrix.remap.A(I_sort,1:200))
    hold on
    ax1 = gca;
    
    colormap(ax1,cmap_blue)
    caxis(ax1,[0 3])
    
    %ax1.YTick = [1,100];
    ax1.XTick = [1,30,90,150];
    ax1.XTickLabel = {'0','1','3','5'};
    ax1.YAxis.TickLength = [0 0];
    
end

for cc =1
    
    subaxis(3,2,subplot_order(cc,2),'SpacingHorizontal', 0.015,...
        'SpacingVertical',0.01,...
        'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
    imagesc(mean_dff_matrix.remap.B(I_sort,1:200))
    %title('5A5B')
    hold on
    ax2 = gca;
    %cbar= colorbar;
    %cbar.Label.String = 'Normalized activity';
    %cbar.Ticks = [0 1];
    colormap(ax2,cmap_red)
    caxis(ax2,[0 3])

        ax2.YTick = [];%[1,100];
        ax2.XTickLabel = {'0','1'};
        ax2.YAxis.TickLength = [0 0];

            %ax1.YTick = [1,100];
    ax2.XTick = [1,30,90,150];
    ax2.XTickLabel = {'0','1','3','5'};
    ax2.YAxis.TickLength = [0 0];

end


%% Plot mean/sem each remapping neuron 
figure
sp_idx=1;
for aa=1:11
    for rr=1:size(remap_traces_idx{aa}.remap_idx_traces,2)
        subplot(8,8,sp_idx)
        hold on
        %find max and add 0.5
        max_dFF = max([max(mean(remap_traces_idx{aa}.remap_idx_traces{rr}.A,1)), max(mean(remap_traces_idx{aa}.remap_idx_traces{rr}.B,1))]);
        max_frame = max([max(mean(remap_traces_idx{aa}.remap_idx_traces{rr}.A,2)), max(mean(remap_traces_idx{aa}.remap_idx_traces{rr}.B,2))]);
        ylim([0 round(max_dFF,1)+0.5])
        xlim([0 200])
        plot_dFF_mean_trace(remap_traces_idx{aa}.remap_idx_traces{rr}.A, remap_traces_idx{aa}.remap_idx_traces{rr}.B)
        sp_idx = sp_idx +1;
        hold off
        %pause
    end
end


end

