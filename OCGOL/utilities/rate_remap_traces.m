function [outputArg1,outputArg2] = rate_remap_traces(path_dir)

%% Load in data
for aa=1:size(path_dir,2)
    load_data_path{aa} = fullfile(path_dir{aa},'cumul_analysis','remap_traces.mat');
    remap_traces_idx{aa} = load(string(load_data_path{aa}));
end

%% Assemble into more compact representation

%combine into 1 big matrix
mean_dff_matrix.A = zeros(64,300);
mean_dff_matrix.B = zeros(64,300);

sp_idx =1;
for aa=1:11
    for rr=1:size(remap_traces_idx{aa}.remap_idx_traces,2)
        
        %find max number of frames for each mean trace
        max_frame = [length(mean(remap_traces_idx{aa}.remap_idx_traces{rr}.A,1)),length(mean(remap_traces_idx{aa}.remap_idx_traces{rr}.B,1))];
        
        mean_dff_matrix.A(sp_idx,1:max_frame(1)) = mean(remap_traces_idx{aa}.remap_idx_traces{rr}.A,1);
        mean_dff_matrix.B(sp_idx,1:max_frame(2)) = mean(remap_traces_idx{aa}.remap_idx_traces{rr}.B,1);
        
        sp_idx = sp_idx + 1;
    end
end

remap_traces_idx{1, 1}.remap_idx_traces{1, 1}.A  
remap_traces_idx{1, 1}.remap_idx_traces{1, 1}.A 

%% Sort and plot
max_A_value = max(mean_dff_matrix.A,[],2);
max_B_value = max(mean_dff_matrix.B,[],2);

[~,I_sort] = sort(max_A_value - max_B_value,'descend');

figure
subplot(1,2,1)
imagesc(mean_dff_matrix.A(I_sort,:))
hold on
colormap('jet')
caxis([0 3])

subplot(1,2,2)
imagesc(mean_dff_matrix.B(I_sort,:))
hold on
colormap('jet')
caxis([0 3])


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
    
    imagesc(mean_dff_matrix.A(I_sort,1:200))
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
    imagesc(mean_dff_matrix.B(I_sort,1:200))
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

