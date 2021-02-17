function [outputArg1,outputArg2] = fig_sup_activity_remap_master(remap_sup_plot_data)

%% Unload data

common_idx_values = remap_sup_plot_data.common_idx_values;
remap_idx_values = remap_sup_plot_data.remap_idx_values;

frames_sec = remap_sup_plot_data.frames_sec;
A_sel_idx = remap_sup_plot_data.A_sel_idx;
B_sel_idx = remap_sup_plot_data.B_sel_idx;
remap_traces_idx = remap_sup_plot_data.remap_traces_idx;
activity_remap_pk_idx_split = remap_sup_plot_data.activity_remap_pk_idx_split;
activity_remap_mean_AUC_min = remap_sup_plot_data.activity_remap_mean_AUC_min;
activity_remap_sem_AUC_min = remap_sup_plot_data.activity_remap_sem_AUC_min;

%% Plot master figure
%paper colors
paper_cmap = return_paper_colormap;

fig = figure;
fig.Units = 'centimeters';
fig.Position(1) = 7;
fig.Position(2) = 0;
fig.Position(3) = 54;
fig.Position(4) = 18;

%master layout
gridSize = [2,6];
t1 = tiledlayout(fig,gridSize(1),gridSize(2),'TileSpacing','compact','Padding','compact','Units','centimeters');

%subplot for activity remapping ecdf
s1 = tiledlayout(t1,3,2,'TileSpacing','compact','Padding','compact','Units','centimeters');
s1.Layout.Tile = 1;
s1.Layout.TileSpan = [2,2];
%s1.PositionConstraint = 'innerposition';

% Plot the cdf of the peak indices for common vs. remapping neurons
nexttile(s1,1,[2,2])
hold on
axis square
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

%blank graphs for scaling
%nexttile(s1,5,[1,2])

% Generate the plot for each neuron and AUC/min values
%top row : 1-3 A > B
%low to high; animal number first, then neuron idx (linear:
%5,6 ; 11,10; 11,11
%bottom row: 5-7 B > A
%4,2; 10,8; 7,1

%get frame to second equivalenet 1,2,3,4,5 s in frames
frames_sec = round([2,4,6]./0.0334);

%animal number first; then ROI number (linear, not absolute index)
A_sel_idx = [5,6;
            11,10;
            11,11];
        
B_sel_idx = [4,2;
            10,8;
            7,1];
A_tile_order = [3,4,5];
B_tile_order = [9,10,11];

        %A > B
        for ii=1:3
            %for top row
            nexttile(t1,A_tile_order(ii))
            hold on
            axis square
            ylim([0 2.2])
            xlim([0 200])
            yticks([0 1 2])
            xticks(frames_sec)
            xticklabels({'2','4','6'})
            aa=A_sel_idx(ii,1); rr=A_sel_idx(ii,2);
            title(num2str(round(activity_remap_pk_idx_split{aa}(rr),2)))
            plot_dFF_mean_trace(remap_traces_idx{aa}.remap_idx_traces{rr}.A, remap_traces_idx{aa}.remap_idx_traces{rr}.B)
            legend({'A','B'})
            if ii ==1
               ylabel('dF/F') 
            end
            set(gca,'FontSize',14)
            set(gca,'LineWidth',1.5)
        end
%B>A
        for ii=1:3
            nexttile(t1,B_tile_order(ii))
            hold on
            axis square
            ylim([0 2.2])
            xlim([0 200])
            yticks([0 1 2])
            xticks(frames_sec)
            xticklabels({'2','4','6'})
            aa=B_sel_idx(ii,1); rr=B_sel_idx(ii,2);
            title(num2str(round(activity_remap_pk_idx_split{aa}(rr),2)))
            plot_dFF_mean_trace(remap_traces_idx{aa}.remap_idx_traces{rr}.A, remap_traces_idx{aa}.remap_idx_traces{rr}.B)
            legend({'A','B'})
            xlabel('Time [s]')
            if ii ==1
               ylabel('dF/F') 
            end
            set(gca,'FontSize',14)
            set(gca,'LineWidth',1.5)
            
        end

%AUC/min comparison
    nexttile(t1,6)
    hold on
    axis square
    title('RUN')
    xlim([0 3])
    xticks([1 2])
    xticklabels({'A','B'})
    ylim([0 7])
    yticks([0 2 4 6])
    ylabel('AUC/min')
    b =bar([1 2],activity_remap_mean_AUC_min','FaceColor','flat');
    %change color of each bar
    b.CData(1,:) = paper_cmap(1,:);
    b.CData(2,:) = paper_cmap(2,:);
    
    errorbar([1 2],activity_remap_mean_AUC_min, activity_remap_sem_AUC_min,'LineStyle','none','LineWidth',1.5,'Color','k')
    %sigstar([1 2])
    
    set(gca,'FontSize',14)
    set(gca,'LineWidth',1.5)   

%set axis font/label and font size
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12, ...
    'FontWeight','normal', 'LineWidth', 1.5,'layer','top')



end

