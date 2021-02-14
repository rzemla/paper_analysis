function [outputArg1,outputArg2] = fig_sup_speed_place_field_master(event_speed_plot,pf_prop_data)

%% Unload data

%speed scatterplot data
Asel_speed_diff = event_speed_plot.Asel_speed_diff;
Bsel_speed_diff = event_speed_plot.Bsel_speed_diff;
Asel_speed_cum = event_speed_plot.Asel_speed_cum;
Bsel_speed_cum = event_speed_plot.Bsel_speed_cum;

edges = [10:2:70];
edge_center = edges(1:end-1)+1;

%place field properties data
cum_prob_AB_A = pf_prop_data.cum_prob_AB_A;
cum_prob_AB_B = pf_prop_data.cum_prob_AB_B;
cum_prob_A_sel = pf_prop_data.cum_prob_A_sel;
cum_prob_B_sel = pf_prop_data.cum_prob_B_sel;

sem_cum_prob_AB_A = pf_prop_data.sem_cum_prob_AB_A;
sem_cum_prob_AB_B = pf_prop_data.sem_cum_prob_AB_B;
sem_cum_prob_A = pf_prop_data.sem_cum_prob_A;
sem_cum_prob_B = pf_prop_data.sem_cum_prob_B;

grouped_norm_mean = pf_prop_data.grouped_norm_mean;
grouped_norm_sem = pf_prop_data.grouped_norm_sem;

%% Master plot
fig = figure;
fig.Units = 'centimeters';
fig.Position(1) = 7;
fig.Position(2) = 0;
fig.Position(3) = 24;
fig.Position(4) = 36;

%master layout
gridSize = [3,2];
t1 = tiledlayout(fig,gridSize(1),gridSize(2),'TileSpacing','normal','Padding','compact','Units','centimeters');

%color scheme
color_mat = [65,105,225; 220,20,60; 135,206,250; 240,128,128]./255;


% Plot scatters 
nexttile(t1,1)
hold on
title('A selective')
axis square
xlim([0 30])
ylim([0 30])
xlabel('A lap speed [cm/s]')
ylabel('B lap speed [cm/s]')
scatter(Asel_speed_cum(:,1),Asel_speed_cum(:,2),6,'filled','MarkerFaceColor',color_mat(1,:))
plot([0 30],[0 30],'k--')

set(gca,'FontSize',12)
set(gca,'LineWidth',1.5)

nexttile(t1,2)
hold on
title('B selective')
axis square
xlim([0 30])
ylim([0 30])
xlabel('A lap speed [cm/s]')
ylabel('B lap speed [cm/s]')
scatter(Bsel_speed_cum(:,1),Bsel_speed_cum(:,2),6,'filled','MarkerFaceColor',color_mat(2,:))
plot([0 30],[0 30],'k--')

set(gca,'FontSize',12)
set(gca,'LineWidth',1.5)

%plot histogram

nexttile(t1,3)
hold on
title('A selective')
axis square
yticks([0 0.1 0.2])
ylabel('Fraction of neurons')
xlim([-20 20])
hA = histogram(Asel_speed_diff,[-20:2:20],'Normalization','probability');
hA.FaceColor = color_mat(1,:);
xlabel('Speed difference (A-B) [cm/s]')
%5 cm/s cutoffs
plot([-5 -5],[0 0.25],'k','LineWidth',1)
plot([5 5],[0 0.25],'k','LineWidth',1)

set(gca,'FontSize',12)
set(gca,'LineWidth',1.5)

nexttile(t1,4)
hold on
title('B selective')
axis square
yticks([0 0.1 0.2])
ylabel('Fraction of neurons')
xlim([-20 20])
hA = histogram(Bsel_speed_diff,[-20:2:20],'Normalization','probability');
hA.FaceColor = color_mat(2,:);
xlabel('Speed difference (A-B) [cm/s]')
%5 cm/s cutoffs
plot([-5 -5],[0 0.25],'k','LineWidth',1)
plot([5 5],[0 0.25],'k','LineWidth',1)

set(gca,'FontSize',12)
set(gca,'LineWidth',1.5)

%plot bar
nexttile(t1,5)
hold on;
axis square
%title('Place fields per neuron - selective and A&B (both crit)');
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
xlabel('Number of place fields');
ylabel('Fraction of neurons');
legend('A selective','B selective','A&B - A', 'A&B - B')

set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)



%linewidth
lw = 2;
nexttile(t1,6)
hold on
axis square
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



%set axis font/label and font size
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12, ...
    'FontWeight','normal', 'LineWidth', 1.5,'layer','top')


end

