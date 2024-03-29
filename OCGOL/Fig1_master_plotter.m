
%% Figure 1 Master plotter 

%import plotting data

%import lick data
load('G:\Figure_1_OCGOL_learning_long_term\shared_plotting_data\lick_data.mat')
%import speed data for fig 1 e/f
load('G:\Figure_1_OCGOL_learning_long_term\shared_plotting_data\mean_vel.mat')

%which animal to generate plot from
aa=2;
%50 bins
edges = (0:0.02:1);

%magenta = [0.85, 0.42, 0.68];
green = [50,205,50]./255;
blue = [65,105,225]./255;
red =  [220,20,60]./255;


%% master layout
fig = figure;
fig.Units = 'centimeters';
fig.Position(1) = 7;
fig.Position(2) = 0;
fig.Position(3) = 30;
fig.Position(4) = 36;

%master layout
gridSize = [2,2];
t1 = tiledlayout(fig,gridSize(1),gridSize(2),'TileSpacing','normal','Padding','compact','Units','centimeters');

%fraction of licks subplot

s1 = tiledlayout(t1,2,1,'TileSpacing','normal','Padding','compact','Units','centimeters');
s1.Layout.Tile = 3;
s1.Layout.TileSpan = [1,1];
%subplot title
%set_subplot_title(s1,'Common', 'Arial', 16,'bold')

nexttile(s1,1)
hold on
%axis square
xlim([0.5 4.5])
xticks(1:4)
xticklabels({'Random\newline foraging' ,'5A5B','3A3B','Random\newline AB'})
ylim([0 1.2])
yticks(0:0.2:1)
ylabel('Fraction of licks in reward zone')
ae = errorbar(mean_zone_A,sem_zone_A,'Color',[65,105,225]./255,'LineWidth',2);
be = errorbar(mean_zone_B,sem_zone_B,'Color',[220,20,60]./255,'LineWidth',2);

legend([ae be],{'A','B'},'Location','northwest')

% set(gca,'FontSize',16)
% set(gca,'LineWidth',2)

%line plot
nexttile(s1,2)
hold on
%axis square
xlim([0.5 4.5])
xticks(1:4)
xticklabels({'Random\newline foraging' ,'5A5B','3A3B','Random\newline AB'})
%xtickangle(45)
ylim([0 1.2])
yticks([0:0.2:1])
ylabel('Fraction of correct trials')
%mean A corr
ae = errorbar(mean_corr(:,1),sem_corr(:,1),'Color',[65,105,225]./255,'LineWidth',2);
%mean B corr
be = errorbar(mean_corr(:,2),sem_corr(:,2),'Color',[220,20,60]./255,'LineWidth',2);
%mean all corr
abe = errorbar(mean_corr(:,3),sem_corr(:,3),'Color',[139, 0, 139]./255,'LineWidth',2);

legend([ae be abe],{'A','B','All'},'Location','northwest')

% set(gca,'FontSize',16)
% set(gca,'LineWidth',2)

%lick histograms
s2 = tiledlayout(t1,4,2,'TileSpacing','normal','Padding','compact','Units','centimeters');
s2.Layout.Tile = 2;
s2.Layout.TileSpan = [1,1];

%top histogram (RF) sub layout
s21 = tiledlayout(s2,1,2,'TileSpacing','normal','Padding','compact','Units','centimeters');
s21.Layout.Tile = 1;
s21.Layout.TileSpan = [1,2];

%top lick histogram on RF
nexttile(s21,1,[1,2])
hold on
title('Random')
plotHisto(cell2mat(lick_lap_norm{aa}{1}'),edges,green);
plot_orderA = [0 3 5 7];
plot_orderB = [0 4 6 8];

%A lick histograms
for ii=2:4
    nexttile(s2,plot_orderA(ii),[1,1])
    hold on
    title('A')
    plotHisto(cell2mat(lick_lap_norm_A{aa}{ii}'),edges,blue);
    if ii~=4
    xlabel('')
    end
end

%B lick histograms
for ii=2:4
nexttile(s2,plot_orderB(ii),[1,1])
hold on
title('B')
plotHisto(cell2mat(lick_lap_norm_B{aa}{ii}'),edges,red);
    if ii~=4
    xlabel('')
    end
    %turn off ylabel
    ylabel('')
end

%lick histograms
s3 = tiledlayout(t1,4,2,'TileSpacing','normal','Padding','compact','Units','centimeters');
s3.Layout.Tile = 4;
s3.Layout.TileSpan = [1,1];

%top histogram (RF) sub layout
s31 = tiledlayout(s3,1,2,'TileSpacing','normal','Padding','compact','Units','centimeters');
s31.Layout.Tile = 1;
s31.Layout.TileSpan = [1,2];

%RF plot
nexttile(s31,1,[1,2])
hold on
ylabel({'Mean Speed' ; '[cm/s]'})
xlabel('Time [s]')
title('Random')
set(gca,'FontSize',12)
set(gca,'linewidth',2)
ylim([0 30])
xlim([0 151])
yticks([0 10 20 30])
xticks([16 46 76 106 136]);
xticklabels({'-2','-1','0','1','2'})
%line plot with std at each spatial bin
%plot(mean(mean_vel.A{2}),'b')
%shaded std
s1 = shadedErrorBar(1:151,mean_vel.A{1},{@mean,@std},'lineprops','-','transparent',true,'patchSaturation',0.20);
set(s1.edge,'LineWidth',1.5,'LineStyle','-','Color',[[65,105,225]/255, 0.2]) %last value add transparency value
s1.mainLine.LineWidth = 2;
s1.mainLine.Color = [65,105,225]/255;
s1.patch.FaceColor = [65,105,225]/255;

%line plot with std at each spatial bin
%plot(mean(mean_vel.ArelB{2}),'r')
%shaded std
s2 = shadedErrorBar(1:151,mean_vel.B{1},{@mean,@std},'lineprops','-','transparent',true,'patchSaturation',0.20);
set(s2.edge,'LineWidth',1.5,'LineStyle','-','Color',[[220,20,60]/255, 0.2]) %last value add transparency value
s2.mainLine.LineWidth = 2;
s2.mainLine.Color = [220,20,60]/255;
s2.patch.FaceColor = [220,20,60]/255;

%plot reward zones start line
plot([76,76],[0 30],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2)

%plot A lick speeds
sub_orderA = [0 3 5 7];

for ii=2:4
    nexttile(s3,sub_orderA(ii),[1,1])
    hold on
    title('A')
    ylabel({'Mean Velocity' ; '[cm/s]'})
    if ii == 4
    xlabel('Time [s]')
    end 
    set(gca,'FontSize',12)
set(gca,'linewidth',2)
    yticks([0 10 20 30])
    ylim([0 30])
    xlim([0 151])
    xticks([16 46 76 106 136]);
    xticklabels({'-2','-1','0','1','2'})
    %line plot with std at each spatial bin
    %plot(mean(mean_vel.A{2}),'b')
    %shaded std
    s1 = shadedErrorBar(1:151,mean_vel.A{ii},{@mean,@std},'lineprops','-','transparent',true,'patchSaturation',0.20);
    set(s1.edge,'LineWidth',1.5,'LineStyle','-','Color',[[65,105,225]/255, 0.2]) %last value add transparency value
    s1.mainLine.LineWidth = 2;
    s1.mainLine.Color = [65,105,225]/255;
    s1.patch.FaceColor = [65,105,225]/255;
    
    %line plot with std at each spatial bin
    %plot(mean(mean_vel.ArelB{2}),'r')
    %shaded std
    s2 = shadedErrorBar(1:151,mean_vel.ArelB{ii},{@mean,@std},'lineprops','-','transparent',true,'patchSaturation',0.20);
    set(s2.edge,'LineWidth',1.5,'LineStyle','-','Color',[[220,20,60]/255, 0.2]) %last value add transparency value
    s2.mainLine.LineWidth = 2;
    s2.mainLine.Color = [220,20,60]/255;
    s2.patch.FaceColor = [220,20,60]/255;
    
        %plot reward zones start line
plot([76,76],[0 30],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2)

end

%B trials
sub_orderB = [0 4 6 8];
for ii=2:4
    nexttile(s3,sub_orderB(ii),[1,1])
    hold on
    title('B')
    yticks([0 10 20 30])
    %ylabel('Mean Velocity [cm/s]')
    if ii == 4
    xlabel('Time [s]')
    end 
    set(gca,'FontSize',12)
    set(gca,'linewidth',2)
    ylim([0 30])
    xlim([0 151])
    xticks([16 46 76 106 136]);
    xticklabels({'-2','-1','0','1','2'})
    %line plot with std at each spatial bin
    %plot(mean(mean_vel.A{2}),'b')
    %shaded std
    s1 = shadedErrorBar(1:151,mean_vel.B{ii},{@mean,@std},'lineprops','-','transparent',true,'patchSaturation',0.20);
    set(s1.edge,'LineWidth',1.5,'LineStyle','-','Color',[[220,20,60]/255, 0.2]) %last value add transparency value
    s1.mainLine.LineWidth = 2;
    s1.mainLine.Color = [220,20,60]/255;
    s1.patch.FaceColor = [220,20,60]/255;
    
    %line plot with std at each spatial bin
    %plot(mean(mean_vel.ArelB{2}),'r')
    %shaded std
    s2 = shadedErrorBar(1:151,mean_vel.BrelA{ii},{@mean,@std},'lineprops','-','transparent',true,'patchSaturation',0.20);
    set(s2.edge,'LineWidth',1.5,'LineStyle','-','Color',[[65,105,225]/255, 0.2]) %last value add transparency value
    s2.mainLine.LineWidth = 2;
    s2.mainLine.Color = [65,105,225]/255;
    s2.patch.FaceColor = [65,105,225]/255;
    
        %plot reward zones start line
plot([76,76],[0 30],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2)
end

% Frame indices corresponding to time relative to reward zone
% 16 = -2s
% 46 = -1s
% 76 = 0s
% 106 = 1s
% 136 = 2s


%set axis font/label and font size
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12, ...
    'FontWeight','normal', 'LineWidth', 1.5,'layer','top')


%% Mean lick rate stats and plot for cell reports reviewer #1

% lick at each position according to normalized position of licks
%each animal, then each training phase, the each lap; then normalized
%position of lick
%mean sem 
%10 cm - 20 bins
%20 cm - 10 bins
lick_lap_norm
lick_lap_norm_A
lick_lap_norm_B

%cluster each vector by number of licks in bin
%animal 1, RF trials, A and B laps, all laps
%for each animal

%RF
for ii =1:4
    merge_laps = sort(cell2mat(lick_lap_norm{ii}{1}(:)));
    
    %10 bins (20 cm wide)
    den_10_RF(ii,:) = histcounts(merge_laps,0:0.1:1,'Normalization','probability');
    %20 bins (10 cm wide)
    den_20_RF(ii,:) = histcounts(merge_laps,0:0.05:1,'Normalization','probability'); 
end

%calculate mean and sem
m_10_den_RF = mean(den_10_RF,1);
sem_10_den_RF = std(den_10_RF,0,1)./sqrt(4);

m_20_den_RF = mean(den_20_RF,1);
sem_20_den_RF = std(den_20_RF,0,1)./sqrt(4);

figure
hold on
ylim([0,0.5])
errorbar(m_10_den_RF,sem_10_den_RF)

figure
hold on
ylim([0,0.5])
errorbar(m_20_den_RF,sem_20_den_RF)

%for remainder of training trials except RF
for jj = 2:4
    %5A5B (A and B trials)
    for ii =1:4
        merge_laps_A = sort(cell2mat(lick_lap_norm_A{ii}{jj}(:)));
        merge_laps_B = sort(cell2mat(lick_lap_norm_B{ii}{jj}(:)));
        
        %A
        %10 bins (20 cm wide)
        den_10_A{jj}(ii,:) = histcounts(merge_laps_A,0:0.1:1,'Normalization','probability');
        %20 bins (10 cm wide)
        den_20_A{jj}(ii,:) = histcounts(merge_laps_A,0:0.05:1,'Normalization','probability');
        %B
         %10 bins (20 cm wide)
        den_10_B{jj}(ii,:) = histcounts(merge_laps_B,0:0.1:1,'Normalization','probability');
        %20 bins (10 cm wide)
        den_20_B{jj}(ii,:) = histcounts(merge_laps_B,0:0.05:1,'Normalization','probability');       
    end
    %A mean and sem
    m_10_den_A{jj} = mean(den_10_A{jj},1);
    sem_10_den_A{jj} = std(den_10_A{jj},0,1)./sqrt(4);
    
    m_20_den_A{jj} = mean(den_20_A{jj},1);
    sem_20_den_A{jj} = std(den_20_A{jj},0,1)./sqrt(4);
    
    
    %b mean and sem
    m_10_den_B{jj} = mean(den_10_B{jj},1);
    sem_10_den_B{jj} = std(den_10_B{jj},0,1)./sqrt(4);
    
    m_20_den_B{jj} = mean(den_20_B{jj},1);
    sem_20_den_B{jj} = std(den_20_B{jj},0,1)./sqrt(4);
end

figure('Position', [1320, 256, 404, 341])
hold on
ylim([0,0.7])
errorbar(m_10_den_A{2},sem_10_den_A{2})


figure
hold on
ylim([0,0.7])
errorbar(m_20_den_B{2},sem_20_den_B{2})

%% Mean lick relative to start of each reward zone bin for A and B trials

pos = [1311.7       41.667       940.33       599.33];
%pos = [1320, 256, 404, 341];

fig = figure('Position',pos);
%master layout
gridSize = [2,2];
t2 = tiledlayout(fig,gridSize(1),gridSize(2),'TileSpacing','normal','Padding','compact','Units','centimeters');

%rotate the vectors so that the bin start is in the center of the vector 
%10 bins - A start - 4 and B start 8
%20 bins - A start 7 and B start - 15
a10rot = -3;
b10rot = 1;

%first cell - nb bins, 10, 20
%second cell - vectors corresponding to training day
%A trials
m_den_rotA.A{1}{1} = circshift(m_10_den_RF,a10rot);
m_den_rotA.A{1}{2} = circshift(m_10_den_A{2},a10rot);
m_den_rotA.A{1}{3} = circshift(m_10_den_A{3},a10rot);
m_den_rotA.A{1}{4} = circshift(m_10_den_A{4},a10rot);

s_den_rotA.A{1}{1} = circshift(sem_10_den_RF,a10rot);
s_den_rotA.A{1}{2} = circshift(sem_10_den_A{2},a10rot);
s_den_rotA.A{1}{3} = circshift(sem_10_den_A{3},a10rot);
s_den_rotA.A{1}{4} = circshift(sem_10_den_A{4},a10rot);

%A trials - b reward zone
m_den_rotB.A{1}{1} = circshift(m_10_den_RF,b10rot);
m_den_rotB.A{1}{2} = circshift(m_10_den_A{2},b10rot);
m_den_rotB.A{1}{3} = circshift(m_10_den_A{3},b10rot);
m_den_rotB.A{1}{4} = circshift(m_10_den_A{4},b10rot);

s_den_rotB.A{1}{1} = circshift(sem_10_den_RF,b10rot);
s_den_rotB.A{1}{2} = circshift(sem_10_den_A{2},b10rot);
s_den_rotB.A{1}{3} = circshift(sem_10_den_A{3},b10rot);
s_den_rotB.A{1}{4} = circshift(sem_10_den_A{4},b10rot);

%B Trial - A reward zone
m_den_rotA.B{1}{1} = circshift(m_10_den_RF,a10rot);
m_den_rotA.B{1}{2} = circshift(m_10_den_B{2},a10rot);
m_den_rotA.B{1}{3} = circshift(m_10_den_B{3},a10rot);
m_den_rotA.B{1}{4} = circshift(m_10_den_B{4},a10rot);

s_den_rotA.B{1}{1} = circshift(sem_10_den_RF,a10rot);
s_den_rotA.B{1}{2} = circshift(sem_10_den_B{2},a10rot);
s_den_rotA.B{1}{3} = circshift(sem_10_den_B{3},a10rot);
s_den_rotA.B{1}{4} = circshift(sem_10_den_B{4},a10rot);

%A trials - b reward zone
m_den_rotB.B{1}{1} = circshift(m_10_den_RF,b10rot);
m_den_rotB.B{1}{2} = circshift(m_10_den_B{2},b10rot);
m_den_rotB.B{1}{3} = circshift(m_10_den_B{3},b10rot);
m_den_rotB.B{1}{4} = circshift(m_10_den_B{4},b10rot);

s_den_rotB.B{1}{1} = circshift(sem_10_den_RF,b10rot);
s_den_rotB.B{1}{2} = circshift(sem_10_den_B{2},b10rot);
s_den_rotB.B{1}{3} = circshift(sem_10_den_B{3},b10rot);
s_den_rotB.B{1}{4} = circshift(sem_10_den_B{4},b10rot);

a20rot = -5;
b20rot = 3;

%20 bin equivalent
%first cell - nb bins, 10, 20
%second cell - vectors corresponding to training day
%A trials
m_den_rotA.A{2}{1} = circshift(m_20_den_RF,a20rot);
m_den_rotA.A{2}{2} = circshift(m_20_den_A{2},a20rot);
m_den_rotA.A{2}{3} = circshift(m_20_den_A{3},a20rot);
m_den_rotA.A{2}{4} = circshift(m_20_den_A{4},a20rot);

s_den_rotA.A{2}{1} = circshift(sem_20_den_RF,a20rot);
s_den_rotA.A{2}{2} = circshift(sem_20_den_A{2},a20rot);
s_den_rotA.A{2}{3} = circshift(sem_20_den_A{3},a20rot);
s_den_rotA.A{2}{4} = circshift(sem_20_den_A{4},a20rot);

%A trials - b reward zone
m_den_rotB.A{2}{1} = circshift(m_20_den_RF,b20rot);
m_den_rotB.A{2}{2} = circshift(m_20_den_A{2},b20rot);
m_den_rotB.A{2}{3} = circshift(m_20_den_A{3},b20rot);
m_den_rotB.A{2}{4} = circshift(m_20_den_A{4},b20rot);

s_den_rotB.A{2}{1} = circshift(sem_20_den_RF,b20rot);
s_den_rotB.A{2}{2} = circshift(sem_20_den_A{2},b20rot);
s_den_rotB.A{2}{3} = circshift(sem_20_den_A{3},b20rot);
s_den_rotB.A{2}{4} = circshift(sem_20_den_A{4},b20rot);

%B Trial - A reward zone
m_den_rotA.B{2}{1} = circshift(m_20_den_RF,a20rot);
m_den_rotA.B{2}{2} = circshift(m_20_den_B{2},a20rot);
m_den_rotA.B{2}{3} = circshift(m_20_den_B{3},a20rot);
m_den_rotA.B{2}{4} = circshift(m_20_den_B{4},a20rot);

s_den_rotA.B{2}{1} = circshift(sem_20_den_RF,a20rot);
s_den_rotA.B{2}{2} = circshift(sem_20_den_B{2},a20rot);
s_den_rotA.B{2}{3} = circshift(sem_20_den_B{3},a20rot);
s_den_rotA.B{2}{4} = circshift(sem_20_den_B{4},a20rot);

%A trials - b reward zone
m_den_rotB.B{2}{1} = circshift(m_20_den_RF,b20rot);
m_den_rotB.B{2}{2} = circshift(m_20_den_B{2},b20rot);
m_den_rotB.B{2}{3} = circshift(m_20_den_B{3},b20rot);
m_den_rotB.B{2}{4} = circshift(m_20_den_B{4},b20rot);

s_den_rotB.B{2}{1} = circshift(sem_20_den_RF,b20rot);
s_den_rotB.B{2}{2} = circshift(sem_20_den_B{2},b20rot);
s_den_rotB.B{2}{3} = circshift(sem_20_den_B{3},b20rot);
s_den_rotB.B{2}{4} = circshift(sem_20_den_B{4},b20rot);

blue_grad = cbrewer('seq','Blues',16);
red_grad = cbrewer('seq','Reds',16);

%plot relative to A reward zone
nexttile(t2,1,[1,1])
hold on
%axis square
title('A trials')
for ii=1:4
ae(ii) = errorbar((2:8),m_den_rotA.A{1}{ii}(2:8),s_den_rotA.A{1}{ii}(2:8),'Color',blue_grad(ii*3,:),'LineWidth',1);
end

for ii=1:4
be(ii) = errorbar((2:8),m_den_rotB.A{1}{ii}(2:8),s_den_rotB.A{1}{ii}(2:8),'Color',red_grad(ii*3,:), 'LineStyle', '--','LineWidth',1)
end
ylim([0,1])
xlim([1.5,8.5])
yticks([0 0.5 1])
ylabel('Lick fraction')
xticks([5])
xticklabels({'Reward bin'})
xlabel('Spatial bin [10 bins]')
al = legend([ae,be],{'RF','5A5B','3A3B','Random',...
    'RF','5A5B','3A3B','Random'},'Location','northeast');

nexttile(t2,2,[1,1])
hold on
%axis square
title('B trials')
for ii=1:4
ae(ii) = errorbar((2:8),m_den_rotA.B{1}{ii}(2:8),s_den_rotA.B{1}{ii}(2:8),'Color',blue_grad(ii*3,:), 'LineStyle', '--','LineWidth',1);
end

for ii=1:4
be(ii) = errorbar((2:8),m_den_rotB.B{1}{ii}(2:8),s_den_rotB.B{1}{ii}(2:8),'Color',red_grad(ii*3,:), 'LineStyle', '-','LineWidth',1)
end
ylim([0,1])
xlim([1.5,8.5])
yticks([0 0.5 1])
ylabel('Lick fraction')
xticks([5])
xticklabels({'Reward bin'})
xlabel('Spatial bin [10 bins]')
al = legend([ae,be],{'RF','5A5B','3A3B','Random',...
    'RF','5A5B','3A3B','Random'},'Location','northeast');

%20 bin plots
nexttile(t2,3,[1,1])
hold on
%axis square
title('A trials')
for ii=1:4
ae(ii) = errorbar((4:16),m_den_rotA.A{2}{ii}(4:16),s_den_rotA.A{2}{ii}(4:16),'Color',blue_grad(ii*3,:),'LineWidth',1);
end

for ii=1:4
be(ii) = errorbar((4:16),m_den_rotB.A{2}{ii}(4:16),s_den_rotB.A{2}{ii}(4:16),'Color',red_grad(ii*3,:), 'LineStyle', '--','LineWidth',1);
end
ylim([0,1])
xlim([3.5,16.5])
yticks([0 0.5 1])
ylabel('Lick fraction')
xticks([10])
xticklabels({'Reward bin'})
xlabel('Spatial bin [20 bins]')
al = legend([ae,be],{'RF','5A5B','3A3B','Random',...
    'RF','5A5B','3A3B','Random'},'Location','northeast');

nexttile(t2,4,[1,1])
hold on
%axis square
title('B trials')
for ii=1:4
ae(ii) = errorbar((4:16),m_den_rotA.B{2}{ii}(4:16),s_den_rotA.B{2}{ii}(4:16),'Color',blue_grad(ii*3,:), 'LineStyle', '--','LineWidth',1);
end

for ii=1:4
be(ii) = errorbar((4:16),m_den_rotB.B{2}{ii}(4:16),s_den_rotB.B{2}{ii}(4:16),'Color',red_grad(ii*3,:), 'LineStyle', '-','LineWidth',1);
end
ylim([0,1])
xlim([3.5,16.5])
yticks([0 0.5 1])
ylabel('Lick fraction')
xticks([10])
xticklabels({'Reward bin'})
xlabel('Spatial bin [20 bins]')
al = legend([ae,be],{'RF','5A5B','3A3B','Random',...
    'RF','5A5B','3A3B','Random'},'Location','northeast');

%set axis font/label and font size
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12, ...
    'FontWeight','normal', 'LineWidth', 1.5,'layer','top')

%% Mean lick plotter code (revision)master layout - mean lick (10 bins)

fig2 = figure('Position',[1311.7       41.667       473.33       599.33]);
% fig.Units = 'centimeters';
% fig.Position(1) = 7;
% fig.Position(2) = 0;
% fig.Position(3) = 30;
% fig.Position(4) = 36;

blue_color = [65,105,225]./255;
red_color = [220,20,60]./255;
green_color = [0,110,51]./255;

%master layout
gridSize = [4,2];
t2 = tiledlayout(fig2,gridSize(1),gridSize(2),'TileSpacing','normal','Padding','compact','Units','centimeters');

%lick histograms
s2 = tiledlayout(t2,1,2,'TileSpacing','compact','Padding','compact','Units','centimeters');
s2.Layout.Tile = 1;
s2.Layout.TileSpan = [1,2];

%RF
nexttile(s2,[1 2])
hold on
title('RF')
%axis square
xlim([0.5 10.5])
xlabel([{'Spatial bins'},{'(N = 10 bins)'}])
xticks([1,10])
%xticklabels({'Spatial bins (N = 10 bins)'})
ylim([0 0.6])
yticks(0:0.2:0.5)
ylabel('Fraction of licks')
errorbar(m_10_den_RF,sem_10_den_RF,'Color',green_color,'LineWidth',1.5);

%B reward zone line
xline(4,'Color',red_color,'LineWidth',1,'LineStyle', '--')
%A reward zone line
xline(8,'Color',blue_color,'LineWidth',1,'LineStyle', '--')

%A trials
seq_order = 2:4;
plot_order = [3,5,7];
for ii=1:3
nexttile(t2,plot_order(ii),[1 1])
hold on
if ii==1
title('A')
end
xlim([0.5 10.5])
if ii == 3
xlabel([{'Spatial bins'},{'(N = 10 bins)'}])
end
xticks([1,10])
ylim([0 1])
yticks(0:0.5:1)
ylabel('Fraction of licks')
errorbar(m_10_den_A{seq_order(ii)},sem_10_den_A{seq_order(ii)},'Color',blue_color,'LineWidth',1.5)

%B reward zone line
xline(4,'Color',red_color,'LineWidth',1,'LineStyle', '--')
%A reward zone line
xline(8,'Color',blue_color,'LineWidth',1,'LineStyle', '--')
end


%B trials
seq_order = 2:4;
plot_order = [4,6,8];
for ii=1:3
nexttile(t2,plot_order(ii),[1 1])
hold on
if ii==1
title('B')
end
xlim([0.5 10.5])
if ii == 3
xlabel([{'Spatial bins'},{'(N = 10 bins)'}])
end
xticks([1,10])
ylim([0 1])
yticks(0:0.5:1)
%ylabel('Fraction of licks')
errorbar(m_10_den_B{seq_order(ii)},sem_10_den_B{seq_order(ii)},'Color',red_color,'LineWidth',1.5)
%B reward zone line
xline(4,'Color',red_color,'LineWidth',1,'LineStyle', '--')
%A reward zone line
xline(8,'Color',blue_color,'LineWidth',1,'LineStyle', '--')

end

%set axis font/label and font size
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12, ...
    'FontWeight','normal', 'LineWidth', 1.5,'layer','top')

%% Bin speed plotter code (revision)master layout - mean lick (10 bins)

fig2 = figure('Position',[1311.7       41.667       473.33       599.33]);
% fig.Units = 'centimeters';
% fig.Position(1) = 7;
% fig.Position(2) = 0;
% fig.Position(3) = 30;
% fig.Position(4) = 36;

blue_color = [65,105,225]./255;
red_color = [220,20,60]./255;
green_color = [0,110,51]./255;

%master layout
gridSize = [4,2];
t2 = tiledlayout(fig2,gridSize(1),gridSize(2),'TileSpacing','normal','Padding','compact','Units','centimeters');

%lick histograms
s2 = tiledlayout(t2,1,2,'TileSpacing','compact','Padding','compact','Units','centimeters');
s2.Layout.Tile = 1;
s2.Layout.TileSpan = [1,2];

%RF
nexttile(s2,[1 2])
hold on
title('RF')
%axis square
xlim([0.5 20.5])
xlabel([{'Spatial bins'},{'(N = 20 bins)'}])
xticks([1,20])
%xticklabels({'Spatial bins (N = 10 bins)'})
ylim([0 0.6])
yticks(0:0.2:0.5)
ylabel('Fraction of licks')
errorbar(m_20_den_RF,sem_20_den_RF,'Color',green_color,'LineWidth',1.5);

%B reward zone line
xline(7,'Color',red_color,'LineWidth',1,'LineStyle', '--')
%A reward zone line
xline(15,'Color',blue_color,'LineWidth',1,'LineStyle', '--')

%A trials
seq_order = 2:4;
plot_order = [3,5,7];
for ii=1:3
nexttile(t2,plot_order(ii),[1 1])
hold on
if ii==1
title('A')
end
xlim([0.5 20.5])
if ii == 3
xlabel([{'Spatial bins'},{'(N = 20 bins)'}])
end
xticks([1,20])
ylim([0 1])
yticks(0:0.5:1)
ylabel('Fraction of licks')
errorbar(m_20_den_A{seq_order(ii)},sem_20_den_A{seq_order(ii)},'Color',blue_color,'LineWidth',1.5)

%B reward zone line
xline(7,'Color',red_color,'LineWidth',1,'LineStyle', '--')
%A reward zone line
xline(15,'Color',blue_color,'LineWidth',1,'LineStyle', '--')
end

%B trials
seq_order = 2:4;
plot_order = [4,6,8];
for ii=1:3
nexttile(t2,plot_order(ii),[1 1])
hold on
if ii==1
title('B')
end
xlim([0.5 20.5])
if ii == 3
xlabel([{'Spatial bins'},{'(N = 20 bins)'}])
end
xticks([1,20])
ylim([0 1])
yticks(0:0.5:1)
%ylabel('Fraction of licks')
errorbar(m_20_den_B{seq_order(ii)},sem_20_den_B{seq_order(ii)},'Color',red_color,'LineWidth',1.5)
%B reward zone line
xline(7,'Color',red_color,'LineWidth',1,'LineStyle', '--')
%A reward zone line
xline(15,'Color',blue_color,'LineWidth',1,'LineStyle', '--')

end

%set axis font/label and font size
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12, ...
    'FontWeight','normal', 'LineWidth', 1.5,'layer','top')
