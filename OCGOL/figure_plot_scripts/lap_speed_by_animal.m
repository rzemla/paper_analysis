function [outputArg1,outputArg2] = lap_speed_by_animal(path_dir)

%% Load in lap speed data

%read relevant data
for ee=1:size(path_dir,2)
    load_data_path{ee} = fullfile(path_dir{ee},'cumul_analysis','lap_and_event_speed.mat');
    speed_data{ee} = load(string(load_data_path{ee}),'mean_bin_speed');
end

%% Generate mean difference between bins

%get mean across each bin for A and B laps
for ii=1:size(path_dir,2)
    %get mean for A and B animal
    %A laps
    mean_speed.A(ii,:) = nanmean(speed_data{ii}.mean_bin_speed.A,1);
    %B laps
    mean_speed.B(ii,:) = nanmean(speed_data{ii}.mean_bin_speed.B,1);
end

%get mean difference between A and B laps
mean_diff = mean_speed.A - mean_speed.B; 


%% Reviewer speed analysis
%for each animal
%split into bin 30-35 - B zone (complete zone - 25-40)
%split into bin 70-75 - A zone ( complete zone - 65-80)

pos = [2205	169.800000000000	940.800000000000	973.600000000000];
%pos = [1320, 256, 404, 341];

fig = figure('Position',pos);
%master layout
gridSize = [11,4];
t2 = tiledlayout(fig,gridSize(1),gridSize(2),'TileSpacing','normal','Padding','compact','Units','centimeters');
%mean +/- sem with color shading plot for each animal
%11 x 4 tabular graph
%create input function to calculate standard error of mean
sem = @(x) nanstd(x,0,1)./sqrt(size(x,1));
list_order = (1:44);
list_order = reshape(list_order,[4,11])';

% A reward zone, A trials
for ii=1:11
    nexttile(t2,list_order_A(ii,1),[1,1])
    s1 = shadedErrorBar(60:85,speed_data{ii}.mean_bin_speed.A(:,60:85),{@nanmean,@nanstd},'lineprops','-','transparent',true,'patchSaturation',0.20);
    set(s1.edge,'LineWidth',0.2,'LineStyle','-','Color',[[65,105,225]/255, 0.2]) %last value add transparency value
    s1.mainLine.LineWidth = 2;
    s1.mainLine.Color = [65,105,225]/255;
    s1.patch.FaceColor = [65,105,225]/255;
    ylim([0 30])
    yticks([0 15 30])
    xline(70,'Color',[65,105,225]/255);
    xline(74,'Color',[65,105,225]/255);
end
%B reward zone, A trials
for ii=1:11
    nexttile(t2,list_order_A(ii,2),[1,1])
    s1 = shadedErrorBar(20:45,speed_data{ii}.mean_bin_speed.A(:,20:45),{@nanmean,@nanstd},'lineprops','-','transparent',true,'patchSaturation',0.20);
    set(s1.edge,'LineWidth',0.2,'LineStyle','-','Color',[[65,105,225]/255, 0.2]) %last value add transparency value
    s1.mainLine.LineWidth = 2;
    s1.mainLine.Color = [65,105,225]/255;
    s1.patch.FaceColor = [65,105,225]/255;
    ylim([0 30])
    yticks([])
    xline(30,'Color',[220,20,60]/255);
    xline(34,'Color',[220,20,60]/255);

end
%B reward zone, B trials
for ii=1:11
    nexttile(t2,list_order_A(ii,3),[1,1])
    s1 = shadedErrorBar(20:45,speed_data{ii}.mean_bin_speed.B(:,20:45),{@nanmean,@nanstd},'lineprops','-','transparent',true,'patchSaturation',0.20);
    set(s1.edge,'LineWidth',0.2,'LineStyle','-','Color',[[220,20,60]/255, 0.2]) %last value add transparency value
    s1.mainLine.LineWidth = 2;
    s1.mainLine.Color = [220,20,60]/255;
    s1.patch.FaceColor = [220,20,60]/255;
    ylim([0 30])
    xline(30,'Color',[220,20,60]/255);
    xline(34,'Color',[220,20,60]/255);
end

%A reward zone, A trials
for ii=1:11
    nexttile(t2,list_order_A(ii,4),[1,1])
    s2 = shadedErrorBar(60:85,speed_data{ii}.mean_bin_speed.B(:,60:85),{@nanmean,@nanstd},'lineprops','-','transparent',true,'patchSaturation',0.20);
    set(s2.edge,'LineWidth',0.2,'LineStyle','-','Color',[[220,20,60]/255, 0.2]) %last value add transparency value
    s2.mainLine.LineWidth = 2;
    s2.mainLine.Color = [220,20,60]/255;
    s2.patch.FaceColor = [220,20,60]/255;
    ylim([0 30])
    xline(70,'Color',[65,105,225]/255);
    xline(74,'Color',[65,105,225]/255);
end

%% Generate subplot for all animals with mean speed on A/B laps across all bins

figure('Position',[2308 118 1148 861])
for ii=1:size(path_dir,2)
    subplot(3,4,ii)
    hold on
    axis square
    title(num2str(ii))
    ylim([0 40])
    yticks([0 10 20 30]);
    xlabel('Normalized position');
    xticks([1 50 100])
    xticklabels({'0','0,5','1'});
    ylabel('Speed [cm/s]');
    
    %create input function to calculate standard error of mean
    sem = @(x) nanstd(x,0,1)./sqrt(size(x,1));
    
    s1 = shadedErrorBar(1:100,speed_data{ii}.mean_bin_speed.A,{@nanmean,@nanstd},'lineprops','-','transparent',true,'patchSaturation',0.20);
    set(s1.edge,'LineWidth',0.2,'LineStyle','-','Color',[[65,105,225]/255, 0.2]) %last value add transparency value
    s1.mainLine.LineWidth = 2;
    s1.mainLine.Color = [65,105,225]/255;
    s1.patch.FaceColor = [65,105,225]/255;
    
    s2 = shadedErrorBar(1:100,speed_data{ii}.mean_bin_speed.B,{@nanmean,@nanstd},'lineprops','-','transparent',true,'patchSaturation',0.20);
    set(s2.edge,'LineWidth',0.2,'LineStyle','-','Color',[[220,20,60]/255, 0.2]) %last value add transparency value
    s2.mainLine.LineWidth = 2;
    s2.mainLine.Color = [220,20,60]/255;
    s2.patch.FaceColor = [220,20,60]/255;
    
    %plot legend only for first animal
    if ii==1
    legend('A laps','B laps','location','northeast','AutoUpdate','off');
    end
    
    %plot reward zone A and B markers
    %reward zone A
    plot([70 70],[0 40],'--','Color',[65,105,225]/255)
    %reward zone B
    plot([30 30],[0 40],'--','Color',[220,20,60]/255)
    
    set(gca,'FontSize',11)
    set(gca,'LineWidth',1.5)
end

subplot(3,4,12)
hold on
axis square
%title(num2str(ii))
ylim([-15 15])
yticks([-15 -10 -5 0 5 10 15]);
xlabel('Normalized position');
xticks([1 50 100])
xticklabels({'0','0,5','1'});
ylabel('Mean Speed Difference (A-B) \newline  [cm/s]');
s1 = shadedErrorBar(1:100,mean_diff,{@nanmean,@nanstd},'lineprops','-','transparent',true,'patchSaturation',0.20);
set(s1.edge,'LineWidth',1.5,'LineStyle','-','Color',[[34,139,34]/255, 0.2]) %last value add transparency value
s1.mainLine.LineWidth = 2;
s1.mainLine.Color = [34,139,34]./255;
s1.patch.FaceColor = [34,139,34]./255;

%plot reward zone A and B markers
%reward zone A
plot([70 70],[-15 15],'--','Color',[65,105,225]/255)
%reward zone B
plot([30 30],[-15 15],'--','Color',[220,20,60]/255)

%plot 5cm/s difference lines
plot([1 100],[5 5],'-','Color','k')
plot([1 100],[-5 -5],'-','Color','k')

set(gca,'FontSize',11)
set(gca,'LineWidth',1.5)


end

