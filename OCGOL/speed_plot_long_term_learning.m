%% File directories
%I45_RT
pathDir{1}{1} = 'G:\OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_rand_d1_052218';
pathDir{1}{2} = 'G:\OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_5A5B_053018';
pathDir{1}{3} = 'G:\OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_3A3B_060518';
pathDir{1}{4} = 'G:\OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_AB_061418';

%I46
pathDir{2}{1} = 'G:\OCGOL_learning_long_term\I46\behavior_only\I46_rand_d1_052918';
pathDir{2}{2} = 'G:\OCGOL_learning_long_term\I46\behavior_only\I46_5A5B_060118';
pathDir{2}{3} = 'G:\OCGOL_learning_long_term\I46\behavior_only\I46_3A3B_060718';
pathDir{2}{4} = 'G:\OCGOL_learning_long_term\I46\behavior_only\I46_AB_061518';

%I47_RS
pathDir{3}{1} = 'G:\OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_rand_d2_051518';
pathDir{3}{2} = 'G:\OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_5AB_d7_052218';
pathDir{3}{3} = 'G:\OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_3AB_d8_052418';
pathDir{3}{4} = 'G:\OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_AB_061418';

%I47_LP
pathDir{4}{1} = 'G:\OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_rand_d2_051518';
pathDir{4}{2} = 'G:\OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_5AB_d1_051718';
pathDir{4}{3} = 'G:\OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_3AB_d8_052418';
pathDir{4}{4} = 'G:\OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_AB_061418';

%% Read in velocity related data for each animal and construct related data

%for each animal
for aa=1:4
    %for each session
    for ss=1:4
        speed_data{aa}{ss} = load(fullfile(pathDir{aa}{ss},'cumul_analysis','speed_data.mat'));
    end
end

%extract lap idx's with A or B trials for each animal (2-4 sessions only)
for aa=1:4
    for ss=2:4
        %lap idx by trial type
        speed_data{aa}{ss}.lap_idx.A = find(speed_data{aa}{ss}.Behavior.lap_id.trial_based == 2);
        speed_data{aa}{ss}.lap_idx.B = find(speed_data{aa}{ss}.Behavior.lap_id.trial_based == 3);
    end
end

%% Get mean velocity from each animal in peri-reward region on A and B trials

%range over which to average velocity
idx_width = [75,75];

%for each animal
for aa=1:4
    %for each OCGOL session
    for ss=2:4
        %for each lap (A laps)
        for ll=1:size(speed_data{aa}{ss}.lap_idx.A,2)
            %plot line plot along range
            A_speed{aa}{ss}(ll,:) = speed_data{aa}{ss}.rewards.A.speed{ll}((speed_data{aa}{ss}.rewards.A.Imin(ll)-idx_width(1)):(speed_data{aa}{ss}.rewards.A.Imin(ll)+idx_width(2)));
            A_speed_relB{aa}{ss}(ll,:) = speed_data{aa}{ss}.rewards.A.speed{ll}((speed_data{aa}{ss}.rewards.A.IminB(ll)-idx_width(1)):(speed_data{aa}{ss}.rewards.A.IminB(ll)+idx_width(2)));
        end
        %for each lap (B laps)
        for ll=1:size(speed_data{aa}{ss}.lap_idx.B,2)
            %plot line plot along range
            B_speed{aa}{ss}(ll,:) = speed_data{aa}{ss}.rewards.B.speed{ll}((speed_data{aa}{ss}.rewards.B.Imin(ll)-idx_width(1)):(speed_data{aa}{ss}.rewards.B.Imin(ll)+idx_width(2)));
            B_speed_relA{aa}{ss}(ll,:) = speed_data{aa}{ss}.rewards.B.speed{ll}((speed_data{aa}{ss}.rewards.B.IminA(ll)-idx_width(1)):(speed_data{aa}{ss}.rewards.B.IminA(ll)+idx_width(2)));
        end
    end
end


%% Calculate the mean velocities for the RF day relative to mean reward location on the 5A5B session for each animal

%find the min index corresponding to each RF lap to mean A reward onset
%position and mean B reward onset position
for aa=1:4
    for ss=1
        for ll = 1:size(speed_data{aa}{ss}.Behavior.lap,2)
            %target A reward zones
            rf{aa}.A.res_idx{ll} = find(speed_data{aa}{ss}.lap_idx_resampled == ll);
            
            [~,rf{aa}.A.Imin(ll)] = min(abs(speed_data{aa}{ss}.position(rf{aa}.A.res_idx{ll}) - speed_data{aa}{2}.rewards.A.pos_mean));
            
            rf{aa}.A.lap_position{ll} = speed_data{aa}{ss}.position(rf{aa}.A.res_idx{ll});
            
            rf{aa}.A.speed{ll} = speed_data{aa}{ss}.speed(rf{aa}.A.res_idx{ll});
                        
            % B reward zones
            [~,rf{aa}.B.Imin(ll)] = min(abs(speed_data{aa}{ss}.position(rf{aa}.A.res_idx{ll}) - speed_data{aa}{2}.rewards.B.pos_mean));
            %extract binned indices corresponding to A laps
            %rewards.A.pos_bins{ll} = pos_bins(rewards.A.res_idx{ll});
        end
    end
end

%extract speed data in peri-reward A zone and peri-reward B zone
%for each animal
for aa=1:4
    %for each OCGOL session
    for ss=1
        %for each lap
        for ll = 1:size(speed_data{aa}{ss}.Behavior.lap,2)
            A_speed{aa}{ss}(ll,:) = speed_data{aa}{ss}.split_lap.speed{ll}((rf{aa}.A.Imin(ll)-idx_width(1)):...
                (rf{aa}.A.Imin(ll)+idx_width(2)));
            B_speed{aa}{ss}(ll,:) = speed_data{aa}{ss}.split_lap.speed{ll}((rf{aa}.B.Imin(ll)-idx_width(1)):...
                (rf{aa}.B.Imin(ll)+idx_width(2)));
        end
    end
end

%get mean velocities in per-reward area
for aa=1:4
    %for each OCGOL session
    for ss=1:4
        %get mean velocity at each timepoint for each animal
        %A trials
        %mean
        mean_vel.A{ss}(aa,:) = mean(A_speed{aa}{ss},1);
        mean_vel.ArelB{ss}(aa,:) = mean(A_speed_relB{aa}{ss},1);
        %std
        std_vel.A{ss}(aa,:) = std(A_speed{aa}{ss},0,1);
        std_vel.ArelB{ss}(aa,:) = std(A_speed_relB{aa}{ss},0,1);
        %B trials
        mean_vel.B{ss}(aa,:) = mean(B_speed{aa}{ss},1);
        mean_vel.BrelA{ss}(aa,:) = mean(B_speed_relA{aa}{ss},1);
        %std
        std_vel.B{ss}(aa,:) = std(B_speed{aa}{ss},0,1);
        std_vel.BrelA{ss}(aa,:) = std(B_speed_relA{aa}{ss},0,1);
    end
end

%% Plot
figure('Position', [1150 85 534 830])

%RF plot
subplot(4,2,[1 2])
hold on
ylabel({'Mean Velocity' ; '[cm/s]'})
xlabel('Time [s]')
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

%A trials
sub_orderA = [0 3 5 7];

for ii=2:4
    subplot(4,2,sub_orderA(ii))
    hold on
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
    subplot(4,2,sub_orderB(ii))
    hold on
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

% 
% 16 = -2s
% 46 = -1s
% 76 = 0s
% 106 = 1s
% 136 = 2s


%%
figure
hold on
title('A zone speed')
for aa=1:4
    for ss=2
        for ll=1:size(speed_data{aa}{ss}.lap_idx.A,2)
            %plot line plot along range
            plot(speed_data{aa}{ss}.rewards.A.speed{ll}((speed_data{aa}{ss}.rewards.A.Imin(ll)-idx_width(1)):(speed_data{aa}{ss}.rewards.A.Imin(ll)+idx_width(2))));
            %plot line showing start of reward range
            plot([idx_width(1) idx_width(1)],[0 25],'k')
        end
    end
end




