%% File directories
%I45_RT  (1,1,1,1)
pathDir{1}{1} = 'G:\Figure_1_OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_rand_d1_052218';
pathDir{1}{2} = 'G:\Figure_1_OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_5A5B_053018';
pathDir{1}{3} = 'G:\Figure_1_OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_3A3B_060518';
pathDir{1}{4} = 'G:\Figure_1_OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_AB_061418';

%I46 (1,1,1,1)
pathDir{2}{1} = 'G:\Figure_1_OCGOL_learning_long_term\I46\behavior_only\I46_rand_d1_052918';
pathDir{2}{2} = 'G:\Figure_1_OCGOL_learning_long_term\I46\behavior_only\I46_5A5B_060118';
pathDir{2}{3} = 'G:\Figure_1_OCGOL_learning_long_term\I46\behavior_only\I46_3A3B_060718';
pathDir{2}{4} = 'G:\Figure_1_OCGOL_learning_long_term\I46\behavior_only\I46_AB_061518';

%I47_RS (1,0,1,1)
pathDir{3}{1} = 'G:\Figure_1_OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_rand_d2_051518';
pathDir{3}{2} = 'G:\Figure_1_OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_5AB_d7_052218';
pathDir{3}{3} = 'G:\Figure_1_OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_3AB_d8_052418';
pathDir{3}{4} = 'G:\Figure_1_OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_AB_061418';

%I47_LP (1,1,1,1)
pathDir{4}{1} = 'G:\Figure_1_OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_rand_d2_051518';
pathDir{4}{2} = 'G:\Figure_1_OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_5AB_d1_051718';
pathDir{4}{3} = 'G:\Figure_1_OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_3AB_d8_052418';
pathDir{4}{4} = 'G:\Figure_1_OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_AB_061418';

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

%1 - random reward
%2 - 5A5B
%3 - 3A3B
%4 - random AB

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
ylabel({'Mean Speed' ; '[cm/s]'})
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

% Frame indices corresponding to time relative to reward zone
% 16 = -2s
% 46 = -1s
% 76 = 0s
% 106 = 1s
% 136 = 2s


%% Generate mean velocity 2s prior and 2s post entry to reward zone

%reward/peri-reward zone speeds
%animal index,then session; then lap x reward zone velocity - event split - 76 fr - 0 s
%76  30  = 46
%2s prior range: 17-76   
%2s post range: 77-136

%for each animal, for each animal
for aa=1:4
    for ss=1:4
        pre_speed{aa}{ss}.A = A_speed{aa}{ss}(:,17:76);
        post_speed{aa}{ss}.A = A_speed{aa}{ss}(:,77:136);

        pre_speed{aa}{ss}.B = B_speed{aa}{ss}(:,17:76);
        post_speed{aa}{ss}.B = B_speed{aa}{ss}(:,77:136);
        
        %no first session (RF) in relative zones (1 lap category)
        if ss~=1
            pre_speed{aa}{ss}.ArelB = A_speed_relB{aa}{ss}(:,17:76);
            post_speed{aa}{ss}.ArelB = A_speed_relB{aa}{ss}(:,77:136);
            
            pre_speed{aa}{ss}.BrelA = B_speed_relA{aa}{ss}(:,17:76);
            post_speed{aa}{ss}.BrelA = B_speed_relA{aa}{ss}(:,77:136);
            
        end
        
        %mean on each lap
        pre_post_speed_mean{aa}{ss}.A = [mean(pre_speed{aa}{ss}.A,2),mean(post_speed{aa}{ss}.A,2)];
        
        pre_post_speed_mean{aa}{ss}.B = [mean(pre_speed{aa}{ss}.B,2),mean(post_speed{aa}{ss}.B,2)];
        
        if ss~=1 %no input for RF
            pre_post_speed_mean{aa}{ss}.ArelB = [mean(pre_speed{aa}{ss}.ArelB,2),mean(post_speed{aa}{ss}.ArelB,2)];
            
            pre_post_speed_mean{aa}{ss}.BrelA = [mean(pre_speed{aa}{ss}.BrelA,2),mean(post_speed{aa}{ss}.BrelA,2)];
        end
        
        %get across lap mean and sem for each animal session, pre,post,
        %animal
        pre_post_speed_mean_lap.A(ss,:,aa) = mean(pre_post_speed_mean{aa}{ss}.A,1);
        pre_post_speed_sem_lap.A(ss,:,aa) = std(pre_post_speed_mean{aa}{ss}.A,0,1)./sqrt(size(pre_post_speed_mean{aa}{ss}.A,1));

        pre_post_speed_mean_lap.B(ss,:,aa) = mean(pre_post_speed_mean{aa}{ss}.B,1);
        pre_post_speed_sem_lap.B(ss,:,aa) = std(pre_post_speed_mean{aa}{ss}.B,0,1)./sqrt(size(pre_post_speed_mean{aa}{ss}.B,1));   
        
        if ss~=1 %no input for RF
            pre_post_speed_mean_lap.ArelB(ss,:,aa) = mean(pre_post_speed_mean{aa}{ss}.ArelB,1);
            pre_post_speed_sem_lap.ArelB(ss,:,aa) = std(pre_post_speed_mean{aa}{ss}.ArelB,0,1)./sqrt(size(pre_post_speed_mean{aa}{ss}.ArelB,1));
            
            pre_post_speed_mean_lap.BrelA(ss,:,aa) = mean(pre_post_speed_mean{aa}{ss}.BrelA,1);
            pre_post_speed_sem_lap.BrelA(ss,:,aa) = std(pre_post_speed_mean{aa}{ss}.BrelA,0,1)./sqrt(size(pre_post_speed_mean{aa}{ss}.BrelA,1));
        end
        
    end
end

%% Statistics - paired wilcoxon on animal means - pre vs post for all trial (all animals)

%fewer than 5 pairs - don't use wilcoxon; reasonable to assume normality
%and use simple paired t-test
%for each session type
for ss=1:4
    [h.A(ss),p.A(ss),ci.A(:,ss),stats.A(ss)] = ttest(squeeze(pre_post_speed_mean_lap.A(ss,1,:)), squeeze(pre_post_speed_mean_lap.A(ss,2,:)));
    [h.B(ss),p.B(ss),ci.B(:,ss),stats.B(ss)] = ttest(squeeze(pre_post_speed_mean_lap.B(ss,1,:)), squeeze(pre_post_speed_mean_lap.B(ss,2,:)));
    if ss~=1
        [h.ArelB(ss),p.ArelB(ss),ci.ArelB(:,ss),stats.ArelB(ss)] = ttest(squeeze(pre_post_speed_mean_lap.ArelB(ss,1,:)), squeeze(pre_post_speed_mean_lap.ArelB(ss,2,:)));
        [h.BrelA(ss),p.BrelA(ss),ci.BrelA(:,ss),stats.BrelA(ss)] = ttest(squeeze(pre_post_speed_mean_lap.BrelA(ss,1,:)), squeeze(pre_post_speed_mean_lap.BrelA(ss,2,:)));
    end
end

%get 0,1,2,3 star significance level
%-1 - ns - >= 0.05
%1 * - <0.05, >= 0.01
%2 ** - <0.01 >= 0.001
%3 *** - <0.001
p_sig_level.A = get_star_sig(p.A);
p_sig_level.B = get_star_sig(p.B);
p_sig_level.ArelB = get_star_sig(p.ArelB);
p_sig_level.BrelA = get_star_sig(p.BrelA);

%set RF session (1) to nan in zone sig on opposing laps (since there are
%non)
p_sig_level.ArelB(1) = nan;
p_sig_level.BrelA(1) = nan;

%% Plot zone data

%A zones in A laps
%plot
figure('Position',[769   223   758   684])
subplot(2,2,1)
hold on
title('A laps - A reward zone')
idx_skip= 0;
xticks([1 2 3.5 4.5 6.0 7.0 8.5 9.5]);
xticklabels({'Pre','Post','Pre','Post','Pre','Post','Pre','Post'})
xtickangle(45)
for ss=1:4
    
    
    ylabel('Speed [cm/s]')
    ylim([0  25])
    xlim([0 10.5])
    %pre A rew on A laps
    errorbar([0.85, 0.95, 1.05, 1.15]+idx_skip,squeeze(pre_post_speed_mean_lap.A(ss,1,:)),squeeze(pre_post_speed_sem_lap.A(ss,1,:))'.','LineStyle','none','Color',[65,105,225]./255,'LineWidth',1.5)
    %post A rew on A laps
    errorbar([0.85, 0.95, 1.05, 1.15]+1+idx_skip,squeeze(pre_post_speed_mean_lap.A(ss,2,:)),squeeze(pre_post_speed_sem_lap.A(ss,2,:))'.','LineStyle','none','Color',[65,105,225]./255,'LineWidth',1.5)
    %connecting line
    plot([0.85, 0.95, 1.05, 1.15; [0.85, 0.95, 1.05, 1.15]+1]+idx_skip,[squeeze(pre_post_speed_mean_lap.A(ss,1,:)), squeeze(pre_post_speed_mean_lap.A(ss,2,:))]','Color',[65,105,225]./255,'LineWidth',1.5)
    
    idx_skip=  idx_skip+2.5;
end
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%B zones in B laps
subplot(2,2,2)
hold on
title('B laps - B reward zone')
idx_skip= 0;
xticks([1 2 3.5 4.5 6.0 7.0 8.5 9.5]);
xticklabels({'Pre','Post','Pre','Post','Pre','Post','Pre','Post'})
xtickangle(45)
for ss=1:4
    
    
    ylabel('Speed [cm/s]')
    ylim([0  25])
    xlim([0 10.5])
    %pre A rew on A laps
    errorbar([0.85, 0.95, 1.05, 1.15]+idx_skip,squeeze(pre_post_speed_mean_lap.B(ss,1,:)),squeeze(pre_post_speed_sem_lap.B(ss,1,:))'.','LineStyle','none','Color',[220,20,60]./255,'LineWidth',1.5)
    %post A rew on A laps
    errorbar([0.85, 0.95, 1.05, 1.15]+1+idx_skip,squeeze(pre_post_speed_mean_lap.B(ss,2,:)),squeeze(pre_post_speed_sem_lap.B(ss,2,:))'.','LineStyle','none','Color',[220,20,60]./255,'LineWidth',1.5)
    %connecting line
    plot([0.85, 0.95, 1.05, 1.15; [0.85, 0.95, 1.05, 1.15]+1]+idx_skip,[squeeze(pre_post_speed_mean_lap.B(ss,1,:)), squeeze(pre_post_speed_mean_lap.B(ss,2,:))]','Color',[220,20,60]./255,'LineWidth',1.5)
    
    idx_skip=  idx_skip+2.5;
end
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%A zones in B laps
subplot(2,2,3)
hold on
title('B laps - A reward zone')
idx_skip= 2.5;
xticks([3.5 4.5 6.0 7.0 8.5 9.5]);
xticklabels({'Pre','Post','Pre','Post','Pre','Post','Pre','Post'})
xtickangle(45)
for ss=2:4
    
    ylabel('Speed [cm/s]')
    ylim([0  25])
    xlim([0 10.5])
    %pre A rew on A laps
    errorbar([0.85, 0.95, 1.05, 1.15]+idx_skip,squeeze(pre_post_speed_mean_lap.ArelB(ss,1,:)),squeeze(pre_post_speed_sem_lap.ArelB(ss,1,:))'.','LineStyle','none','Color',[65,105,225]./255,'LineWidth',1.5)
    %post A rew on A laps
    errorbar([0.85, 0.95, 1.05, 1.15]+1+idx_skip,squeeze(pre_post_speed_mean_lap.ArelB(ss,2,:)),squeeze(pre_post_speed_sem_lap.ArelB(ss,2,:))'.','LineStyle','none','Color',[65,105,225]./255,'LineWidth',1.5)
    %connecting line
    plot([0.85, 0.95, 1.05, 1.15; [0.85, 0.95, 1.05, 1.15]+1]+idx_skip,[squeeze(pre_post_speed_mean_lap.ArelB(ss,1,:)), squeeze(pre_post_speed_mean_lap.ArelB(ss,2,:))]','Color',[65,105,225]./255,'LineWidth',1.5)
    
    idx_skip=  idx_skip+2.5;
end
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%B zones in A laps
subplot(2,2,4)
hold on
title('A laps - B reward zone')
idx_skip= 2.5;
xticks([3.5 4.5 6.0 7.0 8.5 9.5]);
xticklabels({'Pre','Post','Pre','Post','Pre','Post','Pre','Post'})
xtickangle(45)
for ss=2:4
    
    ylabel('Speed [cm/s]')
    ylim([0  25])
    xlim([0 10.5])
    %pre A rew on A laps
    errorbar([0.85, 0.95, 1.05, 1.15]+idx_skip,squeeze(pre_post_speed_mean_lap.BrelA(ss,1,:)),squeeze(pre_post_speed_sem_lap.BrelA(ss,1,:))'.','LineStyle','none','Color',[220,20,60]./255,'LineWidth',1.5)
    %post A rew on A laps
    errorbar([0.85, 0.95, 1.05, 1.15]+1+idx_skip,squeeze(pre_post_speed_mean_lap.BrelA(ss,2,:)),squeeze(pre_post_speed_sem_lap.BrelA(ss,2,:))'.','LineStyle','none','Color',[220,20,60]./255,'LineWidth',1.5)
    %connecting line
    plot([0.85, 0.95, 1.05, 1.15; [0.85, 0.95, 1.05, 1.15]+1]+idx_skip,[squeeze(pre_post_speed_mean_lap.BrelA(ss,1,:)), squeeze(pre_post_speed_mean_lap.BrelA(ss,2,:))]','Color',[220,20,60]./255,'LineWidth',1.5)
    
    idx_skip=  idx_skip+2.5;
end
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)



%% 


%% Generate supplementary plots with slope change -1s to 1s in reward zone (try -1 to 0s)

%1s prior range: 47-76   
%1s post range: 77-106
%1 - RF; %2 - 5A5B; %3 - 3A3B; %4 - rand AB

%reward/peri-reward zone speeds
%animal index,then session; then lap x reward zone velocity - event split - 76 fr - 0 s
%76  30  = 46
%2s prior range: 17-76   
%2s post range: 77-136

%for each animal, for each animal
for aa=1:4
    for ss=1:4
        speed_sl{aa}{ss}.A = A_speed{aa}{ss}(:,47:76);
        %post_speed_sl{aa}{ss}.A = A_speed{aa}{ss}(:,77:136);

        speed_sl{aa}{ss}.B = B_speed{aa}{ss}(:,47:76);
        %post_speed_sl{aa}{ss}.B = B_speed{aa}{ss}(:,77:136);
        
        %no first session (RF) in relative zones (1 lap category)
        if ss~=1
            speed_sl{aa}{ss}.ArelB = A_speed_relB{aa}{ss}(:,47:76);
            %speed_sl{aa}{ss}.ArelB = A_speed_relB{aa}{ss}(:,77:136);
            
            speed_sl{aa}{ss}.BrelA = B_speed_relA{aa}{ss}(:,47:76);
            %post_speed_sl{aa}{ss}.BrelA = B_speed_relA{aa}{ss}(:,77:136);
            
        end
    end
end

% Get slopes for each lap with interval
for aa=1:4
    for ss=1:4
        %for each lap
        for ll=1:size(speed_sl{aa}{ss}.A,1)
            slopes{aa}{ss}.A(ll,:) =  polyfit(linspace(-1,0,30),speed_sl{aa}{ss}.A(ll,:),1);
        end
        
        for ll=1:size(speed_sl{aa}{ss}.B,1)
            slopes{aa}{ss}.B(ll,:) =  polyfit(linspace(-1,0,30),speed_sl{aa}{ss}.B(ll,:),1);
        end
        
    end
    
    for ss=2:4
        %for each lap
        for ll=1:size(speed_sl{aa}{ss}.ArelB,1)
            slopes{aa}{ss}.ArelB(ll,:) =  polyfit(linspace(-1,0,30),speed_sl{aa}{ss}.ArelB(ll,:),1);
        end
        
        for ll=1:size(speed_sl{aa}{ss}.BrelA,1)
            slopes{aa}{ss}.BrelA(ll,:) =  polyfit(linspace(-1,0,30),speed_sl{aa}{ss}.BrelA(ll,:),1);
        end
        
    end
end

%for each animal and each session, get mean slope
for aa=1:4
    for ss=1:4
        mean_slope.A(ss,aa) = mean(slopes{aa}{ss}.A(:,1));
        mean_slope.B(ss,aa) = mean(slopes{aa}{ss}.B(:,1));
    end
    
    for ss=2:4
        mean_slope.ArelB(ss,aa) = mean(slopes{aa}{ss}.ArelB(:,1));
        mean_slope.BrelA(ss,aa) = mean(slopes{aa}{ss}.BrelA(:,1));
    end
end



%         %mean on each lap
          %speed_sl_mean{aa}{ss}.A = mean(pre_speed{aa}{ss}.A,2),mean(post_speed{aa}{ss}.A,2)];
%         
%         pre_post_speed_mean{aa}{ss}.B = [mean(pre_speed{aa}{ss}.B,2),mean(post_speed{aa}{ss}.B,2)];
%         
%         if ss~=1 %no input for RF
%             pre_post_speed_mean{aa}{ss}.ArelB = [mean(pre_speed{aa}{ss}.ArelB,2),mean(post_speed{aa}{ss}.ArelB,2)];
%             
%             pre_post_speed_mean{aa}{ss}.BrelA = [mean(pre_speed{aa}{ss}.BrelA,2),mean(post_speed{aa}{ss}.BrelA,2)];
%         end
%         
%         %get across lap mean and sem for each animal session, pre,post,
%         %animal
%         pre_post_speed_mean_lap.A(ss,:,aa) = mean(pre_post_speed_mean{aa}{ss}.A,1);
%         pre_post_speed_sem_lap.A(ss,:,aa) = std(pre_post_speed_mean{aa}{ss}.A,0,1)./sqrt(size(pre_post_speed_mean{aa}{ss}.A,1));
% 
%         pre_post_speed_mean_lap.B(ss,:,aa) = mean(pre_post_speed_mean{aa}{ss}.B,1);
%         pre_post_speed_sem_lap.B(ss,:,aa) = std(pre_post_speed_mean{aa}{ss}.B,0,1)./sqrt(size(pre_post_speed_mean{aa}{ss}.B,1));   
%         
%         if ss~=1 %no input for RF
%             pre_post_speed_mean_lap.ArelB(ss,:,aa) = mean(pre_post_speed_mean{aa}{ss}.ArelB,1);
%             pre_post_speed_sem_lap.ArelB(ss,:,aa) = std(pre_post_speed_mean{aa}{ss}.ArelB,0,1)./sqrt(size(pre_post_speed_mean{aa}{ss}.ArelB,1));
%             
%             pre_post_speed_mean_lap.BrelA(ss,:,aa) = mean(pre_post_speed_mean{aa}{ss}.BrelA,1);
%             pre_post_speed_sem_lap.BrelA(ss,:,aa) = std(pre_post_speed_mean{aa}{ss}.BrelA,0,1)./sqrt(size(pre_post_speed_mean{aa}{ss}.BrelA,1));
%         end
        
  %  end
%end

%% Plot individual traces by lap
%30 frames/s 

%-1 to 1s
idx_width_custom = [30,30];

figure
hold on
title('A zone speed')
for aa=2%1:4
    for ss=4
        %for A laps
        for ll=1:size(speed_data{aa}{ss}.lap_idx.A,2)
            %plot line plot along range
            plot(speed_data{aa}{ss}.rewards.A.speed{ll}((speed_data{aa}{ss}.rewards.A.Imin(ll)-idx_width_custom(1)):(speed_data{aa}{ss}.rewards.A.Imin(ll)+idx_width_custom(2))));
            %plot line showing start of reward range
            plot([idx_width_custom(1) idx_width_custom(1)],[0 25],'k')
        end
    end
end




