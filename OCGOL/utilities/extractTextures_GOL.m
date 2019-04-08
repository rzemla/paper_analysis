function [Behavior] = extractTextures_GOL(CSV, Behavior,options)

%% Assign variables from Behavior struct

position = Behavior.position;
position_norm = Behavior.normalizedposition;
time = Behavior.time;
lap = Behavior.lap;

%% find idx's of start of first complete lap and end of last complete lap
start_lap_idx = find(time == lap{1}(1));
end_lap_idx = find(time == lap{end}(2));

%% Channel and voltage ranges for each microtexture - remove eventually

%RF day 0 for I55_RTLS
% textureRanges =  [1.20 1.26;...    % Cue 2
%     1.80 1.86;...    % Cue 3
%     2.40 2.47;];     % Cue 4

%channel which corresponds to texture voltage marks
tagLocationCh = CSV(:,3);

%% Texture signal discovery

%how many voltage bins; 20 = 0.25 V, 10 = 0.50 V
voltageBins = 10;

%discover voltage peaks
%voltage 0.4 --> 0.8 V (to avoid discovery of double tagged tex cue 0.6 mixed with 1.2V signal)
[pks, pks_idx] = findpeaks(tagLocationCh, 'MinPeakHeight',0.8);

%exclude voltage peaks based on positional proximity
peak_diff_pos = [position(pks_idx), [nan ;diff(position(pks_idx))]];

%find voltage signals with positional diff < 0.3 cm
exclude_peaks_idx = find(peak_diff_pos(:,2) < 0.3 & peak_diff_pos(:,2) >= 0.0);
keep_peaks = false(size(peak_diff_pos,1),1);
keep_peaks(exclude_peaks_idx) = 1;
keep_peaks = ~keep_peaks;

%final peaks to include
pks_final_idx = pks_idx(keep_peaks);
pks_final = pks(keep_peaks);


%bin the voltage signals associated with textures
[N_voltage_bins,~] = histcounts(pks_final,voltageBins);

%find unique bins (= # of textures for input into k means)
uniqueTexBins = length(find(N_voltage_bins ~= 0));

%k-means cluster - find each texture and voltage centroid
%C - how many voltage classes
[idx_C,C] = kmeans(pks_final,uniqueTexBins);

%find the index of each cluster class (without sorting)
for kk=1:size(C,1)
    %get indices for each class of signals
    class_C_idx{kk} = find(idx_C == kk);
    %get position, norm_position, time of each class of signal
    %index from behavioral recording
    pks_final_indexes{kk} = pks_final_idx(class_C_idx{kk});
    %position
    pks_final_pos{kk} = position(pks_final_idx(class_C_idx{kk}));
    %normalized position
    pks_final_pos_norm{kk} = position_norm(pks_final_idx(class_C_idx{kk}));
    %time
    pks_final_time{kk} = time(pks_final_idx(class_C_idx{kk}));
end

%% Check signals that are detected

x=1;
figure;
hold on


%% Parse the signals into their respective classes

%split the reward tags (alternate between two reward zones
%split the first cue tags and sound on tags

%preallocate empty index assignment
lap_cue_tex_idx = [];
reward_cue_tex_idx = [];

for kk=1:size(pks_final_pos,2)
    %first cue tag vs. sound tag
    first_cue_idx = find(pks_final_pos{kk} > 17 & pks_final_pos{kk} < 26);
    %when discovered:
    if ~isempty(first_cue_idx)
        %within this texture set, which idxs correspond to first cue
        %sound_first_cue_split_idx = first_cue_idx;
        %make a logical with first cue signals = 1
        first_cue_on = false(size(pks_final_pos{kk},1),1);
        first_cue_on(first_cue_idx) = 1;
        
        %mark which set of textures corresponds to first cue vs audio on cue
        first_cue_tex_idx = kk;
    end
    
    %check for presence of early reward
    reward_cue_idx = find(pks_final_pos{kk} > 52 & pks_final_pos{kk} < 62);
    %when discovered:
    if ~isempty(reward_cue_idx)
        %within this cue, get the indices of the early rewards
        %reward_early_idx = reward_cue_idx;
        %make a logical with first cue signals = 1
        reward_early_on = false(size(pks_final_pos{kk},1),1);
        reward_early_on(reward_cue_idx) = 1;
        
        reward_cue_tex_idx = kk;
    end
    
    %lap cue idx (redundant signal that animal crossed the lap)
    lap_cue_idx = find(pks_final_pos{kk} > 0 & pks_final_pos{kk} < 5);
    %when discovered:
    if ~isempty(lap_cue_idx)
        %which set of signals corresponds to lap idx
        lap_cue_tex_idx = kk;
    end
    
end

%define remaining idxs of texture signals
tex_signals_idx = 1:size(pks_final_pos,2);

%check if reward lap flag is empty
if ~isempty(lap_cue_tex_idx)
tex_signals_idx([first_cue_tex_idx, reward_cue_tex_idx, lap_cue_tex_idx]) = [];
%if no lap texture and there is a reward texture
elseif isempty(lap_cue_tex_idx) && ~isempty(reward_cue_tex_idx)
    tex_signals_idx([first_cue_tex_idx, reward_cue_tex_idx]) = [];
end


%isolate first cue and sound signal in separate cells by absolute indices
first_cue_onsets = pks_final_indexes{first_cue_tex_idx}(first_cue_on);

sound_onsets = pks_final_indexes{first_cue_tex_idx}(~first_cue_on);

early_reward_onsets = pks_final_indexes{reward_cue_tex_idx}(reward_early_on);

late_reward_onsets = pks_final_indexes{reward_cue_tex_idx}(~reward_early_on);

%% Get the indices of each class of signal

%place all texture indices on behavior into cells
%first tex in  first cell
textures_idx{1} = first_cue_onsets;

for ii=2:(size(tex_signals_idx,2)+1)
    textures_idx{ii} =  pks_final_indexes{tex_signals_idx(ii-1)};
end

%rewards
%early (B trial = 3)
rewards_idx{1} = early_reward_onsets;
%late (A trial = 2)
rewards_idx{2} = late_reward_onsets;

%sound
sound_idx = sound_onsets;

%sound, 4 microtextures, 2 reward zones
%% Get position, normalized postion, and time for each signal + store indices

%microtextures
for ii=1:size(textures_idx,2)
    textures{ii}.position = position(textures_idx{ii});
    textures{ii}.position_norm = position_norm(textures_idx{ii});
    textures{ii}.time = time(textures_idx{ii});
    textures{ii}.idx = textures_idx{ii};
end

%rewards
for ii=1:size(rewards_idx,2)
    rewards{ii}.position = position(rewards_idx{ii});
    rewards{ii}.position_norm = position_norm(rewards_idx{ii});
    rewards{ii}.time = time(rewards_idx{ii});
    rewards{ii}.idx = rewards_idx{ii};    
end

%make copy of reward_idx
rewards_full_lap_idx = rewards_idx;

for ii =1:size(rewards_full_lap_idx,2)
    rewards_full_lap_idx{ii}((rewards_full_lap_idx{ii} < start_lap_idx)) = [];
    rewards_full_lap_idx{ii}((rewards_full_lap_idx{ii} > end_lap_idx)) = [];
end

%complete lap rewards
for ii=1:size(rewards_idx,2)
    rewards{ii}.full_laps.position = position(rewards_full_lap_idx{ii});
    rewards{ii}.full_laps.position_norm = position_norm(rewards_full_lap_idx{ii});
    rewards{ii}.full_laps.time = time(rewards_full_lap_idx{ii});
    rewards{ii}.full_laps.idx = rewards_full_lap_idx{ii};    
end

%sound
    sound.position = position(sound_idx);
    sound.position_norm = position_norm(sound_idx);
    sound.time = time(sound_idx);
    sound.idx = sound_idx;

%% (Not significant) - check how the distance/time changes for each cue if shift to the start of signal    
    
%% Determine voltages ranges for each class of signal
%sort k-means cluster centroids in ascending voltage
C = sort(C);

%assign voltage range
voltRange = 0.05;

%assign discovered texture ranges
%preallocate
textureRanges = zeros(size(C,1),2);

for ii = 1:size(C,1)
    %low voltage thres
    textureRanges(ii,1) = C(ii) - voltRange;
    %high voltage thres
    textureRanges(ii,2) = C(ii) + voltRange;
end

%% Plot each class of textures, rewards and sound

figure;
hold on
xlabel('Time [s]');
ylabel('Position [cm]');
plot(time,position, 'k')

%A trial rewards (blue)
stem(rewards{1}.time, rewards{1}.position, 'b','LineStyle','none')
%B trial rewards (red)
stem(rewards{2}.time, rewards{2}.position, 'r','LineStyle','none')

%colormap to iterate through
colorRange = cbrewer('qual', 'Dark2',8);

%textures
for ii=1:size(textures,2)
stem(textures{ii}.time, textures{ii}.position,'Color',colorRange(ii+2,:),'LineStyle','none');
end

%sound
stem(sound.time, sound.position, 'g','LineStyle','none')

%% Trial type discovery - early or late block run for GOL - work on this

%assign trial type channel to separate variable
trialTypeCh = CSV(:,2);

%run discovery of trial types
[trialRanges] = discoverTrialType(trialTypeCh, position);

%get time and position of trial textures (may change to findpeaks approach
%like done above depending on stability)
[trialType] = defineTrialSignal(trialRanges, trialTypeCh, Behavior);

%% Lick signal onsets, position, norm_position, time, idx

%place signal in separate channel
lickCh = CSV(:,8); 
%digitize signal     
lick_high_idx = find(lickCh > 3.5);
lick_digital = zeros(size(lickCh,1),1);
lick_digital(lick_high_idx) = 1;

%find lick onsets
lick_diff = diff(lick_digital);
lick_on = find(lick_diff == 1);
lick_on = lick_on + 1;

%place into struct with postion, norm_position, 
lick.position = position(lick_on);
lick.position_norm = position_norm(lick_on);
lick.time = time(lick_on);
lick.idx = lick_on;

% %plot - QC - make sure indices correspond to lick onsets
% figure;
% hold on;
% ylim([0 2]);
% plot(lick_digital)
% stem(lick_on, 2*ones(1,size(lick_on,1)), 'r');
    
%% Reward collected signal onsets, position, norm_position, time, idx

%place signal in separate channel
rewardCollectedCh = CSV(:,7); 

%digitize signal
reward_high_idx = find(rewardCollectedCh > 4.5);
reward_digital = zeros(size(rewardCollectedCh,1),1);
reward_digital(reward_high_idx) = 1;

%find lick onsets
reward_diff = diff(reward_digital);
reward_on = find(reward_diff == 1);
reward_on = reward_on + 1;

%place into struct with postion, norm_position, 
reward_coll.position = position(reward_on);
reward_coll.position_norm = position_norm(reward_on);
reward_coll.time = time(reward_on);
reward_coll.idx = reward_on;

% %plot - QC - make sure indices correspond to lick onsets
% figure;
% hold on;
% ylim([0 2]);
% plot(reward_digital)
% stem(reward_on, 2*ones(1,size(reward_on,1)), 'r');

%% OCGOL quality control
if strcmpi(options.BehaviorType,'OCGOL')
    
    %runs several QC checks to make sure that important signals were not missed
    %displays alerts when signals are missed or don't match to what is
    %expected
    [lap_id] = OCGOL_QC(Behavior,textures, rewards, trialType, sound);
    
elseif strcmpi(options.BehaviorType,'GOL')
    
    [lap_id] = GOL_QC(Behavior,textures, rewards, trialType);
end

%% Plot summary plot

figure;
%all textures and rewards
subplot(4,1,1)
hold on
title('Textures, reward zones, and sound');
xlabel('Time [s]');
ylabel('Position [cm]');
plot(time,position, 'k')

%B trial rewards (red) - early
stem(rewards{1}.time, rewards{1}.position, 'r','LineStyle','none')
%A trial rewards (blue) - far
stem(rewards{2}.time, rewards{2}.position, 'b','LineStyle','none')

%colormap to iterate through
colorRange = cbrewer('qual', 'Dark2',8);
%textures
for ii=1:size(textures,2)
stem(textures{ii}.time, textures{ii}.position,'Color',colorRange(ii+2,:),'LineStyle','none');
end
%sound
stem(sound.time, sound.position, 'g','LineStyle','none')

%reward zone and trial type
subplot(4,1,2)
hold on
title('Reward zones and trial types');
xlabel('Time [s]');
ylabel('Position [cm]');
plot(time,position, 'k')

%reward zone signal
%B trial rewards (red) - early
stem(rewards{1}.time, rewards{1}.position, 'r','LineStyle','none')
%A trial rewards (blue) - late
stem(rewards{2}.time, rewards{2}.position, 'b','LineStyle','none')

%trial type signal
%A trial rewards (blue) (2 - high voltage, far, A, blue) - check name
%assignment
stem(trialType.time{2}, trialType.position{2}, 'b','LineStyle','none')
%A trial rewards (blue) (1 - high voltage, far, A, blue)
stem(trialType.time{1}, trialType.position{1}, 'r','LineStyle','none')

%correct reward zones and licks
subplot(4,1,3)
hold on
title('Reward zones licks (green)');
xlabel('Time [s]');
ylabel('Position [cm]');
plot(time,position, 'k')

%licks
stem(lick.time, lick.position, 'g','LineStyle','none')

%reward zone signal
%B trial rewards (red) - early
stem(rewards{1}.time, rewards{1}.position, 'r','LineStyle','-')
%A trial rewards (blue) - late
stem(rewards{2}.time, rewards{2}.position, 'b','LineStyle','-')

%licks vs. reward collected
subplot(4,1,4)
hold on
title('Reward zones licks (green) vs reward collected (red)');
xlabel('Time [s]');
ylabel('Position [cm]');
plot(time,position, 'k')

%licks
stem(lick.time, lick.position, 'g','LineStyle','none')

%reward collected
stem(reward_coll.time, reward_coll.position, 'r','LineStyle','--')
 
 %% Save to textures and reward substruct
 
%sound
Behavior.sound = sound;
 
%texture signals
Behavior.textures = textures;
 
 %reward signals
 Behavior.rewards = rewards;
 
 %lick signals
 Behavior.lick = lick; 
 
 %reward collected signal
 Behavior.reward_coll = reward_coll;
 
 %laps identified based on reward and texture signal
 Behavior.lap_id = lap_id;
 
 %trial type signals
 Behavior.trialType = trialType;
  
 %laps with respective A or B rewards/trials
 %B trials (3); early;
 Behavior.reward.reward_early = find(reward_early_on);
 %A trials; (2) late reward indices
 Behavior.reward.reward_late =  find(~reward_early_on);
 
 %% Testing code that is related
 
 
%check all tags
% figure
% for ii=1:size(pks_final_idx,1)
%     hold on
%     title(num2str(ii));
%     plot(tagLocationCh(pks_final_idx(ii)-10: pks_final_idx(ii)+10));
%     pause;
%      clf;
% end

% ii=5;
% 
%     figure;
%     hold on
%     title(num2str(ii));
%     plot(tagLocationCh(pks_final_idx(ii)-10: pks_final_idx(ii)+10));
%     
%     figure;
%     [double_pk_temp,~] = findpeaks(diff(tagLocationCh(pks_final_idx(ii)-10: pks_final_idx(ii)+10)),'MinPeakHeight',0.3);
%     
% %deal with 
% %if
% for ee=1:size(sig_2_discharge,2)
%     ii= sig_2_discharge(ee);
%     %temp{ee} = find(diff(tagLocationCh(pks_final_idx(ii)-10: pks_final_idx(ii)+10)) > 0.3);
%     [double_pk_temp,~] = findpeaks(diff(tagLocationCh(pks_final_idx(ii)-10: pks_final_idx(ii)+10)),'MinPeakHeight',0.3);
%     temp(ee) = length(double_pk_temp);
% end
% 
% for ee=1:size(all_rest,2)
%     ii= all_rest(ee);
%     [double_pk_temp_rest,~] = findpeaks(diff(tagLocationCh(pks_final_idx(ii)-10: pks_final_idx(ii)+10)),'MinPeakHeight',0.3);
%     temp_rest(ee) = length(double_pk_temp_rest);
% end
% 
% sig_2_discharge = [4 11 18 25 32 39 46 53 60 67 74 81 88 95 102 109 116 123 130 137 144 151 158];
% 
% all_rest = 1:158;
% all_rest(sig_2_discharge) = [];


%check the each voltage class has the same position; if divergent position,
%separate into additional classes
 
 
% %% Plot peak filtered
% figure
% subplot(2,1,1)
% hold on;
% title('Before peak filter');
% plot(time,position, 'k')
% stem(time(pks_idx),position(pks_idx),'r')
% hold off
% 
% subplot(2,1,2)
% hold on;
% title('After peak filter');
% plot(time,position, 'k')
% stem(time(pks_final_idx),position(pks_final_idx),'g')
% hold off


%% Plot selected textures with signal clustering

% figure
% hold on;
% title('Before peak filter');
% plot(time,position, 'k')
% for tt=1:size(pks_final_pos,2)
%     stem(pks_final_time{tt},pks_final_pos{tt},'Color',colorRange(tt,:))
% end
% 
% hold off

 %% Plot the textures on the restricted laps
        
%  figure;
%  hold on
%  title('Texture locations on each lap');
%  %plot time (seconds) against normalized position
%  plot(Behavior.time,Behavior.normalizedposition,'k');
%  xlabel('Time [s]');
%  ylabel('Position');
%  
%  %plot RFID crossing indices against time
%  plot(Behavior.time,Behavior.lap_binary, 'm');
%  
%  %plot each microtexture (should be 4)
%  for ii=1:size(time_tex,2)
%      stem(time_tex{ii}, position_norm_tex{ii}, '*r');
%  end
%  
%  %plot the lap index - lap # scaled down by a factor of 20 - each
%  %step 0.05
%  plot(Behavior.time,Behavior.lapNb./20,'b');

end

