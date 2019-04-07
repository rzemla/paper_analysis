function [Behavior] = extractTextures(CSV, Behavior,options)

%% Assign variables from Behavior struct

position = Behavior.position;
position_norm = Behavior.normalizedposition;
time = Behavior.time;
lap = Behavior.lap;

%% Channel and voltage ranges for each microtexture
%I42_
% textureRanges = [0.6, 0.65;...  % Cue 1 (odor stop tag)
%     1.20 1.26;...    % Cue 2
%     1.80 1.86;...    % Cue 3
%     2.40 2.47;];     % Cue 4

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

%split the reward tags (alternate between two reward zones
%split the first cue tags and sound on tags

for kk=1:size(pks_final_pos,2)
    %first cue tag vs. sound tag
    first_cue_idx = find(pks_final_pos{kk} > 17 & pks_final_pos{kk} < 24);
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

%remaining idxs of texture signals
tex_signals_idx = 1:size(pks_final_pos,2);
tex_signals_idx([first_cue_tex_idx, reward_cue_tex_idx, lap_cue_tex_idx]) = [];


%isolate first cue and sound signal in separate cells by absolute indices
first_cue_onsets = pks_final_indexes{first_cue_tex_idx}(first_cue_on);

sound_onsets = pks_final_indexes{first_cue_tex_idx}(~first_cue_on);

early_reward_onsets = pks_final_indexes{reward_cue_tex_idx}(reward_early_on);

late_reward_onsets = pks_final_indexes{reward_cue_tex_idx}(~reward_early_on);

%textures
%te


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

%% Plot peak filtered
figure
subplot(2,1,1)
hold on;
title('Before peak filter');
plot(time,position, 'k')
stem(time(pks_idx),position(pks_idx),'r')
hold off

subplot(2,1,2)
hold on;
title('After peak filter');
plot(time,position, 'k')
stem(time(pks_final_idx),position(pks_final_idx),'g')
hold off

%% Plot selected textures with signal clustering

%colormap to iterate through
colorRange = cbrewer('qual', 'Dark2',8);

figure
hold on;
title('Before peak filter');
plot(time,position, 'k')
for tt=1:size(pks_final_pos,2)
    stem(pks_final_time{tt},pks_final_pos{tt},'Color',colorRange(tt,:))
end

hold off

%% Trial type discovery

%assign trial type channel to separate variable
trialTypeCh = CSV(:,2);

%run discovery of trial types
[trialRanges] = discoverTrialType(trialTypeCh, position);

%get time and position of trial textures
[trialType] = defineTrialSignal(trialRanges, trialTypeCh, Behavior);

%% Reward collected signal - work on this


            
%% Assign laps, lap index, and trial type to corresponding frames (time and position),          

%find each microtexture within the assigned voltage range 
for ii= 1:size(textureRanges,1)
    texOnIdx{ii} = find(tagLocationCh >= textureRanges(ii,1) & tagLocationCh <= textureRanges(ii,2)); %OK
    texOnsetOnlyIdx{ii} = find(diff(texOnIdx{ii}) >  1) + 1; %OK
    filterTexIndices{ii} = [1; texOnsetOnlyIdx{ii}];
    rawTexIndices{ii} =  texOnIdx{ii}(filterTexIndices{ii});
end

%scan the previous 2 indices and forward 2 indices to see that voltage
%does not exceed or fall below the range for that voltage - remove indices
%that have any values that exceed that range

removeTexIdx = cell(1,size(rawTexIndices,2));
%for each microtexture
for ii=1:size(rawTexIndices,2)
    for jj=1:size(rawTexIndices{ii},1)
        %get neighboring voltages for the each index associated with the
        %texture
        neighborVoltages{ii}{jj} = tagLocationCh((rawTexIndices{ii}(jj)-2):(rawTexIndices{ii}(jj)+2));
        %look at previous 2 voltage samples - should be less than 0.3V
        if sum((neighborVoltages{ii}{jj}(1:2) < textureRanges(ii,1))) ~=2
            removeTexIdx{ii} =  [removeTexIdx{ii}, jj];
        end
        %look at next 1 voltage samples should not be out of range
        if ((neighborVoltages{ii}{jj}(end-1) < textureRanges(ii,1)) || (neighborVoltages{ii}{jj}(end-1) > textureRanges(ii,2)))
            %add index for removal
            removeTexIdx{ii} =  [removeTexIdx{ii}, jj];
        end

    end
end

%% Plot the raw from the voltage threshold voltage threshold (possible duplicates not removed)
% 
% 
% %plot lap and first cue
% figure
% hold on;
% plot(time,position,'k')
% for tt = 1:size(rawTexIndices,2)
%     stem(time(rawTexIndices{tt}), position(rawTexIndices{tt}), 'r')
% end
% 
% figure;
% hold on
% subplot(2,1,1)
% plot(time,position, 'k')
% subplot(2,1,2)
% plot(time,tagLocationCh,'r');
% 
% figure
% hold on;
% plot(time,position, 'k')
% %stem(time(pks_idx),(pks./5)*196,'r')
% stem(time(pks_idx),position(pks_idx),'r')
% 
% figure;
% plot(diff(position(pks_idx)))

%% Clear additional signals

%copy the texture indices
rawTexIdxClean = rawTexIndices;

%remove false positive textures
for ii=1:size(rawTexIndices,2)
    %if found false positve textures for those indices
    if isempty(removeTexIdx{ii}) == 0
        rawTexIdxClean{ii}(removeTexIdx{ii}) = [];
    end
end

%% Find time and position of each texture
 
 %get norm position and offset position
 %for each texture
 for ii = 1:size(rawTexIndices,2)
    position_norm_tex{ii} = position_norm(rawTexIdxClean{ii});
    position_tex{ii} = position(rawTexIdxClean{ii});
    %median of position
    position_med(ii) = median(position(rawTexIdxClean{ii}));
    time_tex{ii} = time(rawTexIdxClean{ii});
 end
 
 %% Exclude textures based on position proximity that do not get filtered in voltage section above

%sort in ascending order and get sorting index
[sort_pos, sort_idx] = sort(position_med,'ascend');

%if difference less than 4cm
dup_idx = find(diff(sort_pos) <= 4);

%which indices to remove
removeIdx = [];

%for each duplicate position
for ii = 1:length(dup_idx)
     %remove the one with higher amount of reads (higher voltage signal)
     size1st = size(position_norm_tex{sort_idx(dup_idx(ii))},1);
     size2nd = size(position_norm_tex{sort_idx(dup_idx(ii)+1)},1);
        
     if size1st > size2nd
         keep_idx = dup_idx(ii);
     else
         keep_idx = dup_idx(ii)+1;
     end
     
     removeIdx = [removeIdx, keep_idx]; 
end

%remove false positive duplicates (from high voltage reads)
position_med(sort_idx(removeIdx)) = [];
position_norm_tex(sort_idx(removeIdx)) = [];
position_tex(sort_idx(removeIdx)) = [];
time_tex(sort_idx(removeIdx)) = [];


 %% Check if any of the textures RFIDs is reward related (after exclusion)

%check the normalized position of each RFID tag relative to middle of track
rewardIdx =[];

 for ii=1:size(position_norm_tex,2)
    %get indices of each and late 
    earlySignal = find(position_norm_tex{ii} < 0.5);
    lateSignal = find(position_norm_tex{ii} > 0.5);
    
    %if both are not empty
    if ~isempty(earlySignal) && ~isempty(lateSignal)
        %index of texture group associated with reward
        rewardIdx = [rewardIdx ii];
        reward_early = earlySignal;
        reward_late = lateSignal;
    end
    
 end
 
 
 %% Quality control check for OCGOL
 %if behavior is OCGOL
 if strcmpi(options.BehaviorType,'OCGOL')
     
     %blue = Odor A trials, late reward (high voltage trial signal)
     %red = Odor B trials, early reward (low voltage trials signal)
     
     %% # of complete lap should equal number of total reward zones
     %total # of complete laps
     nb_complete_laps = size(Behavior.lap,2);
     % # of total reward zones (should equal to number of complete laps)
     nb_reward_zones = sum([size(reward_early,1),size(reward_late,1)]);
     
     %check if number of reward zones (total) equals number of complete laps 
     if (nb_complete_laps == nb_reward_zones)
         disp('Number of complete laps equals number of reward zones');
     else
         disp('Number of complete laps does NOT equal number of reward zones !!!');
     end
     
     %% first cue should occur on each lap
     %check that first cue (between 17 - 24 cm) discharged on every lap
     %responsible for turning off odor/reseting reward allowance/marking
     %first microtexture
     %find texture signal corresponding to this cue
     first_cue_idx = find(position_med > 17 & position_med < 24);
     
     %position of first texture
     first_cue_pos = position_tex{first_cue_idx};
     %time of first texture
     first_cue_time = time_tex{first_cue_idx};
     
     %preallocated list of whether cue was registered in that lap
     cue_registered = false(1,nb_complete_laps);
     
     %take the time of each cue and see if it occured inside of 1 of the
     %laps
     for cc=1:size(first_cue_time,1)
         first_cue_time(cc)
         %check each lap range
         for ll=1:size(lap,2)
             cue_check_temp = find(first_cue_time(cc) > lap{ll}(1) &  first_cue_time(cc) < lap{ll}(2));
             if ~isempty(cue_check_temp)
                 cue_registered(ll) = 1;
             end
         end
     end
     
     %check if cue was registered on each lap and print result
     nb_first_cue_registered = sum(cue_registered);
     
     
     %diplay if the first cue (odor off was registered on each lap) 
     if (nb_first_cue_registered == nb_complete_laps)
         disp('First cues (odor off/reward reset) registered on all complete laps.');
     else
         disp('First cue not registered on all all laps !!!'); 
     end
     
     %% 
     
     
     
     
     
     
 end



 %% Plot the textures on the restricted laps
        
 figure;
 hold on
 title('Texture locations on each lap');
 %plot time (seconds) against normalized position
 plot(Behavior.time,Behavior.normalizedposition,'k');
 xlabel('Time [s]');
 ylabel('Position');
 
 %plot RFID crossing indices against time
 plot(Behavior.time,Behavior.lap_binary, 'm');
 
 %plot each microtexture (should be 4)
 for ii=1:size(time_tex,2)
     stem(time_tex{ii}, position_norm_tex{ii}, '*r');
 end
 
 %plot the lap index - lap # scaled down by a factor of 20 - each
 %step 0.05
 plot(Behavior.time,Behavior.lapNb./20,'b');
 

 %% Save to textures and reward substruct
 %normalized position of each texture
 Behavior.textures.position_norm = position_norm_tex;
 %cm position of each texture
 Behavior.textures.position = position_tex;
 %behavior timestamp of each texture
 Behavior.textures.time = time_tex;
 
 %median of the texture positions
 Behavior.textures.position_med = position_med;
  
 %reward location for OCGOL
 if ~isempty(rewardIdx)
     %early reward indices
    Behavior.reward.reward_early = reward_early;
    %late reward indices
    Behavior.reward.reward_late = reward_late;
    %index of within texture cell that corresponds to reward
    Behavior.reward.rewardIdx = rewardIdx;
 end
 
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
 
 
end

