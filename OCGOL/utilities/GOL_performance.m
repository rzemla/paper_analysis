function [Behavior] = GOL_performance(Behavior,CSV)
%fraction of licks in reward zone
%fraction of lick on anticipatory zone
%plot histogram across bins

%% Import behavior variables

position = Behavior.position;
time = Behavior.time;
%first lap start time
first_lap_start_time = Behavior.lap{1}(1);
%last lap end time
final_lap_end_time = Behavior.lap{end}(2);
%textures
textures = Behavior.textures;


%% Asssign lick signal and binarize

%raw lick signal channel
lickSignal = CSV(:,8);
%reward signal - high whenever reward droplet is dispensed
rewardSignal = CSV(:,7);

%% Digitize lick signal
%select signal only from left lickport and binarize
digi_idx = find(lickSignal > 3.4);
%create blank 0 vector
binary_lick = zeros(size(lickSignal,1),1);
%insert licks on signal on left port
binary_lick(digi_idx) = 1;

%% Digitize reward signal
%select signal only from left lickport and binarize
digi_rew_idx = find(rewardSignal > 4.5);
%create blank 0 vector
binary_rew = zeros(size(rewardSignal,1),1);
%insert licks on signal on left port
binary_rew(digi_rew_idx) = 1;

%% Get onset of each lick (low to high transition)

lick_on_idx = find(diff(binary_lick) == 1);
lick_on_idx = lick_on_idx + 1;

lick_on_timestamps = time(lick_on_idx);

%create final binary onset lick vector
binary_lick_on = zeros(size(lickSignal,1),1);
binary_lick_on(lick_on_idx) = 1;

%% Get onset of each reward (low to high transition)

rew_on_idx = find(diff(binary_rew) == 1);
rew_on_idx = rew_on_idx + 1;

rew_on_timestamps = time(rew_on_idx);

%create final binary onset reward vector
binary_rew_on = zeros(size(rewardSignal,1),1);
binary_rew_on(rew_on_idx) = 1;

%% Get position,time of each lick onset

%time (copied to another variable) and position of each lick
lick_time = lick_on_timestamps;
lick_position = position(lick_on_idx);

%restrict licks to only complete laps
lick_time_R_idx = find(lick_time >= first_lap_start_time & lick_time <= final_lap_end_time);
lick_time_R = lick_time(lick_time_R_idx);
lick_position_R = lick_position(lick_time_R_idx);

%% Get position,time of each reward onset

%time (copied to another variable) and position of each reward
rew_time = rew_on_timestamps;
rew_position = position(rew_on_idx);

%restrict licks to only complete laps
rew_time_R_idx = find(rew_time >= first_lap_start_time & rew_time <= final_lap_end_time);
rew_time_R = rew_time(rew_time_R_idx);
rew_position_R = rew_position(rew_time_R_idx);

%% Find reward location either for Block 1 or Block 2

%determine whether block 1 or block 2 GOL
if ~isempty(find(textures.position_med >= 54 & textures.position_med <= 70))
    disp('GOL Block 1 reward location discovered')
    reward_idx = find(textures.position_med >= 54 & textures.position_med <= 70);
    %reward loc;
    reward_start_loc = textures.position_med(reward_idx);
    
elseif ~isempty(find(textures.position_med >= 136 & textures.position_med <= 148))
    disp('GOL Block 2 reward location discovered')
    reward_idx = find(textures.position_med >= 136 & textures.position_med <= 148);
    %reward loc;
    reward_start_loc = textures.position_med(reward_idx);
else
    disp('No reward location discovered - random foraging')
    %give block 1 reward as location (RF)
    reward_start_loc = 58.6;
end

%% Find number of licks in respective zones

reward_zone_licks_idx = find(lick_position_R >= reward_start_loc & lick_position_R <= (reward_start_loc + 10));
ant_zone_licks_idx = find(lick_position_R >= (reward_start_loc-10) & lick_position_R < reward_start_loc);

%fraction of licks in reward zone
frac_reward_zone_licks = length(reward_zone_licks_idx)/length(lick_position_R);
%fraction of lick in anticipatory zone
frac_ant_zone_licks = length(ant_zone_licks_idx)/length(lick_position_R);

%dislay fraction of licks in anticipatory zone and reward zone
disp(sprintf('Fraction of licks in ranticipatory zone: %f', frac_ant_zone_licks))
disp(sprintf('Fraction of licks in reward zone: %f', frac_reward_zone_licks))

%% Plot lick distributions

%how many spatial bins
nb_spatial_bins = 100;

figure;
subplot(2,1,1)
hold on
xlim([0 200])
title('Distribution of licks across entire session');
h = histogram(lick_position,nb_spatial_bins,'Normalization','probability');
%ylabel('Lick Count')
ylabel('Normalized density');
ylim([0 0.5]);
xlabel('Binned position [cm]')
%reward start
stem(reward_start_loc, 0.5,'g');
%reward end
stem(reward_start_loc + 10, 0.5,'g');
%ant start
stem(reward_start_loc-10, 0.5,'m');
hold off

subplot(2,1,2)
hold on
xlim([0 200])
title('Distribution of licks across complete laps');
histogram(lick_position_R,nb_spatial_bins,'Normalization','probability');
%ylabel('Lick Count')
ylabel('Normalized density');
ylim([0 0.5]);
xlabel('Binned position [cm]')
%reward start
stem(reward_start_loc, 0.5,'g');
%reward end
stem(reward_start_loc + 10, 0.5,'g');
%ant start
stem(reward_start_loc-10, 0.5,'m');
hold off

%% Plot signals

figure;
subplot(3,1,1);
hold on
plot(time,position,'k');
hold off

subplot(3,1,2);
hold on
%check align of selected onset with entire lick signal
%onsets
%stem(time, 2*binary_lick_on, 'r');
stem(lick_on_timestamps,2*ones(size(lick_on_timestamps,1),1),'r')
%all
plot(time, binary_lick,'b');
hold off

subplot(3,1,3);
hold on
%plot reward signal against time
plot(time,position,'k');
%stem(rew_on_timestamps,2*ones(size(rew_on_timestamps,1),1),'r')
scatter(rew_time_R,rew_position_R,'r')
hold off


%% Save output to Behavior structure

%reward location
Behavior.performance.reward_loc = reward_start_loc;
%fraction of licks in reward zone
Behavior.performance.frac_rew = frac_reward_zone_licks;
%fraction of licks in the anticipatory zone
Behavior.performance.frac_ant = frac_ant_zone_licks;


end

