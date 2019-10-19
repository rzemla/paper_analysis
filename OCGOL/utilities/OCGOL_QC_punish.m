function [lap_id] = OCGOL_QC_punish(Behavior,textures, rewards, trialType,sound)

%% Define variables

lap = Behavior.lap;

 %% Quality control check for OCGOL
 %if behavior is OCGOL
 %blue = Odor A trials, late reward (high voltage trial signal)
 %red = Odor B trials, early reward (low voltage trials signal)
 
     %% # of complete lap should equal number of total reward zones
     %total # of complete laps
     nb_complete_laps = size(Behavior.lap,2);
     % # of total reward zones (should equal to number of complete laps)
     %4/5/19 - make sure that the rewards are from complete laps!
     nb_reward_zones = sum([size(rewards{1}.full_laps.position,1),size(rewards{2}.full_laps.position,1)]);
     
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
     
     %preallocated list of whether cue was registered in that lap
     cue_registered = false(1,nb_complete_laps);
     
     %take the time of each cue and see if it occured inside of 1 of the
     %laps
     for cc=1:size(textures{1}.time,1)
         %first_cue_time(cc)
         %check each lap range
         for ll=1:size(lap,2)
             cue_check_temp = find(textures{1}.time(cc) > lap{ll}(1) &  textures{1}.time(cc) < lap{ll}(2));
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
          
%% Check the signal reward signal went off on each lap and that it matched the trial type signal on that lap
     %preallocated list of whether cue was registered in that lap
     reward_registered = zeros(1,nb_complete_laps);
     
     %take the time of each cue and see if it occured inside of 1 of the
     %laps
     %for each reward
     for rr=1:size(rewards,2)
         %for each reward time
     for cc=1:size(rewards{rr}.time,1)
         %check within each lap time range
         for ll=1:size(lap,2)
             reward_check_temp = find(rewards{rr}.time(cc) > lap{ll}(1) &  rewards{rr}.time(cc) < lap{ll}(2));
             if ~isempty(reward_check_temp)
                 %if B trial
                 if rewards{rr}.position < 100
                 reward_registered(ll) = 3;
                 %if A trial
                 elseif rewards{rr}.position > 100
                     reward_registered(ll) = 2;
                 end
             end
         end
     end
     end
     
     %check if reward signal was registered on each lap and print result
     all_laps_register_reward = find(reward_registered == 0);
     
     %trial type identity based on reward signal
     lap_id_reward_based = reward_registered;

     %diplay if the first cue (odor off was registered on each lap) 
     if isempty((all_laps_register_reward))
         disp('Reward cue registered on all laps');
         A_laps_nb = length(find(reward_registered == 2));
         B_laps_nb = length(find(reward_registered == 3));
         
         strA = sprintf('A laps run: %d', A_laps_nb);
         strB = sprintf('B laps run: %d', B_laps_nb); 
         disp(strA);
         disp(strB);

     else
         disp('Reward cue not registered on all all laps !!!');
     end
     
     %% Check that trial type went off on each lap and that is matches the reward signals
     %preallocated list of whether cue was registered in that lap
     trial_registered = zeros(1,nb_complete_laps);
     
     %take the time of each cue and see if it occured inside of 1 of the
     %laps
     %for each trial type
     for rr=1:size(trialType.time,2)
         %for each reward time
     for cc=1:size(trialType.time{rr},1)
         %check within each lap time range
         for ll=1:size(lap,2)
             trial_check_temp = find(trialType.time{rr}(cc) > lap{ll}(1) &  trialType.time{rr}(cc) < lap{ll}(2));
             if ~isempty(trial_check_temp)
                 %if B trial
                 if (mean(trialType.voltage{rr}) < 3.5) && (mean(trialType.voltage{rr}) > 2)
                 trial_registered(ll) = 3;
                 %if A trial
                 elseif mean(trialType.voltage{rr}) > 3.5
                     trial_registered(ll) = 2;
                 %if P trial
                 elseif mean(trialType.voltage{rr}) < 2
                      trial_registered(ll) = 1;
                 end
             end
         end
     end
     end
     
     %check if reward signal was registered on each lap and print result
     all_laps_register_trial= find(trial_registered == 0);
     
     %trial type identity based on reward signal
     lap_id_trial_based = trial_registered;

     %diplay if the first cue (odor off was registered on each lap) 
     if isempty((all_laps_register_trial))
         disp('Trial signal registered on all laps');
         A_laps_nb_trial = length(find(trial_registered == 2));
         B_laps_nb_trial = length(find(trial_registered == 3));
         P_laps_nb_trial = length(find(trial_registered == 1));
         
         strA_tr = sprintf('A laps run (trials): %d', A_laps_nb_trial);
         strB_tr = sprintf('B laps run (trials): %d', B_laps_nb_trial);
         strP_tr = sprintf('P laps run (trials): %d', P_laps_nb_trial); 
         all_tr = sprintf('All complete laps run (trials): %d', P_laps_nb_trial+B_laps_nb_trial+A_laps_nb_trial); 
         disp(strA_tr);
         disp(strB_tr);
         disp(strP_tr);
         disp(all_tr)

     else
         disp('Trial cue not registered on all all laps !!!');
     end

     % check if reward based and type signal based trial ID agree
     if isequal(lap_id_trial_based, lap_id_reward_based)
         disp('Reward signals agree with trials signals');
     else
         disp('Reward signals do NOT agree with trial signals !!!');
     end
          
     %% Check that sound went off on each lap
     %preallocated list of whether cue was registered in that lap
     sound_registered = false(1,nb_complete_laps);
     
     %take the time of each cue and see if it occured inside of 1 of the
     %laps
     %for each sound time
     for cc=1:size(sound.time,1)
         %check within each lap time range
         for ll=1:size(lap,2)
             sound_check_temp = find(sound.time(cc) > lap{ll}(1) &  sound.time(cc) < lap{ll}(2));
             if ~isempty(sound_check_temp)
                 sound_registered(ll) = 1;
             end
         end
     end

     %check if cue was registered on each lap and print result
     nb_sound_registered = sum(sound_registered);

     %diplay if the first cue (odor off was registered on each lap) 
     if (nb_sound_registered == nb_complete_laps)
         disp('Sound registered on all laps.');
     else
         disp('Sound NOT registered on all all laps !!!'); 
     end
     %% Export the lap IDs based trial and reward signal
     lap_id.trial_based = lap_id_trial_based;
     lap_id.reward_based = lap_id_reward_based;
     
end

