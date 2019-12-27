function [reward_positions] = extract_reward_position(session_vars)

%normalized position
mean_rew_B_norm = mean(session_vars{1}.Behavior.rewards{1}.position_norm);
mean_rew_A_norm = mean(session_vars{1}.Behavior.rewards{2}.position_norm);

mean_rew_B_pos = mean(session_vars{1}.Behavior.rewards{1}.position);
mean_rew_A_pos = mean(session_vars{1}.Behavior.rewards{2}.position);

mean_track_end_cm = mean(session_vars{1}.Behavior.position_lap(:,2));

reward_positions.mean_rew_A_norm = mean_rew_A_norm;
reward_positions.mean_rew_B_norm = mean_rew_B_norm;

reward_positions.mean_rew_A_pos = mean_rew_A_pos;
reward_positions.mean_rew_B_pos = mean_rew_B_pos;

reward_positions.mean_track_end_cm = mean_track_end_cm;

end

