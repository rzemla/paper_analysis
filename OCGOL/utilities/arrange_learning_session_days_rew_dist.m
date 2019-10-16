function [merge_theta_learn_days] = arrange_learning_session_days_rew_dist(theta_reward_zone_learn)

%FIRST 6 days for now until all of data is processed

%% 
%merge recall as a function relative to d1 - d6
for dd=2:6 %for day 1 vs 2 3 4... or equivalent day distance
    
    switch dd
        case 2
            %for each animal, do same session extraction
            for aa=1:3 %all animals
                %extract all the days that correspond to the time duration
                %abs difference in angles relative to start of reward zone
                %A for A laps and B for B laps
                merge_theta_learn_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.A.rewA{1,2},1,2));
                merge_theta_learn_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.B.rewB{1,2},1,2));
                
                 
            end
        case 3
            for aa=1:3 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_learn_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.A.rewA{1,3},1,2));
                merge_theta_learn_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.B.rewB{1,3},1,2));
            end
        case 4 
            for aa=1:3 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_learn_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.A.rewA{1,4},1,2));
                merge_theta_learn_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.B.rewB{1,4},1,2));
            end
        case 5 
            for aa=1:3 %all animals
                %extract all the days that correspond to the time duration
                if aa==1 %skip
                merge_theta_learn_days.A.rewA{aa,dd} = [];
                merge_theta_learn_days.B.rewB{aa,dd} = [];
                elseif aa==2 %skip
                    merge_theta_learn_days.A.rewA{aa,dd} = [];
                    merge_theta_learn_days.B.rewB{aa,dd} = [];
                elseif aa==3
                merge_theta_learn_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.A.rewA{1,5},1,2));
                merge_theta_learn_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.B.rewB{1,5},1,2));
                end

            end    
        case 6 %session 2 vs. 6; session 3 vs. 5
            for aa=1:3 %all animals
                if aa==1 %skip
                merge_theta_learn_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.A.rewA{1,5},1,2));
                merge_theta_learn_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.B.rewB{1,5},1,2));
                elseif aa==2 %skip
                merge_theta_learn_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.A.rewA{1,5},1,2));
                merge_theta_learn_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.B.rewB{1,5},1,2));
                elseif aa==3
                merge_theta_learn_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.A.rewA{1,6},1,2));
                merge_theta_learn_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.B.rewB{1,6},1,2));
                end
            end   
        case 7 %session 2 vs. 6; session 3 vs. 5
            for aa=1:3 %all animals
                if aa==1 %skip
                merge_theta_learn_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.A.rewA{1,6},1,2));
                merge_theta_learn_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.B.rewB{1,6},1,2));
                elseif aa==2 %skip
                merge_theta_learn_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.A.rewA{1,6},1,2));
                merge_theta_learn_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.B.rewB{1,6},1,2));
                elseif aa==3
                merge_theta_learn_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.A.rewA{1,7},1,2));
                merge_theta_learn_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.B.rewB{1,7},1,2));
                end
            end
        case 8 %session 2 vs. 6; session 3 vs. 5
            for aa=1:3 %all animals
                if aa==1 %skip
                merge_theta_learn_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.A.rewA{1,7},1,2));
                merge_theta_learn_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.B.rewB{1,7},1,2));
                elseif aa==2 %skip
                merge_theta_learn_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.A.rewA{1,7},1,2));
                merge_theta_learn_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.B.rewB{1,7},1,2));
                elseif aa==3
                merge_theta_learn_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.A.rewA{1,8},1,2));
                merge_theta_learn_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_learn{aa}.B.rewB{1,8},1,2));
                end
            end            
%         case 9 %session 2 vs. 6; session 3 vs. 5
%             for aa=1:3 %all animals
%                 %extract all the days that correspond to the time duration
%                 merge_theta_learn_days.A{aa,dd} = theta_learn{aa}.A{1, 7};
%                 merge_theta_learn_days.B{aa,dd} = theta_learn{aa}.B{1, 7};
%             end  
%             
            
    end
    
end


end

