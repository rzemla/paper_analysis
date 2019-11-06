function [merge_theta_recall_days] = arrange_recall_session_rew_dist_days(theta_reward_zone_recall)


%% 
%merge recall as a function relative to d1 - d6
for dd=2:9 %for day 1 vs 2 3 4... or equivalent day distance
    
    switch dd
        case 2
            %for each animal
            for aa=1:5 %all animals
                %extract all the days that correspond to the time duration                
                merge_theta_recall_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_recall{aa}.A.rewA{1,2},1,2));
                merge_theta_recall_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_recall{aa}.B.rewB{1,2},1,2));
                
            end
        case 3
            for aa=1:5 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_recall_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_recall{aa}.A.rewA{1,3},1,2));
                merge_theta_recall_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_recall{aa}.B.rewB{1,3},1,2));
            end
        case 4 %session 4 vs. 7; session 3 vs. 4 - equivalent day 1 distance (not from start day 1 for recall)
            for aa=1:5 %all animals
                %extract all the days that correspond to the time duration              
                merge_theta_recall_days.A.rewA{aa,dd} = [abs(diff(theta_reward_zone_recall{aa}.A.rewA{4,7},1,2));...
                                                        abs(diff(theta_reward_zone_recall{aa}.A.rewA{3,4},1,2))];
                                                    
                merge_theta_recall_days.B.rewB{aa,dd} = [abs(diff(theta_reward_zone_recall{aa}.B.rewB{4,7},1,2));...
                                                          abs(diff(theta_reward_zone_recall{aa}.B.rewB{3,4},1,2))];
            end
        case 5 %session 2 vs. 6; session 3 vs. 5
            for aa=1:5 %all animals
                %extract all the days that correspond to the time duration
                
                merge_theta_recall_days.A.rewA{aa,dd} = [abs(diff(theta_reward_zone_recall{aa}.A.rewA{2,6},1,2));...
                    abs(diff(theta_reward_zone_recall{aa}.A.rewA{3,5},1,2))];
                
                merge_theta_recall_days.B.rewB{aa,dd} = [abs(diff(theta_reward_zone_recall{aa}.B.rewB{2,6},1,2));...
                    abs(diff(theta_reward_zone_recall{aa}.B.rewB{3,5},1,2))];
                
            end    
        case 6 %session 2 vs. 6; session 3 vs. 5
            for aa=1:5 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_recall_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_recall{aa}.A.rewA{1,4},1,2));
                merge_theta_recall_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_recall{aa}.B.rewB{1,4},1,2));
            end   
        case 7 %session 2 vs. 6; session 3 vs. 5
            for aa=1:5 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_recall_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_recall{aa}.A.rewA{1,5},1,2));
                merge_theta_recall_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_recall{aa}.B.rewB{1,5},1,2));
            end
        case 8 %session 2 vs. 6; session 3 vs. 5
            for aa=1:5 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_recall_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_recall{aa}.A.rewA{1,6},1,2));
                merge_theta_recall_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_recall{aa}.B.rewB{1,6},1,2));
            end            
        case 9 %session 2 vs. 6; session 3 vs. 5
            for aa=1:5 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_recall_days.A.rewA{aa,dd} = abs(diff(theta_reward_zone_recall{aa}.A.rewA{1,7},1,2));
                merge_theta_recall_days.B.rewB{aa,dd} = abs(diff(theta_reward_zone_recall{aa}.B.rewB{1,7},1,2));
            end  
            
            
    end
    
end


end

