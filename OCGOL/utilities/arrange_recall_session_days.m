function [merge_theta_recall_days] = arrange_recall_session_days(theta_recall)


%% 
%merge recall as a function relative to d1 - d6
for dd=2:9 %for day 1 vs 2 3 4... or equivalent day distance
    
    switch dd
        case 2
            %for each animal
            for aa=1:4 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_recall_days.A{aa,dd} = theta_recall{aa}.A{1, 2};
                merge_theta_recall_days.B{aa,dd} = theta_recall{aa}.B{1, 2};
            end
        case 3
            for aa=1:4 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_recall_days.A{aa,dd} = theta_recall{aa}.A{1, 3};
                 merge_theta_recall_days.B{aa,dd} = theta_recall{aa}.B{1, 3};
            end
        case 4 %session 4 vs. 7; session 3 vs. 4
            for aa=1:4 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_recall_days.A{aa,dd} = [theta_recall{aa}.A{4, 7},theta_recall{aa}.A{3, 4}];
                merge_theta_recall_days.B{aa,dd} = [theta_recall{aa}.B{4, 7},theta_recall{aa}.B{3, 4}];
            end
        case 5 %session 2 vs. 6; session 3 vs. 5
            for aa=1:4 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_recall_days.A{aa,dd} = [theta_recall{aa}.A{2, 6},theta_recall{aa}.A{3, 5}];
                merge_theta_recall_days.B{aa,dd} = [theta_recall{aa}.B{2, 6},theta_recall{aa}.B{3, 5}];
            end    
        case 6 %session 2 vs. 6; session 3 vs. 5
            for aa=1:4 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_recall_days.A{aa,dd} = theta_recall{aa}.A{1, 4};
                merge_theta_recall_days.B{aa,dd} = theta_recall{aa}.B{1, 4};
            end   
        case 7 %session 2 vs. 6; session 3 vs. 5
            for aa=1:4 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_recall_days.A{aa,dd} = theta_recall{aa}.A{1, 5};
                merge_theta_recall_days.B{aa,dd} = theta_recall{aa}.B{1, 5};
            end
        case 8 %session 2 vs. 6; session 3 vs. 5
            for aa=1:4 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_recall_days.A{aa,dd} = theta_recall{aa}.A{1, 6};
                merge_theta_recall_days.B{aa,dd} = theta_recall{aa}.B{1, 6};
            end            
        case 9 %session 2 vs. 6; session 3 vs. 5
            for aa=1:4 %all animals
                %extract all the days that correspond to the time duration
                merge_theta_recall_days.A{aa,dd} = theta_recall{aa}.A{1, 7};
                merge_theta_recall_days.B{aa,dd} = theta_recall{aa}.B{1, 7};
            end  
            
            
    end
    
end


end

