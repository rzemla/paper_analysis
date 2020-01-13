function [performance,session_lap_ct,session_lap_ct_corr] = session_performance(session_vars, options)

%% TODO: Adjust for punish trials on punish days and excluded trials due to technical issues
%should already be done
%technical exclusion adjustment


%A correct (2) divided by all A trials (2 - correct A trials, 20 - incorrect A trials))
for ss = options.sessionSelect
    %get total number of laps after adjustment for technical exclusion
    %only include correct and incorrect trials (not punish or technically
    %excluded trials)
    nb_laps_adj_tech(ss) = sum((session_vars{ss}.Behavior.performance.trialCorrect == 0) | (session_vars{ss}.Behavior.performance.trialCorrect == 1));
    
    %all performance
    performance(1,ss) = length(find(session_vars{ss}.Behavior.performance.trialCorrect ==1))/nb_laps_adj_tech(ss);
    
    %fraction of A trials correct - 2nd row
    performance(2,ss) = length(find(session_vars{ss}.Behavior.performance.trialOrder == 2))/...
        size(find(session_vars{ss}.Behavior.performance.trialOrder == 2 | session_vars{ss}.Behavior.performance.trialOrder == 20),1);
    
    %fraction of B trials correct - 3rd row
    performance(3,ss) = length(find(session_vars{ss}.Behavior.performance.trialOrder == 3))/...
        size(find(session_vars{ss}.Behavior.performance.trialOrder == 3 | session_vars{ss}.Behavior.performance.trialOrder == 30),1);
    
    %all lap count
    session_lap_ct(1,ss) =nb_laps_adj_tech(ss);
    %all A lap count
    session_lap_ct(2,ss) = size(find(session_vars{ss}.Behavior.performance.trialOrder == 2 | session_vars{ss}.Behavior.performance.trialOrder == 20),1);
    %all B lap count
    session_lap_ct(3,ss) = size(find(session_vars{ss}.Behavior.performance.trialOrder == 3 | session_vars{ss}.Behavior.performance.trialOrder == 30),1);
   
    %corr counts
    %all A lap count
    session_lap_ct_corr(2,ss) = size(find(session_vars{ss}.Behavior.performance.trialOrder == 2),1);
    %all B lap count
    session_lap_ct_corr(3,ss) = size(find(session_vars{ss}.Behavior.performance.trialOrder == 3),1);
    
end


end
