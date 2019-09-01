function [performance] = session_performance(session_vars, options)

%% TODO: Adjust for punish trials on punish days and excluded trials due to technical issues

%A correct (2) divided by all A trials (2 - correct A trials, 20 - incorrect A trials))
for ss = options.sessionSelect
    %all performance
    performance(1,ss) = length(find(session_vars{ss}.Behavior.performance.trialCorrect ==1))/size(session_vars{ss}.Behavior.performance.trialCorrect,1);
    
    %fraction of A trials correct - 2nd row
    performance(2,ss) = length(find(session_vars{ss}.Behavior.performance.trialOrder == 2))/...
        size(find(session_vars{ss}.Behavior.performance.trialOrder == 2 | session_vars{ss}.Behavior.performance.trialOrder == 20),1);
    
    %fraction of B trials correct - 3rd row
    performance(3,ss) = length(find(session_vars{ss}.Behavior.performance.trialOrder == 3))/...
        size(find(session_vars{ss}.Behavior.performance.trialOrder == 3 | session_vars{ss}.Behavior.performance.trialOrder == 30),1);
    
end


end

