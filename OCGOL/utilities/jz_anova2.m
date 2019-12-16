function [pGroup, pAll, table] = jz_anova2(Y1, Y2, X1, X2)
    %[pGroup, pAll, table] = jz_anova(Y1, Y2, X1, X2)
    % Y1: AUC of group 1
    % Y2: AUC of group 2
    % X1: Speed of group 1
    % X2: Speed of group 2
    
    %dependent variable
    Y = [Y1(:); Y2(:)];   
    %nuisance variable
    X = [X1(:); X2(:)];
    
    %group as 1 vs. 2
    group = [ones(length(Y1),1); 2*ones(length(Y2),1)];
    
    [pAll, table] = anovan(Y, {group X}, 'continuous', 2, 'display', 'off', 'model','full','varnames', char('Group', 'X-factor'));
    pGroup = pAll(1);
    
end


