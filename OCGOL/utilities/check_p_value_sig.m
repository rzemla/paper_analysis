function [outputArg1,outputArg2] = check_p_value_sig(p_val)
%enter p-value and return star significance

if p_val >= 0.05
    disp('P value not significant')
elseif (p_val < 0.05 && p_val >= 0.01)
    disp('1 star significance')
elseif (p_val < 0.01 && p_val >= 0.001)
    disp('2 star significance')
elseif p_val < 0.001
    disp('3 star significance')
end


end

