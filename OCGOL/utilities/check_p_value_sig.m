function [star_level] = check_p_value_sig(p_val_list)
%enter p-value list and return star significance

for ii=1:numel(p_val_list)
    
    %temp p-value being checked 
    p_val = p_val_list(ii);
    
    if p_val >= 0.05
        disp('P value not significant')
        star_level{ii,1} = 'n.s.'; 
    elseif (p_val < 0.05 && p_val >= 0.01)
        disp('1 star significance')
        star_level{ii,1} = '*'; 
    elseif (p_val < 0.01 && p_val >= 0.001)
        disp('2 star significance')
        star_level{ii,1} = '**'; 
    elseif p_val < 0.001
        disp('3 star significance')
        star_level{ii,1} = '***'; 
    end
    
end

end

