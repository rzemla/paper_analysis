function get_star_sig_single(p_vec)
%get star significance
%-1 - ns - >= 0.05
%1 * - <0.05, >= 0.01
%2 ** - <0.01 >= 0.001
%3 *** - <0.001

%% Return label for single value when input 

if p_vec >= 0.05
    disp('n.s.')
elseif  p_vec < 0.05 && p_vec >= 0.01
    disp('*')
elseif p_vec < 0.01 && p_vec >= 0.001
    disp('**')
elseif p_vec < 0.001
    disp('***')
end

end

