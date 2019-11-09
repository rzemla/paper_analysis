function [sig_lvl] = get_star_sig(p_vec)
%get star significance
%-1 - ns - >= 0.05
%1 * - <0.05, >= 0.01
%2 ** - <0.01 >= 0.001
%3 *** - <0.001

%%
sig_lvl = [0 0 0 0];

sig_lvl(find(p_vec >= 0.05)) = -1;
sig_lvl(find(p_vec < 0.05 & p_vec >= 0.01)) = 1;
sig_lvl(find(p_vec < 0.01 & p_vec >= 0.001)) = 2;
sig_lvl(find(p_vec < 0.001 )) = 3;

end

