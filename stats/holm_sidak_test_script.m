%p_in = [0.027068745, 0.036846062];
%p_in = [0.072951305,0.073562893];
p_in = [0.031465513,0.633431891];
%p_in = [0.020629411, 0.034312945];


%set of p values to compare with Prism output
p_adj = holm_sidak_p_adj(p_in,2,0.05)