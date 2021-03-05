function [outputArg1,outputArg2] = kstest2_mult_compare(input_data)

nb_entries = 1;

fig_num = repmat(3,nb_entries,1);
fig_sub = string(repmat('h',nb_entries,1));
data_agg = string(repmat('pooled',nb_entries,1));
comp_descrip = {'Distribution of partially remapping fields between A and B trials'};
n_sample = string([num2str(numel(partial_A(:,2))),' vs ', num2str(numel(partial_B(:,2)))]);
test_name = repmat({'2-sample Kolmogorov–Smirnov test'},nb_entries,1);

n_dof = string(repmat('N/A', nb_entries,1));
test_statistic = [ks2stat]';
adj_method = string(repmat('N/A', nb_entries,1));
p_all = [p_ks2]';
p_adj = string(repmat('N/A', nb_entries,1));
sig_level = check_p_value_sig(p_all);

%create table
t_2ks_partial_pf_dist = table(fig_num, fig_sub, data_agg, comp_descrip, n_sample,...
            test_name, n_dof, test_statistic, p_all, p_adj, adj_method, sig_level,...
            'VariableNames',{'Figure','Subfigure','Data aggregation',...
            'Comparison','N', 'Test', 'Degrees of Freedom', 'Test statistic',...
            'p-value', 'p-value adjusted', 'Adjustment method','Significance'});



end

