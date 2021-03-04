function [t_out] = kruskal_wallis_single_table_entry(num_in,sub_in,data_agg,...
                comp_descrip, stats_in)

%number of entries - support for 1 entry now
nb_entries = 1;

%% Assign vars
fig_num = repmat(num_in,nb_entries,1);
fig_sub = string(repmat(sub_in,nb_entries,1));
data_agg = string(repmat(data_agg,nb_entries,1));
n_sample = strjoin(string(stats_in(3:(end-1))));
test_name = repmat({'Kruskal-Wallis'},nb_entries,1);
n_dof =stats_in(end);
test_statistic = stats_in(2);

adj_method = string(repmat('N/A', nb_entries,1));
p_all = [stats_in(1)]';
p_adj = string(repmat('N/A', nb_entries,1));
sig_level = check_p_value_sig(p_all);

%% create table
t_out = table(fig_num, fig_sub, data_agg, comp_descrip, n_sample,...
            test_name, n_dof, test_statistic, p_all, p_adj, adj_method, sig_level,...
            'VariableNames',{'Figure','Subfigure','Data aggregation',...
            'Comparison','N', 'Test', 'Degrees of Freedom', 'Test statistic',...
            'p-value', 'p-value adjusted', 'Adjustment method','Significance'});

end

