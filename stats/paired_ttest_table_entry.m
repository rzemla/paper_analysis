function [t_out] = paired_ttest_table_entry(data_input,...
        num_in, sub_in, data_agg, comp_descrip_in)

%number of repeat statistical test
nb_entries = size(data_input,1);

%% format entry
fig_num = repmat(num_in,nb_entries,1);
fig_sub = string(repmat(sub_in,nb_entries,1));
data_agg = string(repmat(data_agg,nb_entries,1));
comp_descrip = comp_descrip_in;
n_sample = data_input(:,3);
test_name = repmat({'Paired t-test'},nb_entries,1);
n_dof = data_input(:,4);
test_statistic = data_input(:,2);
adj_method = repmat(['Holm-Sidak (', num2str(nb_entries), '-way)'] , nb_entries,1);
p_all = data_input(:,1);
p_adj = holm_sidak_p_adj(p_all',numel(p_all),0.05);
sig_level = check_p_value_sig(p_adj);


%% create t test table
t_out = table(fig_num, fig_sub, data_agg, comp_descrip, n_sample,...
            test_name, n_dof, test_statistic,p_all,p_adj', adj_method, sig_level,...
            'VariableNames',{'Figure','Subfigure','Data aggregation',...
            'Comparison','N', 'Test', 'Degrees of Freedom', 'Test statistic',...
            'p-value', 'p-value adjusted', 'Adjustment method','Significance'});


end

