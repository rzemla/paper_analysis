function [t_out ] = paired_wilcoxon_table_entry(num_in,sub_in,data_agg_in,...
                comp_descrip, stats_in)

%check if stats_in is a cell
if iscell(stats_in)
    %convert cell input to matrix
    stats_in_mat = cell2mat(stats_in');
    
    %how many entries are generated for the table
    nb_entries = sum(~cellfun(@isempty, stats_in));
    
else
    stats_in_mat = stats_in;
    nb_entries = size(stats_in_mat,1);
end

fig_num = repmat(num_in,nb_entries,1);
fig_sub = string(repmat(sub_in,nb_entries,1));
data_agg = string(repmat(data_agg_in,nb_entries,1));

%number of samples for each test            
n_sample = stats_in_mat(:,3);
%test ran
test_name = string(repmat('Paired Wilcoxon Sign Rank',nb_entries,1));
%degrees of freedom
n_dof = stats_in_mat(:,4);
%test statistic
test_statistic = stats_in_mat(:,2);
adj_method =  string(repmat(['Holm-Sidak (',num2str(nb_entries),'-way)'], nb_entries,1));
p = stats_in_mat(:,1);
p_adj =holm_sidak_p_adj(p',numel(p),0.05);
sig_level = check_p_value_sig(p_adj');

%table
t_out = table(fig_num, fig_sub, data_agg, comp_descrip, n_sample,...
            test_name, n_dof, test_statistic,p,p_adj', adj_method, sig_level,...
            'VariableNames',{'Figure','Subfigure','Data aggregation',...
            'Comparison','N', 'Test', 'Degrees of Freedom', 'Test statistic',...
            'p-value', 'p-value adjusted', 'Adjustment method','Significance'});

end

