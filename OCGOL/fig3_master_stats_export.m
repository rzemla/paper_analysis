%% Master script for generating stats for paper
%source data, import and format Prism analysis, and export to Excel, and
%Word legend data

%import data for each figure

%figure 2 source data
load('G:\Google_drive\task_selective_place_paper\matlab_data\source_data_fig3.mat')

%laptop directory
%load('C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data\source_data_fig3.mat')


%% Figure 3d Fraction of each class of remapping neurons

%unload fractional data
frac_class_mean = source_data_AB_remap.frac_class_mean;
class_names = source_data_AB_remap.class_names;

%Friedman test for all groups
[friedman_stats] = friedman_test(frac_class_mean);

%create single entry table for Friedman output
nb_entries = 1;

fig_num = repmat(3,nb_entries,1);
fig_sub = string(repmat('d',nb_entries,1));
data_agg = string(repmat('by animal',nb_entries,1));
comp_descrip = {'Fraction of remapping neuron subtypes (common, activity, global, partial, unclassified)';};
n_sample = [friedman_stats(3)]';
test_name = repmat({'Friedman test'},nb_entries,1);
n_dof = friedman_stats(4);
test_statistic = [friedman_stats(2)]';
adj_method = string(repmat('N/A', nb_entries,1));
p_all = [friedman_stats(1)]';
p_adj = string(repmat('N/A', nb_entries,1));
sig_level = check_p_value_sig(p_all);

%create table
t_frac_friedman = table(fig_num, fig_sub, data_agg, comp_descrip, n_sample,...
            test_name, n_dof, test_statistic, p_all, p_adj, adj_method, sig_level,...
            'VariableNames',{'Figure','Subfigure','Data aggregation',...
            'Comparison','N', 'Test', 'Degrees of Freedom', 'Test statistic',...
            'p-value', 'p-value adjusted', 'Adjustment method','Significance'});

%paired Wilcoxon for each comparison

%common vs activity
stats.com_act = paired_wilcoxon_signrank(frac_class_mean(:,1),frac_class_mean(:,2));

%common vs. global
stats.com_glo = paired_wilcoxon_signrank(frac_class_mean(:,1),frac_class_mean(:,3));

%common vs. partial
stats.com_par = paired_wilcoxon_signrank(frac_class_mean(:,1),frac_class_mean(:,4));

%create excel importable table data

%create AUC/min table
%Figure, Subfigure, Data aggregation, Comparison, N, Test, Degrees of Freedom, 
%Test statistic, p-value, p-value adjusted, ad. method, Significance

%number of repeat statistical test
nb_entities = 3;
%figure number
fig_num = repmat(3,nb_entities,1);
%subplot
fig_sub = repmat('d',nb_entities,1);
%how was data aggregated
data_agg = repmat('by animal',nb_entities,1);
%what comparison is being made
comp_descrip = {'Common vs. activity';...
                'Common vs. global';...
                'Common vs. partial'};
%number of subjects            
%which statistic to extract for each test
stat_idx =3;            
n_sample = [stats.com_act(stat_idx), stats.com_glo(stat_idx),stats.com_par(stat_idx)]';
%Test name
test_name = repmat('Paired Wilcoxon Sign Rank',nb_entities,1);

%degrees of freedom
stat_idx =4; 
n_dof = [stats.com_act(stat_idx), stats.com_glo(stat_idx),stats.com_par(stat_idx)]';
%test statistic value
stat_idx =2; 
test_statistic = [stats.com_act(stat_idx), stats.com_glo(stat_idx),stats.com_par(stat_idx)]';
%multiple comparison adjustment method
adj_method = repmat('Holm-Sidak (3-way)', nb_entities,1);
%p values
stat_idx = 1;
p = [stats.com_act(stat_idx), stats.com_glo(stat_idx),stats.com_par(stat_idx)]';
%mutiple comparisons adjusted p values
p_adj = holm_sidak_p_adj(p',numel(p),0.05);
%star statistical significance level
sig_level = check_p_value_sig(p_adj);

%create RUN AUC/min table
t_frac_remap = table(fig_num, fig_sub, data_agg, comp_descrip, n_sample,...
            test_name, n_dof, test_statistic,p,p_adj', adj_method, sig_level,...
            'VariableNames',{'Figure','Subfigure','Data aggregation',...
            'Comparison','N', 'Test', 'Degrees of Freedom', 'Test statistic',...
            'p-value', 'p-value adjusted', 'Adjustment method','Significance'});

%% Figure 3e Distribution of common neurons on track
%Rayleigh test of circular uniformity
%unload data
common_pf_dist = source_data_AB_remap.common_pf_dist;
bin_spacing = source_data_AB_remap.bin_spacing;
center_bin_radians = source_data_AB_remap.center_bin_radians;
centroid_count_com = source_data_AB_remap.centroid_count_com;

%use this result since you downbin to 25 bins
%test for A 
%run rayleigh test - similar result if even bin spacing is input or not
[pval_cpf, z_cpf] = circ_rtest(center_bin_radians,centroid_count_com ,bin_spacing);

%Test statistic, p-value, p-value adjusted, ad. method, Significance

nb_entries = 1;

fig_num = repmat(3,nb_entries,1);
fig_sub = string(repmat('e',nb_entries,1));
data_agg = string(repmat('pooled',nb_entries,1));
comp_descrip = {'Distribution of PF centroid locations in common neurons (25 bins)';};
n_sample = [numel(common_pf_dist)]';
test_name = repmat({'Rayleigh test of circular uniformity'},nb_entries,1);
n_dof = string(repmat('N/A', nb_entries,1));
test_statistic = [z_cpf]';
adj_method = string(repmat('N/A', nb_entries,1));
p_all = [pval_cpf]';
p_adj = string(repmat('N/A', nb_entries,1));
sig_level = check_p_value_sig(p_all);

%create table
t_pf_dist_com = table(fig_num, fig_sub, data_agg, comp_descrip, n_sample,...
            test_name, n_dof, test_statistic, p_all, p_adj, adj_method, sig_level,...
            'VariableNames',{'Figure','Subfigure','Data aggregation',...
            'Comparison','N', 'Test', 'Degrees of Freedom', 'Test statistic',...
            'p-value', 'p-value adjusted', 'Adjustment method','Significance'});


%% Figure 3f Place field shifts for global remappers on A vs B trials
%1-sample Wilxocon Signed Rank test against 0

zoneI_shift = source_data_AB_remap.global_shift.zoneI_diff;
zoneII_shift = source_data_AB_remap.global_shift.zoneII_diff;
zoneIII_shift = source_data_AB_remap.global_shift.zoneIII_diff;

%paired wilcoxon A, B, AB sel on A vs B laps during RUN epochs
%I
stats.global_shiftI = paired_wilcoxon_signrank(zoneI_shift,0);

%II
stats.global_shiftII = paired_wilcoxon_signrank(zoneII_shift,0);

%III
stats.global_shiftIII = paired_wilcoxon_signrank(zoneIII_shift,0);

nb_entries = 3;

fig_num = repmat(3,nb_entries,1);
fig_sub = repmat('f',nb_entries,1);
data_agg = repmat('by animal',nb_entries,1);

%description of comparison
comp_descrip = {'Zone I (B before A) vs. 0';...
                'Zone II (B before A) vs. 0';...
                'Zone III (B before A) vs. 0'};
%number of samples for each test            
n_sample = [stats.global_shiftI(3) stats.global_shiftII(3) stats.global_shiftIII(3)]';
%test ran
test_name = repmat('1-sample Wilcoxon Sign Rank',nb_entries,1);
%degrees of freedom
n_dof = [stats.global_shiftI(4) stats.global_shiftII(4) stats.global_shiftIII(4)]';
%test statistic
test_statistic = [stats.global_shiftI(2) stats.global_shiftII(2) stats.global_shiftIII(2)]';
adj_method = repmat('N/A', nb_entries,1);
p = [stats.global_shiftI(1) stats.global_shiftII(1) stats.global_shiftIII(1)];
p_adj = repmat('N/A', nb_entries,1);
sig_level = check_p_value_sig(p);

%create noRUN AUC/min table
t_zone_shift = table(fig_num, fig_sub, data_agg, comp_descrip, n_sample,...
            test_name, n_dof, test_statistic,p',p_adj, adj_method, sig_level,...
            'VariableNames',{'Figure','Subfigure','Data aggregation',...
            'Comparison','N', 'Test', 'Degrees of Freedom', 'Test statistic',...
            'p-value', 'p-value adjusted', 'Adjustment method','Significance'});
        
%% Figure 3g Remapping of global remappers between reward zones on A vs B trials
%unload data
frac_zones = source_data_AB_remap.zone_shifts.frac_zones;
%labels for moving between reward zones for global remappers
frac_shift_labels = source_data_AB_remap.zone_shifts.frac_shift_labels;

%Friedman, 
%Friedman test for all groups
[friedman_stats] = friedman_test(frac_zones);

%create single entry table for Friedman output
nb_entries = 1;

fig_num = repmat(3,nb_entries,1);
fig_sub = string(repmat('g',nb_entries,1));
data_agg = string(repmat('by animal',nb_entries,1));
comp_descrip = {'Global remapping neurons switching fields between A or B laps between any of three zones I,II,III';};
n_sample = [friedman_stats(3)]';
test_name = repmat({'Friedman test'},nb_entries,1);
n_dof = friedman_stats(4);
test_statistic = [friedman_stats(2)]';
adj_method = string(repmat('N/A', nb_entries,1));
p_all = [friedman_stats(1)]';
p_adj = string(repmat('N/A', nb_entries,1));
sig_level = check_p_value_sig(p_all);

%create table
t_frac_zones_friedman = table(fig_num, fig_sub, data_agg, comp_descrip, n_sample,...
            test_name, n_dof, test_statistic, p_all, p_adj, adj_method, sig_level,...
            'VariableNames',{'Figure','Subfigure','Data aggregation',...
            'Comparison','N', 'Test', 'Degrees of Freedom', 'Test statistic',...
            'p-value', 'p-value adjusted', 'Adjustment method','Significance'});
        
%paired Wilcoxon
%  1 {'AI and BI'}  
%  2 {'AII and BII'}  
%  3 {'AIII and BIII'}  
%  4 {'AI and BII'}   
%  5 {'BI and AII'}
%  
%  6 {'AII and BIII'}  
%  7 {'BII and AIII'}   
%  8 {'AI and BIII'}  
%  9 {'BI and AIII'}
 
%1 {'AI and BI'} vs. 2 {'AII and BII'} 
stats.zone1v2 = paired_wilcoxon_signrank(frac_zones(:,1),frac_zones(:,2));

%1 {'AI and BI'}  vs.  3 {'AIII and BIII'}
stats.zone1v3 = paired_wilcoxon_signrank(frac_zones(:,1),frac_zones(:,3));

%2 {'AII and BII'} vs. 3 {'AIII and BIII'}
stats.zone2v3 = paired_wilcoxon_signrank(frac_zones(:,2),frac_zones(:,3));

% 4 {'AI and BII'} vs 5 {'BI and AII'}
stats.zone4v5 = paired_wilcoxon_signrank(frac_zones(:,4),frac_zones(:,5));

% 6 {'AII and BIII'} vs 7 {'BII and AIII'} 
stats.zone6v7 = paired_wilcoxon_signrank(frac_zones(:,6),frac_zones(:,7));

%  8 {'AI and BIII'}  vs 9 {'BI and AIII'}
stats.zone8v9 = paired_wilcoxon_signrank(frac_zones(:,8),frac_zones(:,9));

%create excel importable table data
[stats.zone1v2(stat_idx), stats.zone1v3(stat_idx), stats.zone2v3(stat_idx),...
    stats.zone4v5(stat_idx), stats.zone6v7(stat_idx), stats.zone8v9(stat_idx)];

%create AUC/min table
%Figure, Subfigure, Data aggregation, Comparison, N, Test, Degrees of Freedom, 
%Test statistic, p-value, p-value adjusted, ad. method, Significance

%number of repeat statistical test
nb_entities = 6;
%figure number
fig_num = repmat(3,nb_entities,1);
%subplot
fig_sub = repmat('g',nb_entities,1);
%how was data aggregated
data_agg = repmat('by animal',nb_entities,1);
%what comparison is being made
comp_descrip = {'AI-BI vs. AII-BII';...
                'AI-BI vs. AIII-BIII';...
                'AII-BII vs. AIII-BIII';...
                'AI-BII vs. AII-BI';...
                'AII-BIII vs. AIII-BII';...
                'AI-BIII vs. AIII-BI'};
%number of subjects            
%which statistic to extract for each test
stat_idx =3;            
n_sample = [stats.zone1v2(stat_idx), stats.zone1v3(stat_idx), stats.zone2v3(stat_idx),...
    stats.zone4v5(stat_idx), stats.zone6v7(stat_idx), stats.zone8v9(stat_idx)]';
%Test name
test_name = repmat('Paired Wilcoxon Sign Rank',nb_entities,1);

%degrees of freedom
stat_idx =4; 
n_dof = [stats.zone1v2(stat_idx), stats.zone1v3(stat_idx), stats.zone2v3(stat_idx),...
    stats.zone4v5(stat_idx), stats.zone6v7(stat_idx), stats.zone8v9(stat_idx)]';
%test statistic value
stat_idx =2; 
test_statistic = [stats.zone1v2(stat_idx), stats.zone1v3(stat_idx), stats.zone2v3(stat_idx),...
    stats.zone4v5(stat_idx), stats.zone6v7(stat_idx), stats.zone8v9(stat_idx)]';
%multiple comparison adjustment method
adj_method = repmat('Holm-Sidak (6-way)', nb_entities,1);
%p values
stat_idx = 1;
p = [stats.zone1v2(stat_idx), stats.zone1v3(stat_idx), stats.zone2v3(stat_idx),...
    stats.zone4v5(stat_idx), stats.zone6v7(stat_idx), stats.zone8v9(stat_idx)]';
%mutiple comparisons adjusted p values
p_adj = holm_sidak_p_adj(p',numel(p),0.05);
%star statistical significance level
sig_level = check_p_value_sig(p_adj);

%create RUN AUC/min table
t_frac_zone_remap = table(fig_num, fig_sub, data_agg, comp_descrip, n_sample,...
            test_name, n_dof, test_statistic,p,p_adj', adj_method, sig_level,...
            'VariableNames',{'Figure','Subfigure','Data aggregation',...
            'Comparison','N', 'Test', 'Degrees of Freedom', 'Test statistic',...
            'p-value', 'p-value adjusted', 'Adjustment method','Significance'});


%% Figure 3h Distribution of partial remapping place fields between partial  A vs. partial B place fields
%2-sample KS test

%unload data
partial_A = source_data_AB_remap.partial_shift.partial_A_common_far_input;
partial_B = source_data_AB_remap.partial_shift.partial_B_common_far_input;

%all neurons
[~,p_ks2,ks2stat] = kstest2(partial_A(:,2),partial_B(:,2));

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



%% Combine Figure 3 stat tables
t1 = repmat({' '},1,12);
blank_row = cell2table(t1);

%spreadsheet name
spreadsheet_name = 'statistics_summary.xlsx';

%sheet name
sheet_name = 'Figure 3';

%write to Excel spreadsheet
insert_table_rows(t_frac_friedman,spreadsheet_name,sheet_name,'overwritesheet')
insert_table_rows(t_frac_remap,spreadsheet_name,sheet_name,'append')

insert_table_rows(blank_row,spreadsheet_name,sheet_name,'append')
insert_table_rows(t_pf_dist_com,spreadsheet_name,sheet_name,'append')

insert_table_rows(blank_row,spreadsheet_name,sheet_name,'append')
insert_table_rows(t_zone_shift,spreadsheet_name,sheet_name,'append')

insert_table_rows(blank_row,spreadsheet_name,sheet_name,'append')
insert_table_rows(t_frac_zones_friedman,spreadsheet_name,sheet_name,'append')
insert_table_rows(t_frac_zone_remap,spreadsheet_name,sheet_name,'append')

insert_table_rows(blank_row,spreadsheet_name,sheet_name,'append')
insert_table_rows(t_2ks_partial_pf_dist,spreadsheet_name,sheet_name,'append')

%% 1-way Repeated Measures ANOVA - deal with this later (first figure - do with GG correction - works with MATLAB exchange function


%rm_sample_data = [];

% training_groups = repmat([1:4],4,1);
% training_groups = training_groups(:);
% 
% rm_table = table(rm_sample_data(:),training_groups,'VariableNames', {'score','t_grps'});
% 
% 
% rm_model = fitrm(rm_table,'score~t_grps','WithinModel','separatemeans')
% 
% rm_model = fitrm(rm_table,'score~t_grps','WithinDesign',table(rm_sample_data'))
% 
% rm_anova = ranova(rm_model)
% 
% anova_rm(rm_sample_data)

%repeated measures 1-way ANOVA
%paired t-test

%unpaired t-test
%kruskall wallis test
%write up functions for each test
%import all matlab data here

