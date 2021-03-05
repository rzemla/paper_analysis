%% Export stats for all supplementary figures in paper

%import supplementary data

%figure 2/3 sup data
%load('G:\Google_drive\task_selective_place_paper\matlab_data\source_data_sup_2_3.mat');

%figure 4/5 sup data
%load('G:\Google_drive\task_selective_place_paper\matlab_data\source_data_fig4_5_and_sup.mat');

%ex figure 10 sup data
%load('G:\Google_drive\task_selective_place_paper\matlab_data\source_data_ex10_sup.mat');

%laptop file directory
load('C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data\source_data_sup_2_3.mat');
load('C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data\source_data_fig4_5_and_sup.mat');
load('C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data\source_data_ex10_sup.mat');

%% Figure 4 A-B speed difference for A- and B- selective neurons

%speed difference for A and B- selective neurons
speed_diff_Asel = source_data_sup_2_3.event_speed_plot.Asel_speed_diff;
speed_diff_Bsel = source_data_sup_2_3.event_speed_plot.Bsel_speed_diff;

%remove any nan values from each vector
%A
nan_ct = find(isnan(speed_diff_Asel)== 1);
 if ~isempty(nan_ct)
     speed_diff_Asel(nan_ct) = [];
 end
%B
nan_ct = find(isnan(speed_diff_Bsel)== 1);
 if ~isempty(nan_ct)
     speed_diff_Bsel(nan_ct) = [];
 end
 
%for 1 sample Wilcoxon tests
paired_wilcoxon_stats_Asel = paired_wilcoxon_signrank(speed_diff_Asel,0);

%test comparison description for each table entry \
comp_descrip_in = {'A-B speed difference for A selective neurons'};

%table entries for paired wilcoxon test learning
%A selective
t_paired_wilcox_Asel = paired_wilcoxon_table_entry_no_adj(4,'b','pooled',...
                comp_descrip_in,paired_wilcoxon_stats_Asel);

%for 1 sample Wilcoxon tests
paired_wilcoxon_stats_Bsel = paired_wilcoxon_signrank(speed_diff_Bsel,0);

%test comparison description for each table entry \
comp_descrip_in = {'A-B speed difference for B selective neurons'};            
            
%B selective
t_paired_wilcox_Bsel = paired_wilcoxon_table_entry_no_adj(4,'b','pooled',...
                comp_descrip_in,paired_wilcoxon_stats_Bsel);
            
%% Figure 4 Place field count differences between task selective place cells

%normalized counts
Asel_norm = source_data_sup_2_3.pf_count.Asel_norm;
Bsel_norm = source_data_sup_2_3.pf_count.Bsel_norm;

AB.A_norm = source_data_sup_2_3.pf_count.AB.A_norm;
AB.B_norm = source_data_sup_2_3.pf_count.AB.B_norm;

%assemble into single field matrices (animal x A,B,AB-A, AB-B
single_field = [Asel_norm(:,1), Bsel_norm(:,1), AB.A_norm(:,1), AB.B_norm(:,1)];
double_field = [Asel_norm(:,2), Bsel_norm(:,2), AB.A_norm(:,2), AB.B_norm(:,2)];
triple_field = [Asel_norm(:,3), Bsel_norm(:,3), AB.A_norm(:,3), AB.B_norm(:,3)];

%single paired tests
%A vs B
paired_wilcoxon_stats.single.AvsB = paired_wilcoxon_signrank(single_field(:,1),single_field(:,2));
%A vs AB-A
paired_wilcoxon_stats.single.AvsAB_A = paired_wilcoxon_signrank(single_field(:,1),single_field(:,3));
%B vs AB-B
paired_wilcoxon_stats.single.BvsAB_B = paired_wilcoxon_signrank(single_field(:,2),single_field(:,4));

%double paired tests
paired_wilcoxon_stats.double.AvsB = paired_wilcoxon_signrank(double_field(:,1),double_field(:,2));
%A vs AB-A
paired_wilcoxon_stats.double.AvsAB_A = paired_wilcoxon_signrank(double_field(:,1),double_field(:,3));
%B vs AB-B
paired_wilcoxon_stats.double.BvsAB_B = paired_wilcoxon_signrank(double_field(:,2),double_field(:,4));

%triple paired tests
paired_wilcoxon_stats.triple.AvsB = paired_wilcoxon_signrank(triple_field(:,1),triple_field(:,2));
%A vs AB-A
paired_wilcoxon_stats.triple.AvsAB_A = paired_wilcoxon_signrank(triple_field(:,1),triple_field(:,3));
%B vs AB-B
paired_wilcoxon_stats.triple.BvsAB_B = paired_wilcoxon_signrank(triple_field(:,2),triple_field(:,4));

%single field stats combined
stats_in_single = [paired_wilcoxon_stats.single.AvsB; paired_wilcoxon_stats.single.AvsAB_A; paired_wilcoxon_stats.single.BvsAB_B];

comp_descrip = {'A vs B (1 vs 2) single field';'A vs A&B A (1 vs 3) single field';'B vs A&B (2 vs 4) single field'};

%do the adjustment in the function  for multo comp (3-way Holm-Sidak
%correction)
[t_out.wilcox_singlePF ] = paired_wilcoxon_table_entry(4,'c','by animal',...
                comp_descrip, stats_in_single);

%double field stats combined
stats_in_double = [paired_wilcoxon_stats.double.AvsB; paired_wilcoxon_stats.double.AvsAB_A; paired_wilcoxon_stats.double.BvsAB_B];

comp_descrip = {'A vs B (1 vs 2) double field';'A vs A&B A (1 vs 3) double field';'B vs A&B (2 vs 4) double field'};

%do the adjustment in the function  for multo comp (3-way Holm-Sidak
%correction)
[t_out.wilcox_doublePF ] = paired_wilcoxon_table_entry(4,'c','by animal',...
                comp_descrip, stats_in_double);

%triple field stats combined
stats_in_triple = [paired_wilcoxon_stats.triple.AvsB; paired_wilcoxon_stats.triple.AvsAB_A; paired_wilcoxon_stats.triple.BvsAB_B];

comp_descrip = {'A vs B (1 vs 2) triple field';'A vs A&B A (1 vs 3) triple field';'B vs A&B (2 vs 4) triple field'};

%do the adjustment in the function  for multo comp (3-way Holm-Sidak
%correction)
[t_out.wilcox_triplePF ] = paired_wilcoxon_table_entry(4,'c','by animal',...
                comp_descrip, stats_in_triple);
%% Figure 4 Place field width differences between task selective place cells
%2-sample KS test

%unload the data
width_Asel = source_data_sup_2_3.pf_width_pool.Asel;
width_Bsel = source_data_sup_2_3.pf_width_pool.Bsel;
width_ABA = source_data_sup_2_3.pf_width_pool.AB.A;
width_ABB = source_data_sup_2_3.pf_width_pool.AB.B;

%run all the KS tests

%Asel vs Bsel
kstest2_AB = kstest2_stats(width_Asel,width_Bsel);

%Asel vs ABAsel
kstest2_A_ABA = kstest2_stats(width_Asel,width_ABA);

%Asel vs ABBsel
kstest2_A_ABB = kstest2_stats(width_Asel,width_ABB);

%Bsel vs ABAsel
kstest2_B_ABA = kstest2_stats(width_Bsel,width_ABA);

%Bsel vs ABBsel
kstest2_B_ABB = kstest2_stats(width_Bsel,width_ABB);

%ABAsel vs ABBsel
kstest2_ABA_ABB = kstest2_stats(width_ABA,width_ABB);

%KS test mult compare correction and table entry generation
%setup data input for KS2 table function
input_data = [kstest2_AB; kstest2_A_ABA; kstest2_A_ABB;...
            kstest2_B_ABA; kstest2_B_ABB; kstest2_ABA_ABB];

kstest2_mult_compare(input_data)


        
%6-way Holm-Sidak (Holm-Bonferroni used previously)

%% Figure 7 Activity index score difference and activity rate diff. for activity remapping neurons


%% Figure 9 Eqivalent long term correlation data (Fig 4) using SI criterion

%% Figure 10 Long term recall data SI and TS criteria

%% Create Ex. Figure 4 stats export spreadsheet
%spreadsheet name
spreadsheet_name = 'statistics_summary.xlsx';

%sheet name
sheet_name = 'Extended Figure 4';

%empty row
t1 = repmat({' '},1,12);

%exported Excel spreadsheet
%write to Excel spreadsheet
%Ex Fig 4b
writetable(t_paired_wilcox_Asel,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','overwritesheet')
writetable(t_paired_wilcox_Bsel,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 4c - place field count between selective neurons
writetable(cell2table(t1),spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_out.wilcox_singlePF,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_out.wilcox_doublePF,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')
writetable(t_out.wilcox_triplePF,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode','append')

%Ex Fig 4c - place field width between task selective neurons


