%% Import all table formatted entries here for data generated for each figure

%import data for each figure
laptop_access = 1;

%laptop path directory
laptop_path_dir = 'C:\Users\rzeml\Google Drive\task_selective_place_paper\matlab_data';
%desktop path directory
desktop_path_dir = 'G:\Google_drive\task_selective_place_paper\matlab_data';

%load table data
if laptop_access == 1
cd(laptop_path_dir)
else
    cd(desktop_path_dir)
end

%load in figure 2 data
fig2_data = load('fig2_table_data.mat');

%% Make Word export stats script for each figure legend

%% Start Word document that will contain the formatted stats data

WordFileName='legend_stats_formatted.doc';
CurDir=pwd;
FileSpec = fullfile(CurDir,WordFileName);
%active X handle for manipulating document (ActXWord)
[ActXWord,WordHandle]=StartWord(FileSpec);

fprintf('Document will be saved in %s\n',FileSpec);

%% Figure 2 Stats Write

%get input data for Friedman test:

dof = table2array(fig2_data.table_list.t_frac_si(1,7));
test_stat = table2array(fig2_data.table_list.t_frac_si(1,8));
p_val = table2array(fig2_data.table_list.t_frac_si(1,9));
sample_n = table2array(fig2_data.table_list.t_frac_si(1,5));
%description of comparison
comp_descrip = 'S.I. criterion';

%Fig. 2d output
writeFriedmanTest(ActXWord,WordHandle,comp_descrip,test_stat,p_val, dof, sample_n)


CloseWord(ActXWord,WordHandle,FileSpec);

%%
%write each statisics as text to this Word document

%input variables 
dof1, dof2, f_score, p_val,
writeWordFtest_1way
writeWordFtest_2way
writeWordText




%% 
%chi character code: 1D712 - unicode hex value;
%convert hex to dec: 120594

%order of entry
%test name, test statistic, p value (include rounding and sigstar add),
%test statistic - 2 decimal places
%p-value 3 decimal places and sub to <0.001 for low p values
%nb of samples as a char entry

%activeX.VBA control input into Word
Selection.Font.Name = "Arial"
Selection.Font.Size = 14
actx_word_p.Selection.TypeParagraph; %newline/enter command
actx_word_p.Selection.Font.Italic = true;
actx_word_p.Selection.Style
actx_word_p.Selection.Font.Subscript

actx_word_p.Selection.Font.ColorIndex='wdAuto'
%insertSymbol (I think this is the ascii limited form (no special chars
%beyond this)
actx_word_p.Selection.InsertSymbol(symbol_int_p);

%insert chi symbol macro commands in Word




Style='Normal';
TextString='This is a test ';
WordText(ActXWord,TextString,Style,[0,1],0,0);%enter after text

TextString='F';
WordText(ActXWord,TextString,Style,[0,0],1,0);%enter after text

TextString='(.506, 4.51)';
WordText(ActXWord,TextString,Style,[0,0],0,1);%enter after text

TextString=' = 12.44*';
WordText(ActXWord,TextString,Style,[0,0],0,0);%enter after text

%close the word document

close all;


%% Original function customized
WriteToWordFromMatlab_testing
