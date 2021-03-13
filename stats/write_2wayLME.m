function [outputArg1,outputArg2] = write_2wayLME(actx_word_p,doc_handle,comp_descrip,test_stat,p_val, dof, sample_n)

%word style and font settings
%actx_word_p.Selection.Style = 'Normal';
actx_word_p.Selection.Font.Name = 'Arial';
actx_word_p.Selection.Font.Size = 12;
%set font color 
actx_word_p.Selection.Font.ColorIndex='wdAuto';%set back to default color

%set subscript state
actx_word_p.Selection.Font.Subscript = false;

%% Descriptive test in bold
%set bold font
actx_word_p.Selection.Font.Bold = true;
%desriptor
actx_word_p.Selection.TypeText(comp_descrip)
%turn bold off
actx_word_p.Selection.Font.Bold = false;
%write test name
actx_word_p.Selection.TypeText(': two-way RM mixed effects analysis, ');

%process the data for the F-test statistics here

%extract dofs from test statistic input
%& char(10) & - Excel same cell newline separator
split_Fstat = split(test_stat,{'& char(10) &','=', '"',});
%remove whitespaces (1,2 - dof and test stat, 3,4 - dof and test stat, 5,6
%- dof and test stat for each factor/interaction
split_Fstat = split_Fstat([3,4,7,8,11,12]);

%create empty struct
split_dofs = {};

%get first time effect and round to 2 decimal places
idx_dof = [1,3,5];

%split the dofs for each F test
for ii=1:3
    split_dofs{ii} = split(split_Fstat(idx_dof(ii)),{', ','(',')'});
end

%extract both dofs and covert to numbers for each test
for ii=1:3
    %1st DOF of first F statistic
    dof_in{ii,1} = num2str(round(str2num(split_dofs{ii}(2)),2));
    %2nd DOF
    dof_in{ii,2} = num2str(round(str2num(split_dofs{ii}(3)),2));
end

%isolate the test statistic for each F-test
tstat_format = round(double(split_Fstat([2,4,6])),2);

%split each p-value associated with each F statistic
split_pval = split(p_val,{'& char(10) &','=', '"',});
split_pval = double(split_pval([3,6,9]));

%write out the formatted stats to Word

%F-stat test description (1st)
actx_word_p.Selection.TypeText('effect of time, ');

%insert test statistic
actx_word_p.Selection.Font.Italic = true;
actx_word_p.Selection.TypeText('F')
actx_word_p.Selection.Font.Italic = false;

%with DOF as subscript
%set subscript state
actx_word_p.Selection.Font.Subscript = true;
actx_word_p.Selection.TypeText([dof_in{1,1},', ',dof_in{1,2}])
actx_word_p.Selection.Font.Subscript = false;
actx_word_p.Selection.TypeText(' = ')

%test statistic value (2 decimal round) - 1st F
actx_word_p.Selection.TypeText(num2str(tstat_format(1)));

%insert p value star significance (first F value)
actx_word_p.Selection.TypeText(', '); 
%check if significant to insert significance star
if split_pval(1) < 0.05
actx_word_p.Selection.TypeText(char(check_p_value_sig(split_pval(1))))
end

%turn italic font on
actx_word_p.Selection.Font.Italic = true;
actx_word_p.Selection.TypeText('P')
actx_word_p.Selection.Font.Italic = false;
%insert p-value depending on value
if split_pval(1) < 0.001
actx_word_p.Selection.TypeText(' < 0.001')
else
    actx_word_p.Selection.TypeText(' = ');
    actx_word_p.Selection.TypeText(num2str(round(split_pval(1),3)));
end

%2nd F test stats
%F-stat test description (1st)
actx_word_p.Selection.TypeText(', ');
actx_word_p.Selection.TypeText('effect of behavior, ');

%insert test statistic
actx_word_p.Selection.Font.Italic = true;
actx_word_p.Selection.TypeText('F')
actx_word_p.Selection.Font.Italic = false;

%with DOF as subscript
%set subscript state
actx_word_p.Selection.Font.Subscript = true;
actx_word_p.Selection.TypeText([dof_in{2,1},', ',dof_in{2,2}])
actx_word_p.Selection.Font.Subscript = false;
actx_word_p.Selection.TypeText(' = ')

%test statistic value (2 decimal round) - 1st F
actx_word_p.Selection.TypeText(num2str(tstat_format(2)));

%insert p value star significance (first F value)
actx_word_p.Selection.TypeText(', '); 
%check if significant to insert significance star
if split_pval(2) < 0.05
actx_word_p.Selection.TypeText(char(check_p_value_sig(split_pval(2))))
end

%turn italic font on
actx_word_p.Selection.Font.Italic = true;
actx_word_p.Selection.TypeText('P')
actx_word_p.Selection.Font.Italic = false;
%insert p-value depending on value
if split_pval(2) < 0.001
actx_word_p.Selection.TypeText(' < 0.001')
else
    actx_word_p.Selection.TypeText(' = ');
    actx_word_p.Selection.TypeText(num2str(round(split_pval(2),3)));
end

%3nd F test stats
%F-stat test description (3rd)
actx_word_p.Selection.TypeText(', ');
actx_word_p.Selection.TypeText('interaction between time and behavior, ');

%insert test statistic
actx_word_p.Selection.Font.Italic = true;
actx_word_p.Selection.TypeText('F')
actx_word_p.Selection.Font.Italic = false;

%with DOF as subscript
%set subscript state
actx_word_p.Selection.Font.Subscript = true;
actx_word_p.Selection.TypeText([dof_in{3,1},', ',dof_in{3,2}])
actx_word_p.Selection.Font.Subscript = false;
actx_word_p.Selection.TypeText(' = ')

%test statistic value (2 decimal round) - 
actx_word_p.Selection.TypeText(num2str(tstat_format(3)));

%insert p value star significance 
actx_word_p.Selection.TypeText(', '); 
%check if significant to insert significance star
if split_pval(3) < 0.05
actx_word_p.Selection.TypeText(char(check_p_value_sig(split_pval(3))))
end

%turn italic font on
actx_word_p.Selection.Font.Italic = true;
actx_word_p.Selection.TypeText('P')
actx_word_p.Selection.Font.Italic = false;
%insert p-value depending on value
if split_pval(3) < 0.001
actx_word_p.Selection.TypeText(' < 0.001')
else
    actx_word_p.Selection.TypeText(' = ');
    actx_word_p.Selection.TypeText(num2str(round(split_pval(3),3)));
end

%process number of animals input
split_n = double(split(sample_n,', '));

%number of animals
actx_word_p.Selection.TypeText(', ');
actx_word_p.Selection.Font.Italic = true;
actx_word_p.Selection.TypeText('n');
actx_word_p.Selection.Font.Italic = false;
actx_word_p.Selection.TypeText([' = ', num2str(split_n(1)),' learn cohort, ',...
            num2str(split_n(2)),' recall cohort',' mice']);



end

