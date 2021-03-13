function [outputArg1,outputArg2] = write_1wayLME(actx_word_p,doc_handle,comp_descrip,test_stat,p_val, dof, sample_n)

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
actx_word_p.Selection.TypeText(': one-way RM mixed effects analysis, ');

%F-stat test description
actx_word_p.Selection.TypeText('effect of training day, ');

%insert test statistic W
actx_word_p.Selection.Font.Italic = true;
actx_word_p.Selection.TypeText('F')
actx_word_p.Selection.Font.Italic = false;

%extract dofs from test statistic input
split_tstat = split(test_stat,{'{','(',')','}'});

%get dofs and round to 2 decimal places
split_dofs = split(split_tstat(2),{', '});
%1st DOF
dof_1 = num2str(round(str2num(split_dofs(1)),2));
%2nd DOF
dof_2 = num2str(round(str2num(split_dofs(2)),2));

%with DOF as subscript
%set subscript state
actx_word_p.Selection.Font.Subscript = true;
actx_word_p.Selection.TypeText([dof_1,', ',dof_2])
actx_word_p.Selection.Font.Subscript = false;

actx_word_p.Selection.TypeText(' = ')

%isolate the test statistic
split_tstat_temp = split(split_tstat(3),' = ');
tstat_format = num2str(round(str2num(split_tstat_temp(2)),2));

%test statistic value (2 decimal round)
actx_word_p.Selection.TypeText(tstat_format);

%insert p value star significance
actx_word_p.Selection.TypeText(', '); 
%check if significant to insert significance star
if p_val < 0.05
actx_word_p.Selection.TypeText(char(check_p_value_sig(p_val)))
end

%turn italic font on
actx_word_p.Selection.Font.Italic = true;
actx_word_p.Selection.TypeText('P')
actx_word_p.Selection.Font.Italic = false;
%insert p-value depending on value
if p_val < 0.001
actx_word_p.Selection.TypeText(' < 0.001')
else
    actx_word_p.Selection.TypeText(' = ');
    actx_word_p.Selection.TypeText(num2str(round(p_val,3)));
end

%number of animals
actx_word_p.Selection.TypeText(', ');
actx_word_p.Selection.Font.Italic = true;
actx_word_p.Selection.TypeText('n');
actx_word_p.Selection.Font.Italic = false;
actx_word_p.Selection.TypeText([' = ', num2str(sample_n),' ', 'mice']);



end

