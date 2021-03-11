function [outputArg1,outputArg2] = write2ksTest(actx_word_p,doc_handle,comp_descrip,test_stat,p_val, dof, sample_n)

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
actx_word_p.Selection.TypeText(': 2-sample Kolmogorov-Smirnov test, ');

%insert test statistic...
actx_word_p.Selection.Font.Italic = true;
actx_word_p.Selection.TypeText('D')
actx_word_p.Selection.Font.Italic = false;

%insert subscript samples next to test statistic

%use sample for test statistic dof input
sample_split = split(sample_n, ' ');
sample_formatted = strjoin(sample_split([1,3]),', ');
%turn subscript on
actx_word_p.Selection.Font.Subscript = true;
actx_word_p.Selection.TypeText(sample_formatted);
actx_word_p.Selection.Font.Subscript = false;

actx_word_p.Selection.TypeText(' = ');

%test statistic value (2 decimal round)
actx_word_p.Selection.TypeText(num2str(round(test_stat,2)));
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


%overwrite default n sample here (used before to display number of neurons
%join with vs. 'delimiter'
sample_n_neurons = strjoin(sample_split([1,3]),' vs. ');

sample_n = 11; %set fixed value for now

%number of animals
actx_word_p.Selection.TypeText(', ');
actx_word_p.Selection.Font.Italic = true;
actx_word_p.Selection.TypeText('n');
actx_word_p.Selection.Font.Italic = false;
actx_word_p.Selection.TypeText([' = ', num2str(sample_n_neurons),...
        ' neurons, ', num2str(sample_n),' FOV from ', num2str(sample_n-1),' mice']);


%activeX unit conversion
%https://www.mathworks.com/matlabcentral/answers/165541-how-can-i-set-the-movedown-method-to-move-down-paragraphs-headings


end

