function [outputArg1,outputArg2] = writeFriedmanTest(actx_word_p,doc_handle,test_stat,p_val, dof, sample_n)

%word style and font settings
actx_word_p.Selection.Style = 'Normal';
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
actx_word_p.Selection.TypeText('S.I. criterion')
%turn bold off
actx_word_p.Selection.Font.Bold = false;
%write test name
actx_word_p.Selection.TypeText(': Friedman test, ');

%insert chi symbol as part of a equation

%insert equation
doc_handle.OMaths.Add(actx_word_p.Selection.Range);
eqn = doc_handle.OMaths.Item(1);
%insert char sub and superscript template
eqn.Functions.Add(actx_word_p.Selection.Range,'wdOMathFunctionScrSubSup').ScrSubSup.AlignScripts = false;
%shift insertion box 3 moves to the left (first val is setting for
%character mode)
actx_word_p.Selection.MoveLeft(1, 3)
%insert chi symbol (need both code insertions, true - means it used Unicode encoding for the symbol)
actx_word_p.Selection.InsertSymbol(-10187,'Cambria Math',true)
actx_word_p.Selection.InsertSymbol(-8430,'Cambria Math',true)
%move to the right
actx_word_p.Selection.MoveRight(1, 1);
%insert dof
actx_word_p.Selection.TypeText(num2str(dof))
%move up (5 is move up a line)
actx_word_p.Selection.MoveUp(5, 2);
actx_word_p.Selection.TypeText('2')
%done editting equation
eqn.BuildUp();
%move to the right (twice - get out of equation entry mode)
actx_word_p.Selection.MoveRight(1, 2);
actx_word_p.Selection.TypeText('2')

%activeX unit conversion
%https://www.mathworks.com/matlabcentral/answers/165541-how-can-i-set-the-movedown-method-to-move-down-paragraphs-headings


end

