function [outputArg1,outputArg2] = writeDefaultWordText(actx_word_p,doc_handle,txt)
% set/reset default word style and font settings
%actx_word_p.Selection.Style = 'Normal';
actx_word_p.Selection.Font.Name = 'Arial';
actx_word_p.Selection.Font.Size = 12;
%set font color 
actx_word_p.Selection.Font.ColorIndex='wdAuto';%set back to default color

%write out input text
actx_word_p.Selection.TypeText(txt);

end

