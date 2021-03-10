function [outputArg1,outputArg2] = writeWordEnter(actx_word_p,doc_handle,nb_enter)
% set/reset default word style and font settings
%actx_word_p.Selection.Style = 'Normal';
actx_word_p.Selection.Font.Name = 'Arial';
actx_word_p.Selection.Font.Size = 12;
%set font color 
actx_word_p.Selection.Font.ColorIndex='wdAuto';%set back to default color

%insert nb_enter newlines
for ii=1:nb_enter
    actx_word_p.Selection.TypeParagraph;
end

end

