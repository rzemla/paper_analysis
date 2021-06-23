function [outputArg1,outputArg2] = write_mean_sem(actx_word_p,doc_handle,mean_txt,sem_txt)
% set/reset default word style and font settings
%actx_word_p.Selection.Style = 'Normal';
actx_word_p.Selection.Font.Name = 'Arial';
actx_word_p.Selection.Font.Size = 12;
%set font color 
actx_word_p.Selection.Font.ColorIndex='wdAuto';%set back to default color

%write out mean
actx_word_p.Selection.TypeText(num2str(round(mean_txt,2)));

actx_word_p.Selection.TypeText(' ');
%insert +/- symbol for sem
actx_word_p.Selection.InsertSymbol(177)
actx_word_p.Selection.TypeText(' ');

%insert sem value
actx_word_p.Selection.TypeText(num2str(round(sem_txt,2)));


end

