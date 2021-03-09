function WordText(actx_word_p,text_p,style_p,enters_p,italic_p,sub_p)
	%VB Macro
	%Selection.TypeText Text:="Test!"
	%in Matlab
	%set(word.Selection,'Text','test');
	%this also works
	%word.Selection.TypeText('This is a test');    
    if(enters_p(1))
        actx_word_p.Selection.TypeParagraph; %enter
    end
	actx_word_p.Selection.Style = style_p;
    
    if (italic_p ==1)
        %actx_word_p.Selection.Font.ColorIndex=color_p; 
        actx_word_p.Selection.Font.Italic = true;
    else
        actx_word_p.Selection.Font.Italic = false;
    end
    
    if (sub_p ==1)
    actx_word_p.Selection.Font.Subscript = true;
    else
        actx_word_p.Selection.Font.Subscript = false;
    end 
    
	actx_word_p.Selection.TypeText(text_p);
    actx_word_p.Selection.Font.ColorIndex='wdAuto';%set back to default color
    
    %enter an enter after the end of tex
    for k=1:enters_p(2)
        actx_word_p.Selection.TypeParagraph; %enter
    end
return
