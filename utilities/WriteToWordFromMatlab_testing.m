function WriteToWordFromMatlab_testing
% -------------------------------------------------------------------
% File: WriteToWordFromMatlab
% Descr:  This is an example of how to control MS-Word from Matlab.
%         With the subfunctions below it is simple to automatically
%         let Matlab create reports in MS-Word.
%         This example copies two Matlab figures into MS-Word, writes
%         some text and inserts a table.
%         Works with MS-Word 2003 at least.
%         Reported to work with later versions of MS-Word as well.
%         MS-Word developer reference:
%         https://msdn.microsoft.com/en-us/library/office/ff837519.aspx
% Created: 2005-11-22 Andreas Karlsson
% History:
% 051122  AK  Modification of 'save2word' in Mathworks File Exchange   
% 060204  AK  Updated with WordSymbol, WordCreateTable and "Flying Start" section 
% 060214  AK  Pagenumber, font color and TOC added
% 161024  AK  v1.3 Changed calls to color-method into colorIndex to support MS-Word 2013. Thanks Martin van de Ven!
% -------------------------------------------------------------------

	WordFileName='TestDoc.doc';
	CurDir=pwd;
	FileSpec = fullfile(CurDir,WordFileName);
	[ActXWord,WordHandle]=StartWord(FileSpec);
    
    fprintf('Document will be saved in %s\n',FileSpec);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%% Setup an example report statistic for 2c
% 
% textString = '(fraction of licks in A trials on RF, 5A5B, 3A3B, Random AB; one-way RM ANOVA', ...
% effect of training stage, F1.506, 4.51 = 12.44, *P = 0.017, n = 16 training';
% sessions from 4 mice; RF vs. Random AB: two-tailed paired t-test, t3 = 5.74, 
% *P = 0.031, Holm-Sidak correction; 
    Style='Normal';
    TextString='This is a test ';
    WordText(ActXWord,TextString,Style,[0,1],0,0);%enter after text
    
    TextString='F';
    WordText(ActXWord,TextString,Style,[0,0],1,0);%enter after text

    TextString='(.506, 4.51)';
    WordText(ActXWord,TextString,Style,[0,0],0,1);%enter after text
    
    TextString=' = 12.44*';
    WordText(ActXWord,TextString,Style,[0,0],0,0);%enter after text
    
    
    CloseWord(ActXWord,WordHandle,FileSpec);    
    close all;
return    
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [actx_word,word_handle]=StartWord(word_file_p)
    % Start an ActiveX session with Word:
    actx_word = actxserver('Word.Application');
    actx_word.Visible = true;
    trace(actx_word.Visible);  
    if ~exist(word_file_p,'file');
        % Create new document:
        word_handle = invoke(actx_word.Documents,'Add');
    else
        % Open existing document:
        word_handle = invoke(actx_word.Documents,'Open',word_file_p);
    end           
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WordGoTo(actx_word_p,what_p,which_p,count_p,name_p,delete_p)
    %Selection.GoTo(What,Which,Count,Name)
    actx_word_p.Selection.GoTo(what_p,which_p,count_p,name_p);
    if(delete_p)
        actx_word_p.Selection.Delete;
    end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WordSymbol(actx_word_p,symbol_int_p)
    % symbol_int_p holds an integer representing a symbol, 
    % the integer can be found in MSWord's insert->symbol window    
    % 176 = degree symbol
    actx_word_p.Selection.InsertSymbol(symbol_int_p);
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WordCreateTable(actx_word_p,nr_rows_p,nr_cols_p,data_cell_p,enter_p) 
    %Add a table which auto fits cell's size to contents
    if(enter_p(1))
        actx_word_p.Selection.TypeParagraph; %enter
    end
    %create the table
    %Add = handle Add(handle, handle, int32, int32, Variant(Optional))
    actx_word_p.ActiveDocument.Tables.Add(actx_word_p.Selection.Range,nr_rows_p,nr_cols_p,1,1);
    %Hard-coded optionals                     
    %first 1 same as DefaultTableBehavior:=wdWord9TableBehavior
    %last  1 same as AutoFitBehavior:= wdAutoFitContent
     
    %write the data into the table
    for r=1:nr_rows_p
        for c=1:nr_cols_p
            %write data into current cell
            WordText(actx_word_p,data_cell_p{r,c},'Normal',[0,0]);
            
            if(r*c==nr_rows_p*nr_cols_p)
                %we are done, leave the table
                actx_word_p.Selection.MoveDown;
            else%move on to next cell 
                actx_word_p.Selection.MoveRight;
            end            
        end
    end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WordPageNumbers(actx_word_p,align_p)
    %make sure the window isn't split
    if (~strcmp(actx_word_p.ActiveWindow.View.SplitSpecial,'wdPaneNone')) 
        actx_word_p.Panes(2).Close;
    end
    %make sure we are in printview
    if (strcmp(actx_word_p.ActiveWindow.ActivePane.View.Type,'wdNormalView') | ...
        strcmp(actx_word_p.ActiveWindow.ActivePane.View.Type,'wdOutlineView'))
        actx_word_p.ActiveWindow.ActivePane.View.Type ='wdPrintView';
    end
    %view the headers-footers
    actx_word_p.ActiveWindow.ActivePane.View.SeekView='wdSeekCurrentPageHeader';
    if actx_word_p.Selection.HeaderFooter.IsHeader
        actx_word_p.ActiveWindow.ActivePane.View.SeekView='wdSeekCurrentPageFooter';
    else
        actx_word_p.ActiveWindow.ActivePane.View.SeekView='wdSeekCurrentPageHeader';
    end
    %now add the pagenumbers 0->don't display any pagenumber on first page
     actx_word_p.Selection.HeaderFooter.PageNumbers.Add(align_p,0);
     actx_word_p.ActiveWindow.ActivePane.View.SeekView='wdSeekMainDocument';
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PrintMethods(actx_word_p,category_p)
    style='Heading 3';
    text=strcat(category_p,'-methods');
    WordText(actx_word_p,text,style,[1,1]);           
    
    style='Normal';    
    text=strcat('Methods called from Matlab as: ActXWord.',category_p,'.MethodName(xxx)');
    WordText(actx_word_p,text,style,[0,0]);           
    text='Ignore the first parameter "handle". ';
    WordText(actx_word_p,text,style,[1,3]);           
    
    MethodsStruct=eval(['invoke(actx_word_p.' category_p ')']);
    MethodsCell=struct2cell(MethodsStruct);
    NrOfFcns=length(MethodsCell);
    for i=1:NrOfFcns
        MethodString=MethodsCell{i};
        WordText(actx_word_p,MethodString,style,[0,1]);           
    end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FigureIntoWord(actx_word_p)
	% Capture current figure/model into clipboard:
	print -dmeta
	% Find end of document and make it the insertion point:
	end_of_doc = get(actx_word_p.activedocument.content,'end');
	set(actx_word_p.application.selection,'Start',end_of_doc);
	set(actx_word_p.application.selection,'End',end_of_doc);
	% Paste the contents of the Clipboard:
    %also works Paste(ActXWord.Selection)
	invoke(actx_word_p.Selection,'Paste');
    actx_word_p.Selection.TypeParagraph; %enter    
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CloseWord(actx_word_p,word_handle_p,word_file_p)
    if ~exist(word_file_p,'file')
        % Save file as new:
        invoke(word_handle_p,'SaveAs',word_file_p,1);
    else
        % Save existing file:
        invoke(word_handle_p,'Save');
    end
    % Close the word window:
    invoke(word_handle_p,'Close');            
    % Quit MS Word
    invoke(actx_word_p,'Quit');            
    % Close Word and terminate ActiveX:
    delete(actx_word_p);            
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%