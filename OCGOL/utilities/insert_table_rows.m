function [outputArg1,outputArg2] = insert_table_rows(table_entry,spreadsheet_name, sheet_name,write_mode)

%write table rows to specified Excel spreadsheet and sheet 
writetable(table_entry,spreadsheet_name,'Sheet',sheet_name,'UseExcel', true,'WriteMode',write_mode)

end

