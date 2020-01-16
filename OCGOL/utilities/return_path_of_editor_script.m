function [path_x] = return_path_of_editor_script
% returns the path and name of the editor script that is currently open 

path_x = matlab.desktop.editor.getActiveFilename;
fprintf('%s\n',path_x);

end

