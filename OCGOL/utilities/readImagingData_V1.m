function [CNMF_output] = readImagingData_V1(directory_name,options)

%% Set options
%select calcium imagting data type to read
calcium_data_type = options.calcium_data_input;

%% Get CNMF mat file name and full path

file_info = dir([directory_name,'\input\*.mat']);


%% Load in calcium data
%select which type of calcium signal to read in
switch calcium_data_type
    
    case 'CNMF'
        %load necessary CNMF variables
        CNMF_output = load(fullfile(file_info.folder, file_info.name),...
            'A_keep',... %spatial components after filtering for LQ componenets
            'b',...     %background componenet   
            'FOV',...    %average of imaging stack (x*y matrix)
            'C_full',... %temporal componenets extracted on native frame series
            'f_full',...
            'F_dark',...
            'T',...      %length of imaging series (without first frame) 
            'R_full',... %residual of matrix
            'options');  %options used to run CNMF
        

end

%% Save calcium imaging variables as workspace

%save the read_input directory of experiment directory
save(fullfile(directory_name,'read_inputs','CNMF_output.mat'),'CNMF_output');



end

