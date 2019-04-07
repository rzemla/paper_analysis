function [F_vars,CNMF_output] = loadImagingData_V1(directory_name)

%load imaging data from read_inputs directory in experiment directory

tic;
%load CNMF file
disp('Loading CNMF and dF/F data');
load(fullfile(directory_name,'read_inputs','CNMF_output.mat'),'CNMF_output');
%load F_vars file
load(fullfile(directory_name,'read_inputs','F_vars.mat'),'F_vars');
toc;

end

