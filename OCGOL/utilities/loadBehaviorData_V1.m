function [CSV,XML] = loadBehaviorData_V1(directory_name)

%load behavioral data from read_inputs directory in experiment directory

tic;
disp('Loading behavior and imaging timestamp data...');

%load xml file
load(fullfile(directory_name,'read_inputs','XML.mat'),'XML');
%load csv file
load(fullfile(directory_name,'read_inputs','CSV.mat'),'CSV');

disp('Done');
toc;

end

