function [CSV,XML] = readBehaviorData_V1(directory_name)
% read in and save CSV/XML data

%% Get the filenames of the .csv and .xml files

%get name of csv and xml file 
csv_file_nam = dir([directory_name,'\input\*.csv']);
xml_file_nam = dir([directory_name,'\input\*.xml']);

%% Read CSV with behavior voltage data
disp(['Reading CSV']);
%read voltage recording data into matrix
%each column is a single channel
%first column is time
tic;
CSV = csvread(fullfile(csv_file_nam.folder,csv_file_nam.name),1,0);
toc

%% Read in XML into struct with frame timestamps from imaging
disp(['Reading XML']);
%read XML files (imaging frame timestamps)
tic;
XML = xml2structV2((fullfile(xml_file_nam.folder,xml_file_nam.name)));
toc;

%% Save the CSV and XML as mat files

%get main experiment directory
% split_path = split(directory_name, '\');
% exp_dir  = fullfile(split_path{1:(end-1)});

%try to make directory for saving CSV matrix and XML struct
try
    mkdir(fullfile(directory_name,'read_inputs'));
catch
    disp('Directory already exists');
end

%save XML
disp('Saving XML file')
tic;
save(fullfile(directory_name,'read_inputs','XML.mat'),'XML','-v7.3')
toc;

%save CSV
disp('Saving CSV file')
tic;
save(fullfile(directory_name,'read_inputs','CSV.mat'),'CSV','-v7.3')
toc;

end

