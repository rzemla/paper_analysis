%% Read in raw fluorescense from Fiji selected ROIs

%read in csv with mean raw F values as table
dirRead = uigetdir;

filename = dir([dirRead, '\','*.csv']);
F_table = readtable(fullfile(dirRead,filename.name ));

%get mean F traces for each ROI
%get the column indices of the mean F values for each ROI 
F_idx = [3:4:size(F_table,2)];

%convert to a matrix
F = table2array(F_table(:,F_idx));


%% Plot the raw F traces for each ROI
for ii =1:size(F,2)
ROI = ii;
figure; hold on; plot(F(:,ROI)-6000);

end