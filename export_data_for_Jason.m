%% Directory containing data

data_dir{1} = 'G:\OCGOL_learning_short_term\I56_RTLS\crossSession';
data_dir{2} = 'G:\OCGOL_learning_short_term\I57_RTLS\crossSession';
data_dir{3} = 'G:\OCGOL_learning_short_term\I57_LT\crossSession';
data_dir{4} = 'E:\OCGOL_learning_short_term\I58_RT\crossSession';
data_dir{5} = 'E:\OCGOL_learning_short_term\I58_LT\crossSession';
data_dir{6} = 'E:\OCGOL_learning_short_term\I58_RTLP\crossSession';

%% Create directory with each dataset

%collective export directory
coll_dir = 'G:\Google_drive\task_selective_place_paper\export_data_for_Jason';

%make directories corresponding to animals listed above
cd('G:\Google_drive\task_selective_place_paper\export_data_for_Jason')

for ii=1:6
    mkdir(num2str(ii)) 
end

for ii=1:6
    %copy tuning curve data
    copyfile(fullfile(data_dir{ii},'matching_tun_curves.mat'), fullfile(coll_dir,num2str(ii)))
    %copy performance data
    copyfile(fullfile(data_dir{ii},'perf_lap_tables.mat'), fullfile(coll_dir,num2str(ii)))
end


