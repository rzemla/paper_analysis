%% Load place cell structs for each animal

%not much obvious effect
s{1}{1} = load('G:\I57_LT_RF_bas_sal_792_031419\1\output\15_Mar_2019_ca_analysis.mat','Place_cell');
s{1}{3} = load('G:\I57_LT_RF_bas_sal_792_031419\3\output\15_Mar_2019_ca_analysis.mat','Place_cell');

% significant drop in place cells (116 --> 3) FOV1
s{2}{1} = load('G:\I57_RTLS_RF_bas_sal_792_031419\1\output\15_Mar_2019_ca_analysis.mat','Place_cell');
s{2}{3} = load('G:\I57_RTLS_RF_bas_sal_792_031419\3\output\15_Mar_2019_ca_analysis.mat','Place_cell');

% moderate drop in place cells 143 --> 56
s{3}{1} = load('G:\I56_RTLS_RF_bas_sal_792_031419\1\output\16_Mar_2019_ca_analysis.mat','Place_cell');
s{3}{3} = load('G:\I56_RTLS_RF_bas_sal_792_031419\3\output\18_Mar_2019_ca_analysis.mat','Place_cell');


%% Plot tuned STCs for each animal




%% Match neurons across days and plot their STCs



%% Plot absolute # of tuned neurons in each session
