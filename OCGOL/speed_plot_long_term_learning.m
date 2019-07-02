%% File directories
%I45_RT
pathDir{1}{1} = 'G:\OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_rand_d1_052218';
pathDir{1}{2} = 'G:\OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_5A5B_053018';
pathDir{1}{3} = 'G:\OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_3A3B_060518';
pathDir{1}{4} = 'G:\OCGOL_learning_long_term\I45_RT\behavior_only\I45_RT_AB_061418';

%I46
pathDir{2}{1} = 'G:\OCGOL_learning_long_term\I46\behavior_only\I46_rand_d1_052918';
pathDir{2}{2} = 'G:\OCGOL_learning_long_term\I46\behavior_only\I46_5A5B_060118';
pathDir{2}{3} = 'G:\OCGOL_learning_long_term\I46\behavior_only\I46_3A3B_060718';
pathDir{2}{4} = 'G:\OCGOL_learning_long_term\I46\behavior_only\I46_AB_061518';

%I47_RS
pathDir{3}{1} = 'G:\OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_rand_d2_051518';
pathDir{3}{2} = 'G:\OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_5AB_d7_052218';
pathDir{3}{3} = 'G:\OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_3AB_d8_052418';
pathDir{3}{4} = 'G:\OCGOL_learning_long_term\I47_RS\behavior_only\I47_RS_AB_061418';

%I47_LP
pathDir{4}{1} = 'G:\OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_rand_d2_051518';
pathDir{4}{2} = 'G:\OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_5AB_d1_051718';
pathDir{4}{3} = 'G:\OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_3AB_d8_052418';
pathDir{4}{4} = 'G:\OCGOL_learning_long_term\I47_LP\behavior_only\I47_LP_AB_061418';

%% Read in velocity related data for each animal and construct related data

%for each animal
for aa=1:4
    %for each session
    for ss=1:4
        speed_data{aa}{ss} = load(fullfile(pathDir{aa}{ss},'cumul_analysis','speed_data.mat'));
    end
end

%extract lap idx's with A or B trials for each animal (2-4 sessions only)
for aa=1:4
    for ss=2:4
        %lap idx by trial type
        speed_data{aa}{ss}.lap_idx.A = find(speed_data{aa}{ss}.Behavior.lap_id.trial_based == 2);
        speed_data{aa}{ss}.lap_idx.B = find(speed_data{aa}{ss}.Behavior.lap_id.trial_based == 3);
    end
end

%% Plot

idx_width = [75,100];

figure
hold on
title('A zone speed')
for aa=1:4
    for ss=2:4
        for ll=1:size(speed_data{aa}{ss}.lap_idx.A,2)
            %plot line plot along range
            plot(speed_data{aa}{ss}.rewards.A.speed{ll}((speed_data{aa}{ss}.rewards.A.Imin(ll)-idx_width(1)):(speed_data{aa}{ss}.rewards.A.Imin(ll)+idx_width(2))));
            %plot line showing start of reward range
            plot([idx_width(1) idx_width(1)],[0 25],'k')
        end
    end
end




