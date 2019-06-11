function [tuned_count, fracTuned_si, fracTuned_ts] = fractionTuned(tunedLogical)

%% Define # of tuned neurons in each catergory

%spatial information
tuned_count{1}(1) = size(find(tunedLogical.si.onlyA_tuned ==1),2);
tuned_count{1}(2) = size(find(tunedLogical.si.onlyB_tuned ==1),2);
tuned_count{1}(3) = size(find(tunedLogical.si.AandB_tuned ==1),2);
tuned_count{1}(4) = size(find(tunedLogical.si.neither ==1),2); 

%tuning specificity
tuned_count{2}(1) = size(find(tunedLogical.ts.onlyA_tuned ==1),2);
tuned_count{2}(2) = size(find(tunedLogical.ts.onlyB_tuned ==1),2);
tuned_count{2}(3) = size(find(tunedLogical.ts.AandB_tuned ==1),2);
tuned_count{2}(4) = size(find(tunedLogical.ts.neither ==1),2); 

%% Calculate fraction tuned
%si
fracTuned_si = tuned_count{1}/sum(tuned_count{1});

%ts
fracTuned_ts = tuned_count{2}/sum(tuned_count{2});

%% Plot pie chart for each type of tuning criterion

figure('Position',[2050 520 1140 450])
subplot(1,2,1)
p = pie(fracTuned_si,{['A ', num2str(round(100*fracTuned_si(1))), '%'],...
                        ['B ', num2str(round(100*fracTuned_si(2))), '%'],...
                        ['A&B ', num2str(round(100*fracTuned_si(3))), '%'],...
                        ['     Neither ', num2str(round(100*fracTuned_si(4))), '%']});
hold on
title('Percentage of active neurons tuned (S.I.)');
colormap([0 0 1; 1 0 0; 1 0 1; 0.5 0.5 0.5])

subplot(1,2,2)
p = pie(fracTuned_ts,{['A ', num2str(round(100*fracTuned_ts(1))), '%'],...
                        ['B ', num2str(round(100*fracTuned_ts(2))), '%'],...
                        ['A&B ', num2str(round(100*fracTuned_ts(3))), '%'],...
                        ['     Neither ', num2str(round(100*fracTuned_ts(4))), '%']});
hold on
title('Percentage of active neurons tuned (T.S.)');
colormap([0 0 1; 1 0 0; 1 0 1; 0.5 0.5 0.5])

end

