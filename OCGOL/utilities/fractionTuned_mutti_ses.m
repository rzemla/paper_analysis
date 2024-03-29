function [tuned_fractions] = fractionTuned(tunedLogical,pf_count_filtered_log)

%% Define # of tuned neurons in each catergory

%incoporate place field present and event filter
%SI
tuned_count_filt{1}(1) = size(find((tunedLogical.si.onlyA_tuned & pf_count_filtered_log(1,:)) ==1),2);
tuned_count_filt{1}(2) = size(find((tunedLogical.si.onlyB_tuned & pf_count_filtered_log(2,:)) ==1),2);
tuned_count_filt{1}(3) = size(find((tunedLogical.si.AandB_tuned & (pf_count_filtered_log(1,:) & pf_count_filtered_log(2,:))) ==1),2);
tuned_count_filt{1}(4) = size(pf_count_filtered_log,2) -  sum(tuned_count_filt{1});

%TS
tuned_count_filt{2}(1) = size(find((tunedLogical.ts.onlyA_tuned & pf_count_filtered_log(1,:)) ==1),2);
tuned_count_filt{2}(2) = size(find((tunedLogical.ts.onlyB_tuned & pf_count_filtered_log(2,:)) ==1),2);
tuned_count_filt{2}(3) = size(find((tunedLogical.ts.AandB_tuned & (pf_count_filtered_log(1,:) & pf_count_filtered_log(2,:))) ==1),2);
tuned_count_filt{2}(4) = size(pf_count_filtered_log,2) -  sum(tuned_count_filt{2});

%spatial information (not filtered aside from stat sig outcome)
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
fracTuned_si_filt = tuned_count_filt{1}/sum(tuned_count_filt{1});
%ts
fracTuned_ts = tuned_count{2}/sum(tuned_count{2});
fracTuned_ts_filt = tuned_count_filt{2}/sum(tuned_count_filt{2});

%% Plot pie chart for each type of tuning criterion (no filter)

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

%% Plot pie chart for each type of tuning criterion (filtered)

figure('Position',[2050 520 1140 450])
subplot(1,2,1)
p = pie(fracTuned_si_filt,{['A ', num2str(round(100*fracTuned_si_filt(1))), '%'],...
                        ['B ', num2str(round(100*fracTuned_si_filt(2))), '%'],...
                        ['A&B ', num2str(round(100*fracTuned_si_filt(3))), '%'],...
                        ['     Neither ', num2str(round(100*fracTuned_si_filt(4))), '%']});
hold on
title('Percentage of active neurons tuned (S.I.) - Filtered');
colormap([0 0 1; 1 0 0; 1 0 1; 0.5 0.5 0.5])

subplot(1,2,2)
p = pie(fracTuned_ts_filt,{['A ', num2str(round(100*fracTuned_ts_filt(1))), '%'],...
                        ['B ', num2str(round(100*fracTuned_ts_filt(2))), '%'],...
                        ['A&B ', num2str(round(100*fracTuned_ts_filt(3))), '%'],...
                        ['     Neither ', num2str(round(100*fracTuned_ts_filt(4))), '%']});
hold on
title('Percentage of active neurons tuned (T.S.) - Filtered');
colormap([0 0 1; 1 0 0; 1 0 1; 0.5 0.5 0.5])

%% Place into struct for export

%absolute counts for each subcategory
tuned_fractions.tuned_count = tuned_count;
tuned_fractions.tuned_count_filt = tuned_count_filt;
%fraction tuned by SI criterion
tuned_fractions.fracTuned_si = fracTuned_si;
tuned_fractions.fracTuned_si_filt = fracTuned_si_filt;
%fraction tuned by TS criterion
tuned_fractions.fracTuned_ts = fracTuned_ts;
tuned_fractions.fracTuned_ts_filt = fracTuned_ts_filt;

end

