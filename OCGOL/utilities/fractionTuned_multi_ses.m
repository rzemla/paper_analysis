function [tuned_fractions,tuned_logicals] = fractionTuned_multi_ses(tunedLogical,pf_count_filtered_log,options)


%% Define parameters
sessionSelect = options.sessionSelect;
selectTrial = options.selectTrial;


%% Define # of tuned neurons in each catergory

%for each session
for ss=sessionSelect
    
    %incoporate place field present and event filter
    %SI
    tuned_count_filt{ss}{1}(1) = size(find((tunedLogical(ss).si.onlyA_tuned & pf_count_filtered_log{ss}(selectTrial(1),:)) ==1),2);
    tuned_count_filt{ss}{1}(2) = size(find((tunedLogical(ss).si.onlyB_tuned & pf_count_filtered_log{ss}(selectTrial(2),:)) ==1),2);
    tuned_count_filt{ss}{1}(3) = size(find((tunedLogical(ss).si.AandB_tuned & (pf_count_filtered_log{ss}(selectTrial(1),:) & pf_count_filtered_log{ss}(selectTrial(2),:))) ==1),2);
    tuned_count_filt{ss}{1}(4) = size(pf_count_filtered_log{ss},2) -  sum(tuned_count_filt{ss}{1}(1:3));
    
    %TS
    tuned_count_filt{ss}{2}(1) = size(find((tunedLogical(ss).ts.onlyA_tuned & pf_count_filtered_log{ss}(selectTrial(1),:)) ==1),2);
    tuned_count_filt{ss}{2}(2) = size(find((tunedLogical(ss).ts.onlyB_tuned & pf_count_filtered_log{ss}(selectTrial(2),:)) ==1),2);
    tuned_count_filt{ss}{2}(3) = size(find((tunedLogical(ss).ts.AandB_tuned & (pf_count_filtered_log{ss}(selectTrial(1),:) & pf_count_filtered_log{ss}(selectTrial(2),:))) ==1),2);
    tuned_count_filt{ss}{2}(4) = size(pf_count_filtered_log{ss},2) -  sum(tuned_count_filt{ss}{2}(1:3));
    
end

%% Generate logicals corresponding to A/B/AB/neither classes of neurons

for ss=sessionSelect
    
    %incoporate place field present and event filter
    %SI
    %only either category
    tuned_log_filt_si{ss}.onlyA = (tunedLogical(ss).si.onlyA_tuned & pf_count_filtered_log{ss}(selectTrial(1),:));
    tuned_log_filt_si{ss}.onlyB = (tunedLogical(ss).si.onlyB_tuned & pf_count_filtered_log{ss}(selectTrial(2),:)) ;
    
    %all in either category
    tuned_log_filt_si{ss}.allA = (tunedLogical(ss).si.Atuned & pf_count_filtered_log{ss}(selectTrial(1),:));
    tuned_log_filt_si{ss}.allB = (tunedLogical(ss).si.Btuned & pf_count_filtered_log{ss}(selectTrial(2),:)) ;
        
    
    tuned_log_filt_si{ss}.AB = (tunedLogical(ss).si.AandB_tuned & (pf_count_filtered_log{ss}(selectTrial(1),:) & pf_count_filtered_log{ss}(selectTrial(2),:)));
    tuned_log_filt_si{ss}.neither = ~((tuned_log_filt_si{ss}.onlyA | tuned_log_filt_si{ss}.onlyB) | tuned_log_filt_si{ss}.AB);
    
    %TS
    %only either category
    tuned_log_filt_ts{ss}.onlyA = (tunedLogical(ss).ts.onlyA_tuned & pf_count_filtered_log{ss}(selectTrial(1),:));
    tuned_log_filt_ts{ss}.onlyB = (tunedLogical(ss).ts.onlyB_tuned & pf_count_filtered_log{ss}(selectTrial(2),:)) ;
    
    %all in either category
    tuned_log_filt_ts{ss}.allA = (tunedLogical(ss).ts.Atuned & pf_count_filtered_log{ss}(selectTrial(1),:));
    tuned_log_filt_ts{ss}.allB = (tunedLogical(ss).ts.Btuned & pf_count_filtered_log{ss}(selectTrial(2),:)) ;
    
    
    tuned_log_filt_ts{ss}.AB = (tunedLogical(ss).ts.AandB_tuned & (pf_count_filtered_log{ss}(selectTrial(1),:) & pf_count_filtered_log{ss}(selectTrial(2),:)));
    tuned_log_filt_ts{ss}.neither = ~((tuned_log_filt_ts{ss}.onlyA | tuned_log_filt_ts{ss}.onlyB) | tuned_log_filt_ts{ss}.AB);
    
end

%% 
%code for non-event filtered data

% %spatial information (not filtered aside from stat sig outcome)
% tuned_count{1}(1) = size(find(tunedLogical.si.onlyA_tuned ==1),2);
% tuned_count{1}(2) = size(find(tunedLogical.si.onlyB_tuned ==1),2);
% tuned_count{1}(3) = size(find(tunedLogical.si.AandB_tuned ==1),2);
% tuned_count{1}(4) = size(find(tunedLogical.si.neither ==1),2); 
% 
% %tuning specificity
% tuned_count{2}(1) = size(find(tunedLogical.ts.onlyA_tuned ==1),2);
% tuned_count{2}(2) = size(find(tunedLogical.ts.onlyB_tuned ==1),2);
% tuned_count{2}(3) = size(find(tunedLogical.ts.AandB_tuned ==1),2);
% tuned_count{2}(4) = size(find(tunedLogical.ts.neither ==1),2); 

%% Calculate fraction tuned -only filtered data
%si
%fracTuned_si = tuned_count{1}/sum(tuned_count{1});
%filtered by min. 5 event and PF presence
for ss=sessionSelect
    fracTuned_si_filt(ss,:) = tuned_count_filt{ss}{1}/sum(tuned_count_filt{ss}{1});
    %ts
    %fracTuned_ts = tuned_count{2}/sum(tuned_count{2});
    %filtered by min. 5 event and PF presence
    fracTuned_ts_filt(ss,:) = tuned_count_filt{ss}{2}/sum(tuned_count_filt{ss}{2});
end


%% Plot stacked histograms to check distribution across sessions

figure
subplot(1,2,1)
hold on
title('S.I. tuned neurons')
bar(fracTuned_si_filt,'stacked')
subplot(1,2,2)
hold on
title('T.S. tuned neurons')
bar(fracTuned_ts_filt,'stacked')

%% Plot pie chart for each type of tuning criterion (no filter)

% figure('Position',[2050 520 1140 450])
% subplot(1,2,1)
% p = pie(fracTuned_si,{['A ', num2str(round(100*fracTuned_si(1))), '%'],...
%                         ['B ', num2str(round(100*fracTuned_si(2))), '%'],...
%                         ['A&B ', num2str(round(100*fracTuned_si(3))), '%'],...
%                         ['     Neither ', num2str(round(100*fracTuned_si(4))), '%']});
% hold on
% title('Percentage of active neurons tuned (S.I.)');
% colormap([0 0 1; 1 0 0; 1 0 1; 0.5 0.5 0.5])
% 
% subplot(1,2,2)
% p = pie(fracTuned_ts,{['A ', num2str(round(100*fracTuned_ts(1))), '%'],...
%                         ['B ', num2str(round(100*fracTuned_ts(2))), '%'],...
%                         ['A&B ', num2str(round(100*fracTuned_ts(3))), '%'],...
%                         ['     Neither ', num2str(round(100*fracTuned_ts(4))), '%']});
% hold on
% title('Percentage of active neurons tuned (T.S.)');
% colormap([0 0 1; 1 0 0; 1 0 1; 0.5 0.5 0.5])

%% Plot pie chart for each type of tuning criterion (filtered)

for ss=sessionSelect
    figure('Position',[2050 520 1140 450])
    subplot(1,2,1)
    p = pie(fracTuned_si_filt(ss,:),{['A ', num2str(round(100*fracTuned_si_filt(ss,1))), '%'],...
        ['B ', num2str(round(100*fracTuned_si_filt(ss,2))), '%'],...
        ['A&B ', num2str(round(100*fracTuned_si_filt(ss,3))), '%'],...
        ['     Neither ', num2str(round(100*fracTuned_si_filt(ss,4))), '%']});
    hold on
    title('Percentage of active neurons tuned (S.I.) - Filtered');
    colormap([0 0 1; 1 0 0; 1 0 1; 0.5 0.5 0.5])
    
    subplot(1,2,2)
    p = pie(fracTuned_ts_filt(ss,:),{['A ', num2str(round(100*fracTuned_ts_filt(ss,1))), '%'],...
        ['B ', num2str(round(100*fracTuned_ts_filt(ss,2))), '%'],...
        ['A&B ', num2str(round(100*fracTuned_ts_filt(ss,3))), '%'],...
        ['     Neither ', num2str(round(100*fracTuned_ts_filt(ss,4))), '%']});
    hold on
    title('Percentage of active neurons tuned (T.S.) - Filtered');
    colormap([0 0 1; 1 0 0; 1 0 1; 0.5 0.5 0.5])
end

%% Place into struct for export (only additionally filtered neurons)

%absolute counts for each subcategory and fractional counts
%tuned_fractions.tuned_count = tuned_count;
tuned_fractions.tuned_count_filt = tuned_count_filt;
%fraction tuned by SI criterion
%tuned_fractions.fracTuned_si = fracTuned_si;
tuned_fractions.fracTuned_si_filt = fracTuned_si_filt;
%fraction tuned by TS criterion
%tuned_fractions.fracTuned_ts = fracTuned_ts;
tuned_fractions.fracTuned_ts_filt = fracTuned_ts_filt;

% logicals for each sessions with A/B/AB/neigther tuned neurons (separate
% struct)
tuned_logicals.tuned_log_filt_si = tuned_log_filt_si;
tuned_logicals.tuned_log_filt_ts = tuned_log_filt_ts; 

end

