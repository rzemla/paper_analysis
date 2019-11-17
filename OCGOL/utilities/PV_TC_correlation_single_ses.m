function [correlation] = PV_TC_correlation_single_ses(session_vars,tunedLogical,task_selective_ROIs,options)


%TO-DO:
%1)A,B,A&B,Neither - give inputs that are screening for number of
%events
%2) generate A all and B all from Aonly. Bonly, and A&B inputs




%% Define tuned cell combinations across trials

%make conditional here for si or ts tuned neurons
switch options.tuning_criterion
    case 'si' %spatial information
        %for each session
        for ss =1:size(session_vars,2)
            %spatial information criterion
            Atuned{ss} = tunedLogical(ss).si.Atuned;
            Btuned{ss} = tunedLogical(ss).si.Btuned;
            
            AandB_tuned{ss} =  tunedLogical(ss).si.AandB_tuned;
            AorB_tuned{ss} = tunedLogical(ss).si.AorB_tuned;
            onlyA_tuned{ss} = tunedLogical(ss).si.onlyA_tuned;
            onlyB_tuned{ss} = tunedLogical(ss).si.onlyB_tuned;
            
            neither_tuned{ss} = tunedLogical(ss).si.neither;
        end
   case 'ts' %spatial information 
        for ss =1:size(session_vars,2)
            %spatial information criterion
            Atuned{ss} = tunedLogical(ss).ts.Atuned;
            Btuned{ss} = tunedLogical(ss).ts.Btuned;
            
            AandB_tuned{ss} =  tunedLogical(ss).ts.AandB_tuned;
            AorB_tuned{ss} = tunedLogical(ss).ts.AorB_tuned;
            onlyA_tuned{ss} = tunedLogical(ss).ts.onlyA_tuned;
            onlyB_tuned{ss} = tunedLogical(ss).ts.onlyB_tuned;
        end
        
end

%% Task-selective neuron idx
%add task-selective neurons (additional filters)
%neurons idxs associated with selective filtering for
%task-selectivity
select_filt_ROIs.A = task_selective_ROIs.A.idx;
select_filt_ROIs.B = task_selective_ROIs.B.idx;

%% Extract mean STC map in each spatial bin (not normalized and not occupancy divided) (100 bins)
%for each session
%correct only
for ss =1:size(session_vars,2)
    A_STC{ss} = session_vars{ss}.Place_cell{1}.Spatial_tuning_curve;
    B_STC{ss} = session_vars{ss}.Place_cell{2}.Spatial_tuning_curve;
    
    %Gs smoothed, but not normalized (nn) to itself
    A_STC_nn{ss} = session_vars{ss}.Place_cell{1}.Spatial_Info.rate_map_smooth{8};
    B_STC_nn{ss} = session_vars{ss}.Place_cell{2}.Spatial_Info.rate_map_smooth{8};
    
    A_STC_both{ss} = A_STC{ss}(:,AandB_tuned{ss});
    B_STC_both{ss} = B_STC{ss}(:,AandB_tuned{ss});
    
    A_STC_onlyA{ss} = A_STC{ss}(:,onlyA_tuned{ss});
    B_STC_onlyA{ss} = B_STC{ss}(:,onlyA_tuned{ss});
end


%% Normalize each STC ROI across both trials in non-norm STCs
for ss =1:size(session_vars,2)
        %get max value for each ROIs between trials
        max_STC_across_trials{ss} = max([A_STC_nn{ss};B_STC_nn{ss}]);
        %normalize each matrix to these values (tn = trial normalized)
        A_STC_tn{ss} = A_STC_nn{ss}./max_STC_across_trials{ss};
        B_STC_tn{ss} = B_STC_nn{ss}./max_STC_across_trials{ss};
    %for max value to normalize to by 1
    %in future, do normalization range (0,1)
end


%% Calculate PV and TC correlation matrixes for tuned subcategories
%correlations are done on non normalized STCs
%PV correlation - all neurons
PVcorr = corr(A_STC_nn{1}',B_STC_nn{1}', 'Rows', 'complete');


%TC correlation
TCcorr.all = corr(A_STC_nn{1},B_STC_nn{1});
TCcorr.Aonly = corr(A_STC_nn{1}(:,onlyA_tuned{1}),B_STC_nn{1}(:,onlyA_tuned{1}), 'Rows', 'complete');
TCcorr.Bonly = corr(A_STC_nn{1}(:,onlyB_tuned{1}),B_STC_nn{1}(:,onlyB_tuned{1}), 'Rows', 'complete');
TCcorr.AB = corr(A_STC_nn{1}(:,AandB_tuned{1}),B_STC_nn{1}(:,AandB_tuned{1}), 'Rows', 'complete');

%TC correlations for A-selective and B-selective filtered 
TCcorr.Aselective =  corr(A_STC_nn{1}(:,select_filt_ROIs.A),B_STC_nn{1}(:,select_filt_ROIs.A), 'Rows', 'complete');
TCcorr.Bselective =  corr(A_STC_nn{1}(:,select_filt_ROIs.B),B_STC_nn{1}(:,select_filt_ROIs.B), 'Rows', 'complete');

% nanmean(diag(TCcorr.Aonly))
% nanmean(diag(TCcorr.Bonly))
% nanmean(diag(TCcorr.AB))

%sort TC correlation by ROI
diagTC = diag(TCcorr.all);

[sort_tc_val, I_sort_tc] = sort(diagTC,'ascend');

%TC correlation following reordering by Pearson
TCcorr_sort = corr(A_STC_nn{1}(:,I_sort_tc),B_STC_nn{1}(:,I_sort_tc));

figure; 
subplot(2,1,1)
imagesc(PVcorr)
hold on
title('PV correlation - all neurons')
hold off
subplot(2,1,2)
imagesc(TCcorr.all(I_sort_tc,I_sort_tc))
hold on
title('TC correlation - all neurons')

figure;
%check that sort worked as expected (same result)
plot(diag(TCcorr.all(I_sort_tc,I_sort_tc)),'b')
plot(diag(TCcorr_sort),'r')

%% Export correlation struct

%PV correlation
correlation.PVcorr = PVcorr;

%TC correlation
correlation.TCcorr = TCcorr;

end

