function [outputArg1,outputArg2] = PV_TC_correlation_single_ses(session_vars,tunedLogical,options)



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

%% Extract mean STC map in each spatial bin (not normalized and not occupancy divided) (100 bins)
%for each session
%correct only
for ss =1:size(session_vars,2)
    A_df{ss} = session_vars{ss}.Place_cell{1}.Spatial_tuning_curve;
    B_df{ss} = session_vars{ss}.Place_cell{2}.Spatial_tuning_curve;
    
    %Gs smoothed, but not normalized (nn) to itself
    A_STC_nn{ss} = session_vars{ss}.Place_cell{1}.Spatial_Info.rate_map_smooth{8};
    B_STC_nn{ss} = session_vars{ss}.Place_cell{2}.Spatial_Info.rate_map_smooth{8};
    
    A_df_both{ss} = A_df{ss}(:,AandB_tuned{ss});
    B_df_both{ss} = B_df{ss}(:,AandB_tuned{ss});
    
    A_df_onlyA{ss} = A_df{ss}(:,onlyA_tuned{ss});
    B_df_onlyA{ss} = B_df{ss}(:,onlyA_tuned{ss});
end


%% Normalize eahc STC ROI across both trials in non-norm STCs
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
%all neurons
PVcorr = corr(A_STC_nn{1}',B_STC_nn{1}');


%TC
TCcorr.all = corr(A_STC_nn{1},B_STC_nn{1});
TCcorr.Aonly = corr(A_STC_nn{1}(:,onlyA_tuned{1}),B_STC_nn{1}(:,onlyA_tuned{1}));
TCcorr.Bonly = corr(A_STC_nn{1}(:,onlyB_tuned{1}),B_STC_nn{1}(:,onlyB_tuned{1}));
TCcorr.AB = corr(A_STC_nn{1}(:,AandB_tuned{1}),B_STC_nn{1}(:,AandB_tuned{1}));

nanmean(diag(TCcorr.Aonly))
nanmean(diag(TCcorr.Bonly))
nanmean(diag(TCcorr.AB))

%sort TC correlation by ROI
diagTC = diag(TCcorr);

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
imagesc(TCcorr(I_sort_tc,I_sort_tc))
hold on
title('TC correlation - all neurons')

figure;
%check that sort worked as expected (same result)
plot(diag(TCcorr(I_sort_tc,I_sort_tc)))
plot(diag(TCcorr_sort),'r')

end
