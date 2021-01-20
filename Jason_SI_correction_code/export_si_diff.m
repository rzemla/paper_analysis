function [sum_row] = export_si_diff(orig_si_idx,out_temp, place_struct_input)

%% Make table for export for each animal and session
%index (file name), Original SI #, Corrected SI #, # difference, % difference, %
%percentage of original neurons present in corrected count

%# of original ROIs, number of corrected ROIs
si_tuned_counts = [numel(orig_si_idx), numel(out_temp.Info_Rate.cvc)];
%percentage of neurons detected 
frac_change = 100*(si_tuned_counts(2) - si_tuned_counts(1))./si_tuned_counts(1);
%percent of original neurons present in expanded set
intersect_with_orig = numel(setdiff(orig_si_idx', out_temp.Info_Rate.cvo'));
frac_present = (si_tuned_counts(1) - intersect_with_orig)./si_tuned_counts(1);

%TS tuned neurons
nb_ts_tuned = numel(find(place_struct_input.Tuning_Specificity.significant_ROI ==1));

%total ROIs
nbROI = numel(place_struct_input.Tuned_ROI_mask);

%composite row of params
sum_row = [si_tuned_counts, si_tuned_counts(2)-si_tuned_counts(1),frac_change,nb_ts_tuned,nbROI];

%add percentage of each tuned ROI as fraction of total ROIs
sum_row = [sum_row, (sum_row([1,2,5])./sum_row(6))*100];

%round to 2nd sig digit
sum_row = round(sum_row,2);

end

