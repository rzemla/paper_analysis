function [ts_bin_conv_diff,pf_distance_metric_ts] = place_field_AandB_distribution(session_vars,ROI_idx_tuning_class,remapping_corr_idx,cent_diff,pf_count_filtered,max_transient_peak)



%% TS and SI tuned neurons in each class (min 5 events and single fields)

%get TS tuned neurons(event filtered) for each class
%ts
ts_AB_tuned = ROI_idx_tuning_class.ts.AB;

%si
si_AB_tuned = ROI_idx_tuning_class.si.AB;

%get indexes of single field TS A&B tuned neurons
ts_AB_tuned_singleField_idx = intersect(ts_AB_tuned,find(sum((pf_count_filtered == 1),1) == 2));
%si
si_AB_tuned_singleField_idx = intersect(si_AB_tuned,find(sum((pf_count_filtered == 1),1) == 2));

%extract those with single field
ts_rad_angle_diff = cent_diff.angle_diff(ts_AB_tuned_singleField_idx);

%extract those with single field
si_rad_angle_diff = cent_diff.angle_diff(si_AB_tuned_singleField_idx);

%convert angular difference to bin difference
ts_bin_conv_diff = (ts_rad_angle_diff./(2*pi))*100;
%si
si_bin_conv_diff = (si_rad_angle_diff./(2*pi))*100;

%max field idx for ts as
%A
ts_max_field_idx.A = max_transient_peak{1}{1}(ts_AB_tuned_singleField_idx);
ts_max_field_idx.B = max_transient_peak{1}{2}(ts_AB_tuned_singleField_idx);

%A trials
for rr=1:size(ts_max_field_idx.A,2)
    %field A width
    pf_width.ts.AB_single.comb(1,rr) = session_vars{1}.Place_cell{1}.placeField.width{ts_AB_tuned_singleField_idx(rr)}(ts_max_field_idx.A(rr));
    %field B width
    pf_width.ts.AB_single.comb(2,rr) = session_vars{1}.Place_cell{2}.placeField.width{ts_AB_tuned_singleField_idx(rr)}(ts_max_field_idx.B(rr));
end

%distance/sum of place field metric
pf_distance_metric_ts = ts_bin_conv_diff./sum(pf_width.ts.AB_single.comb,1);

%get place fields

%plot
figure
subplot(2,1,1)
hold on
title('Place field distance') 
histogram(ts_bin_conv_diff,20)

subplot(2,1,2)
hold on
title('Distance metric (D/(sum of place widths))') 
histogram(pf_distance_metric_ts,20)

end

