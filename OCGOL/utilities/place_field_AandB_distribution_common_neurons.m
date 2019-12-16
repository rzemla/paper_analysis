function [bin_conv_diff,pf_distance_metric] = place_field_AandB_distribution_common_neurons(session_vars,ROI_idx_tuning_class,remapping_corr_idx,cent_diff,pf_count_filtered,max_transient_peak)


%% Get common remapping neurons for distance calculation

common_idx = remapping_corr_idx.final.common;

%% TS and SI tuned neurons in each class (min 5 events and single fields)

%get indexes of single field TS A&B tuned neurons
AB_tuned_singleField_idx = intersect(common_idx,find(sum((pf_count_filtered == 1),1) == 2));

%extract those with single field
rad_angle_diff = cent_diff.angle_diff(AB_tuned_singleField_idx);

%convert angular difference to bin difference
bin_conv_diff = (rad_angle_diff./(2*pi))*100;

%max field idx for ts as
%A
max_field_idx.A = max_transient_peak{1}{1}(AB_tuned_singleField_idx);
max_field_idx.B = max_transient_peak{1}{2}(AB_tuned_singleField_idx);

%A trials
for rr=1:size(max_field_idx.A,2)
    %field A width
    pf_width.AB_single.comb(1,rr) = session_vars{1}.Place_cell{1}.placeField.width{AB_tuned_singleField_idx(rr)}(max_field_idx.A(rr));
    %field B width
    pf_width.AB_single.comb(2,rr) = session_vars{1}.Place_cell{2}.placeField.width{AB_tuned_singleField_idx(rr)}(max_field_idx.B(rr));
end

%distance/sum of place field metric
pf_distance_metric = bin_conv_diff./sum(pf_width.AB_single.comb,1);

%% Extract difference for common neurons

%A trials
for rr=1:size(max_field_idx.A,2)
    %field A width
    pf_width.common(1,rr) = session_vars{1}.Place_cell{1}.placeField.width{AB_tuned_singleField_idx(rr)}(max_field_idx.A(rr));
    %field B width
    pf_width.common(2,rr) = session_vars{1}.Place_cell{2}.placeField.width{AB_tuned_singleField_idx(rr)}(max_field_idx.B(rr));
end

%% Get 95% quantile for 
bin_dist_cutoff = quantile(bin_conv_diff,0.95);

%% Plot histogram of distributions

%plot
figure
subplot(2,1,1)
hold on
title('Place field distance') 
histogram(bin_conv_diff,20)

subplot(2,1,2)
hold on
title('Distance metric (D/(sum of place widths))') 
histogram(pf_distance_metric,20)

end

