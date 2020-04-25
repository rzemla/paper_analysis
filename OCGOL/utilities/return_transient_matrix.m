function [transient_mat_all] = return_transient_matrix(event_onset_all,match_list,lap_nb,ses)

%make blank 2D matrix for each each assigned and not assigned neuron for
%give session
transient_mat_all = nan(size(match_list,1),size(lap_nb,1));
%fill indexes
neuron_occup_log = ~isnan(match_list(:,ses));
%ROIs corresponding to the matching indices
neuron_idx_ordered = match_list(neuron_occup_log,ses);
%assign these event maps to the matching transient matrix
transient_mat_all(neuron_occup_log,:) = event_onset_all(neuron_idx_ordered,:);

end

