function [trialRanges] = discoverTrialType(trialTypeCh ,position)

%% Texture signal discovery

%how many voltage bins; 20 = 0.25 V, 10 = 0.50 V
voltageBins = 10;

%discover voltage peaks
[pks, pks_idx] = findpeaks(trialTypeCh, 'MinPeakHeight',0.4);

%exclude voltage peaks based on positional proximity
peak_diff_pos = [position(pks_idx), [nan ;diff(position(pks_idx))]];

%find voltage signals with positional diff < 0.3 cm
exclude_peaks_idx = find(peak_diff_pos(:,2) < 0.3 & peak_diff_pos(:,2) >= 0.0);
keep_peaks = false(size(peak_diff_pos,1),1);
keep_peaks(exclude_peaks_idx) = 1;
keep_peaks = ~keep_peaks;

%final peaks to include
pks_final_idx = pks_idx(keep_peaks);
pks_final = pks(keep_peaks);

%bin the voltage signals associated with textures
[N_voltage_bins,~] = histcounts(pks_final,voltageBins);

%find unique bins (= # of textures for input into k means)
uniqueTexBins = length(find(N_voltage_bins ~= 0));

%k-means cluster - find each texture and voltage centroid
[idx,C] = kmeans(pks_final,uniqueTexBins);

%sort k-means cluster centroids in ascending voltage
C = sort(C);

%assign voltage range
voltRange = 0.05;

%assign discovered texture ranges
%preallocate
trialRanges = zeros(size(C,1),2);

for ii = 1:size(C,1)
    %low voltage thres
    trialRanges(ii,1) = C(ii) - voltRange;
    %high voltage thres
    trialRanges(ii,2) = C(ii) + voltRange;
end

%%info output

%# of trial types
nb_trial_type = size(trialRanges,1);

%print number of trial types found
str_output = sprintf('Number of trial types discovered: %d',nb_trial_type);
disp(str_output);



end

