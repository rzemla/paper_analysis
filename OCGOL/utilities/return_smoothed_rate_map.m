function [rate_map_ev_sm_oc_sm_clipped] = return_smoothed_rate_map(event_map,occupancy_s)
%DEFINE INPUTS
%event_map: ROI x 100 bin matrix
%occup_map: row vector - 100 bins seconds


%% Define gaussian smoothing filter

%sigma =3 for Gaussian filter
options.sigma_filter = 3;
gaussFilter = define_Gaussian_kernel(options);

%% Extend the event and occupancy maps, Gaussian smooth, and clip the ends

%turn into function and return smoothed map
%extend map for circular smoothing
extended_event_map = [event_map(:,51:100),event_map,event_map(:,1:50)];
%extend occupnacy map for cicular smoothing
extended_occup_map = [occupancy_s(51:100),occupancy_s,occupancy_s(1:50)];

%circularly smooth event map with Gaussian
for rr=1:size(extended_event_map,1)
    extend_smooth_event_map(rr,:) = conv(extended_event_map(rr,:),gaussFilter, 'same');
end

%circularly smooth occupancy map with Gaussian
extend_smooth_occup_map = conv(extended_occup_map,gaussFilter, 'same');

%divide smoothed event map by occupancy map
rate_map_ev_sm_oc_sm = extend_smooth_event_map./extend_smooth_occup_map;

%clip the ends and export
rate_map_ev_sm_oc_sm_clipped = rate_map_ev_sm_oc_sm(:,51:150)';


end

