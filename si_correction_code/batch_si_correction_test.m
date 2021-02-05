%% Make list of directories for import


%% notes
%tuning specificity values are not filtered for minimum events
%overall - those ROIs that are both SI and TS tuned with and filtered for
%min events - not just SI adjusted ROIs

%Spatial info ROIs
%out_temp.Info_Content.cvo
%measuring the spatial info rate in paper
%out_temp.Info_Rate.

%% Extract corrected vs. original ROIs + stats
%original SI tuned indices
%orig_ts_idx = find(Place_cell{1, 1}.Tuning_Specificity.significant_ROI ==1);

%input place struct; 1 - A corr, 2 - B corr, 4 - all A, 5 - all B
place_struct = Place_cell{1,1} 

%
orig_si_idx = find(Place_cell{1, 1}.Spatial_Info.significant_ROI ==1);
%corrected indices of tuned ROIs
out_temp = whichSignificant2(place_struct);

[sum_row] = export_si_diff(orig_si_idx,out_temp)

%% Compare the ROIs

%these should be the same 
isequal(orig_si_idx, out_temp.Info_Rate.cvo)
%original ROIs: out_temp.Info_Rate.cvo
%corrected ROIs: out_temp.Info_Rate.cvc

%consistent with additional ROIs being adding to the exisiting ones -
%metric being less stringent
isequal(orig_si_idx,intersect(out_temp.Info_Rate.cvc, out_temp.Info_Rate.cvo))






