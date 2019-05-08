%place field tester script

%1) 
%Zaremba 2017
% For all cells, rate maps were formed by dividing the number of transients
% initiated in each spatial bin by the occupancy of that bin.
% We calculated rate maps with 100 position bins and smoothed 
% with a Gaussian kernel (? = 3 bins). To define place fields for 
% cells that were identified as containing significant spatial information,
% we fit each local maximum in the rate map with a Gaussian, merged 
% overlapping putative fields and then discarded any with an area less than 
% 50% of the largest.