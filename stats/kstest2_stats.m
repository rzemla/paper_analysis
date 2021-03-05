function [stats_out] = kstest2_stats(x,y)

%run the KS test
[~,p_ks2,ks2stat] = kstest2(x,y);

%check if any nans in data
if sum([sum(isnan(x)), sum(isnan(y))]) ~=0
   disp('Nans in data set - write code to adjust (preprocess or change fxn)'); 
end

%total number of neurons for each comp.
n_samples = [numel(x),numel(y)];
%format the stats data for output
stats_out = [p_ks2,ks2stat, n_samples, nan];

end

