function [outputArg1,outputArg2] = centroid_dist(tunedLogicals, max_bin_rate,options)
%distribution of max centroids across length of track (bin space)



%% Generate count of at each bin (1-100)

%for each trial
for tt=1:2
    [centroid_count(tt,:),edges] = histcounts(max_bin_rate{1}{tt}, 0.5:1:100.5);
end




%% Plot histogram of max rate for each trial
% for all neurons
figure;
subplot(2,1,1)
hold on
title('A trials')
histogram('BinEdges', edges,'BinCounts', centroid_count(1,:))
subplot(2,1,2)
hold on
title('B trials')
histogram('BinEdges', edges,'BinCounts', centroid_count(2,:))

end

