function [ttest_stats] = unpaired_ttest(x,y)

%unpaired t-test code
[~,p,~,stats] = ttest2(x,y);
n_sample_each = [sum(~isnan(x)),sum(~isnan(y)) ];
ttest_stats = [p,stats.tstat, n_sample_each,stats.df];

end

