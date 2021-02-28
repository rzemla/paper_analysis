function [t_stats] = paired_ttest(x,y)
%run paired t-test and output formatted stats data

[~,p,~,stats] = ttest(x,y);
t_stats = [p, stats.tstat,stats.df+1, stats.df];

end

