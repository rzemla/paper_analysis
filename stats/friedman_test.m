function [friedman_stats] = friedman_test(x)
%runs Friedman test and outputs relevant stats

[p,tbl,stats] = friedman(x);
friedman_stats = cell2mat([p, tbl(2,5), stats.n, tbl(2,3)]);

end

