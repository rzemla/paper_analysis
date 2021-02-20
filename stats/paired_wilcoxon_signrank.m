function [stats] = paired_wilcoxon_signrank(x,y)
%no nan data inputs - modify later
%run paired Wilcoxon signed rank and return p val, signed ranks

%run test, get p-value
[p,~,~] = signrank(x,y);
%use modified function to get the total rank to report in paper
[ranks] = signrank_all_ranks(x,y);
%number of samples
n = numel(x);

%output the variable
%p val, sum of ranks, nb samples, DOF
stats = [p, ranks.total,n,n-1];

end

