function [sortOrder] = sortSTC(STC,binRange)
%input: STC concatenated STCs from each set of trials = ROI x binRange to sort

[~,maxBin] = max(STC(:,binRange)', [], 1);
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,sortOrder] = sort(maxBin,'ascend');

end

