function [p_corr_orig_order] = holm_sidak_p_adj(p,c,alpha)
%return Holm-Sidak adjusted p values

%sort p values in ascending order (smallest to biggest)
[~,I]=sort(p); %sorted p-value vector
p_sort=p(I);

%list of correction p-values
J=1:c; %How many corrected alpha value?

%Sidak alpha corrected values
alphacorr=1-((1-alpha).^(1./(c-J+1)));

%p value adjustments here (Holm-Sidak - prism correction - same result)
%https://personal.utdallas.edu/~herve/abdi-Holm2010-pretty.pdf

%find at which test fail to reject null hypothesis
test_p = p_sort<alphacorr;
%test_p = [true true false false] ;

%index at which test fails
idx_test_fail = find(test_p == 0,1,'first');

%adjust all p-values
p_corrected = 1-(1-p_sort).^(c-J+1);

%crop at test which failed - i.e. don't keep adjusting p-values but copy
%the previous adjusted p-value at the test which failed
%check if test found an early faiure to reject
if ~isempty(idx_test_fail)
    p_corrected(idx_test_fail:end) = p_corrected(idx_test_fail);    
end

%rearrange corrected p-values in order of input
p_corr_orig_order = p_corrected(I);

end

