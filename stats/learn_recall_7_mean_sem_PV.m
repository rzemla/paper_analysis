function [mean_lr7,sem_lr7] = learn_recall_7_mean_sem_PV(s1, s2)

%extract learn data from day
l7 = s1(7,~isnan(s1(7,:)));

%extract recall data from day 7
r7 = s2(7,~isnan(s2(7,:)));

%mean and sem
mean_lr7 = [mean(l7),mean(r7)];
sem_lr7 = [std(l7,0,2)./sqrt(numel(l7)),std(r7,0,2)./sqrt(numel(r7))];

end

